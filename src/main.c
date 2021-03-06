#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

#include "newick_tree.h"

#include "lca.h"
#include "traits.h"
#include "getopt.h"
#include "hash_table.h"
#include "monitored_memory.h"

//#define NDEBUG_DEVEL

real_t get_distance_from_root(newick_node *n) {
    real_t dist = n->dist;
    if (n->parent != NULL) dist += get_distance_from_root(n->parent);
    return dist;
}

mpfr_ptr get_distance_from_root_mpfr(newick_node *n) {
    mpfr_ptr dist = (mpfr_ptr) malloc(sizeof(mpfr_t));
    mpfr_init2(dist, MPFR_PRECISION);
    mpfr_set_zero(dist, 0);
    
    mpfr_set_d(dist, n->dist, MPFR_RNDN);
    if (n->parent != NULL) {
        mpfr_ptr parent_dist = get_distance_from_root_mpfr(n->parent);
        mpfr_add(dist, dist, parent_dist, MPFR_RNDN);
        mpfr_clear(parent_dist);
        free(parent_dist);
    }
    return dist;
}

real_t max_dist_from_root(newick_node *t) {
    newick_child *p;
    real_t dist = 0.f;
    for (p = t->child; p != NULL; p = p->next) {
        dist = fmax(dist, max_dist_from_root(p->node));
    }
    return dist + t->dist;
}

mpfr_ptr max_dist_from_root_mpfr(newick_node *n) {
    newick_child *p;
    mpfr_ptr dist = (mpfr_ptr) malloc(sizeof(mpfr_t));
    mpfr_init2(dist, MPFR_PRECISION);
    mpfr_set_zero(dist, 0);
    
    for (p = n->child; p != NULL; p = p->next) {
        mpfr_ptr rec_dist = max_dist_from_root_mpfr(p->node);
        mpfr_max(dist, dist, rec_dist, MPFR_RNDN);
        mpfr_clear(rec_dist);
        free(rec_dist);
    }
    mpfr_add_d(dist, dist, n->dist, MPFR_RNDN);
    return dist;
}

int combinations(int n, int k) {
    assert(n >= 0);
    assert(k >= 0);
    if (k > n) return 0;
    if (k == 0 || n == k || n <= 1) return 1;
    else return combinations(n - 1, k - 1) + combinations(n - 1, k);
}

real_t H(int l_1) {
    int x;
    real_t res = 1.;
    for (x = 1; x < l_1; ++x) {
        res *= x * x * 0.5;
    }
    return res * l_1;
}

real_t theorem3_func(int i, int l_i, int n, int_array_array *lambda, real_array *spec_dists, node2real_array_array *m) {
    int j, k;
    real_t P = 0;
    printf("\ttheorem3 called with i=%d, l_i=%d, n=%d, m[%d]=%d\n", i, l_i,n,i,array_size(m->array[i]));
    for (j = 0; j < array_size(m->array[i]); ++j) {
        real_t nom = expf(-lambda->array[i].array[j] * (spec_dists->array[i - 1] - spec_dists->array[i]));
        real_t denom = 1.0;
        for (k = 0; k < array_size(m->array[i]); ++k) {
            if (k == j) continue;
            denom *= (lambda->array[i].array[k] - lambda->array[i].array[j]);
        }
        P += nom / denom;
        
        printf("\t\tnominator:%f\n",nom);
        printf("\t\tdenominator:%f\n",denom);
        
    }
    printf("\t\tth3 ret: %f\n",P);
    return P;
}

mpfr_ptr theorem3_func_mpfr(int i, int l_i, int n, int_array_array *lambda, mpfr_t_array *spec_dists, node2mpfr_t_array_array *m) {
    int j, k;
    mpfr_ptr P, nom, denom;
    P = (mpfr_ptr) malloc(sizeof(mpfr_t));
    nom = (mpfr_ptr) malloc(sizeof(mpfr_t));
    denom = (mpfr_ptr) malloc(sizeof(mpfr_t));
    mpfr_inits2(MPFR_PRECISION, P, nom, denom, (mpfr_ptr) 0);
    
    mpfr_set_zero(P, 0);
    for (j = 0; j < array_size(m->array[i]); ++j) {
        mpfr_set_d(nom, -lambda->array[i].array[j] * (spec_dists->array[i - 1] - spec_dists->array[i]), MPFR_RNDN);
        mpfr_exp(nom, nom, MPFR_RNDN);

        mpfr_set_d(denom, 1.0, MPFR_RNDN);
        for (k = 0; k < array_size(m->array[i]); ++k) {
            if (k == j) continue;
            mpfr_mul_d(denom, denom, (lambda->array[i].array[k] - lambda->array[i].array[j]), MPFR_RNDN);
        }
        mpfr_div(nom, nom, denom, MPFR_RNDN);
        mpfr_add(P, P, nom, MPFR_RNDN);
    }
    mpfr_clears(nom, denom, (mpfr_ptr) 0);
    free(nom);
    free(denom);
    
    return P;
}

real_t theorem2_func(int i, int l_i, int n, int_array *g, int_array_array *lambda, real_array *spec_dists, node2real_array_array *m) {
    
    if (i == n - 2 && l_i == n - 1) {
#ifndef NDEBUG_DEVEL
        puts("\tbase");
#endif
        return 1;
    }

#ifndef NDEBUG_DEVEL
    printf("theorem2 called:\n");
    printf("\ti:%d l_i:%d, n:%d\n", i, l_i, n);
#endif
    
    real_t res = 0;
    int l_i_plus_1;
    //printf("l_i_plus_1 inited with %f\n", fmax(l_i, g->array[i + 1]));
    for (l_i_plus_1 = fmax(l_i, g->array[i + 1]); l_i_plus_1 < n; ++l_i_plus_1) {
        res += theorem3_func(i + 1, l_i_plus_1, n, lambda, spec_dists, m) * theorem2_func(i + 1, l_i_plus_1, n, g, lambda, spec_dists, m);
    }
    printf("\tth2 ret: %f\n", res);
    return res;
}

mpfr_ptr theorem2_func_mpfr(int i, int l_i, int n, int_array *g, int_array_array *lambda, real_array *spec_dists, node2mpfr_t_array_array *m) {
    mpfr_ptr x = NULL;
    return x;
}

void do_get_speciation_distances(newick_node *t, real_t distance, real_array *speciation_distances, real_t max_dist_from_root) {
    newick_child *p;
    if (t->childNum > 0) // speciation happened
        append_real_array(speciation_distances, max_dist_from_root - (distance + t->dist));
    for (p = t->child; p != NULL; p = p->next) {
        do_get_speciation_distances(p->node, distance + t->dist, speciation_distances, max_dist_from_root);
    }
}

real_array *get_speciation_distances(newick_node *t, real_t max_dist_from_root) {
    real_array *speciation_distances = (real_array *)malloc(sizeof(real_array));
    init_real_array(speciation_distances);
    do_get_speciation_distances(t, 0.0f, speciation_distances, max_dist_from_root);
    append_real_array(speciation_distances, FLT_MAX);

    qsort(speciation_distances->array, speciation_distances->last - speciation_distances->array,
        sizeof(real_t), (int(*)(const void*,const void*))flt_cmp_desc);
    return speciation_distances;
}

int do_get_gene_lineages(newick_node *t, real_t limit, real_t distance, real_t max_dist_from_root) {
    int lineages = 0;
    newick_child *p;

#ifndef NDEBUG_DEVEL
    printf("\t%f <= %f? ", max_dist_from_root - distance - t->dist, limit);
#endif

//    if ((max_dist_from_root - distance - t->dist) <= limit) {
    if (real_cmp((max_dist_from_root - distance - t->dist), limit) <= 0) {
#ifndef NDEBUG_DEVEL
        puts("y");
#endif
        return lineages + 1;
    }
#ifndef NDEBUG_DEVEL
    else puts("n");
#endif

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist, max_dist_from_root);
    }

    return lineages;
}

int_array *get_gene_lineages(real_array *speciation_distances, newick_node *t, real_t max_dist_from_root) {
    real_t *p;
    int_array *lineages = (int_array *)malloc(sizeof(int_array));

    init_int_array(lineages);
    for (p = speciation_distances->array; p != speciation_distances->last; ++p) {
//        if (max_dist_from_root < *p) {
        if (real_cmp(max_dist_from_root, *p) < 0) {
            append_int_array(lineages, 0);
        } else {
            append_int_array(lineages, do_get_gene_lineages(t, *p, 0.f, max_dist_from_root));
        }
    }
    return lineages;
}

int do_get_gene_lineages_for_k(newick_node *t, real_t limit, real_t distance, real_t max_dist_from_root, int_array *species_turns, int_array *gene_turns, int turn) {
    int lineages = 0;
    newick_child *p;
    int *gp, *sp;

    if (turn != -1) append_int_array(gene_turns, turn);

    // check that gene_turns begins with_species_turns
    int same_topology = 1;
    for (gp = gene_turns->array, sp = species_turns->array; gp != gene_turns->last && sp != species_turns->last; ++gp, ++sp) {
        if (*gp != *sp) {
            same_topology = 0;
            break;
        }
    }

    if (!same_topology) return 0;

#ifndef NDEBUG_DEVEL
    printf("\t%f <= %f? ", max_dist_from_root - distance - t->dist, limit);
#endif

//    if ((max_dist_from_root - distance - t->dist) <= limit) {
    if (real_cmp((max_dist_from_root - distance - t->dist), limit) <= 0) {
#ifndef NDEBUG_DEVEL
        puts("yes");
#endif
        return lineages + 1;
    }
#ifndef NDEBUG_DEVEL
    else puts("no");
#endif

    int child_index = 0;
    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages_for_k(p->node, limit, distance + t->dist, max_dist_from_root, species_turns, gene_turns, child_index++);
    }

    // remove added turn
    if (gene_turns->last - gene_turns->array > 0) --(gene_turns->last);

    return lineages;
}

int get_gene_lineages_for_k(real_t *speciation_distance, newick_node *t, real_t max_dist_from_root, int_array *species_turns) {
    int lineages = 0;
    int_array *gene_turns = (int_array *) malloc(sizeof(int_array));
    init_int_array(gene_turns);
    lineages = do_get_gene_lineages_for_k(t, *speciation_distance, 0.f, max_dist_from_root, species_turns, gene_turns, -1);
    return lineages;
}

/*
 * distances from the root for each node
 */
void do_get_distance_array(newick_node *t, real_t distance, node2real_array *dist_array) {
    newick_child *p;
    node2real pair;
    pair.val = distance + t->dist;
    pair.node = t;
    append_node2real_array(dist_array, pair);
    for (p = t->child; p != NULL; p = p->next) {
        do_get_distance_array(p->node, distance + t->dist, dist_array);
    }
}

int get_tree_coalescence_count(newick_node *t) {
    int count = 0;
    newick_child *p;

    if (t == NULL || t->child == NULL) return 0;
    count += 1;

    for (p = t->child; p != NULL; p = p->next) {
        count += get_tree_coalescence_count(p->node);
    }

    return count;
}

/*
 * coalescence array here is 0-indexed (it's 1-indexed in the original paper)
 * nodes are sorted by distance from root
 */
node2real_array *get_coalescence_array(newick_node *t) {
    node2real *n2f;
    node2real_array *dist_array = (node2real_array *)malloc(sizeof(node2real_array));
    node2real_array *coalescence_array = (node2real_array *)malloc(sizeof(node2real_array));

    init_node2real_array(dist_array);
    init_node2real_array(coalescence_array);

    do_get_distance_array(t, 0.f, dist_array);
    qsort(dist_array->array, dist_array->last - dist_array->array,
        sizeof(node2real), (int(*)(const void*,const void*))node2real_cmp);

    //int c = 0;
    for (n2f = dist_array->array; n2f != dist_array->last; ++n2f) {
        if (n2f->node->childNum > 0) {
            node2real pair;
            pair.node = n2f->node;
            pair.val = n2f->val;

            append_node2real_array(coalescence_array, pair);
        }
    }

    return coalescence_array;
}

int do_get_topology_prefix(int_array *x, newick_node *n, newick_node *target) {
    if (n == target) return 1;

    int child_idx;
    newick_child *p;
    for (p = n->child, child_idx = 0; p != NULL; p = p->next, ++child_idx) {
        int found = do_get_topology_prefix(x, p->node, target);
        if (found) {
            append_int_array(x, child_idx);
            return 1;
        }
    }
    return 0;
}

int_array get_topology_prefix(newick_node *n, newick_node *target) {
    int_array topology_prefix;
    init_int_array(&topology_prefix);
    do_get_topology_prefix(&topology_prefix, n, target);
    return topology_prefix;
}

int do_get_exit_branches(newick_node *n, real_t limit, real_t distance, int_array *topology_prefix, int_array *self_topology_prefix, real_t max_dist_from_root) {
    int exit_branches_count = 0;
    newick_child *p;

#ifndef NDEBUG_DEVEL
    printf("\t%f <= %f? ", max_dist_from_root - distance - n->dist, limit);
#endif

//    if ((max_dist_from_root - distance - n->dist) <= limit) {
    if (real_cmp((max_dist_from_root - distance - n->dist), limit) <= 0) {

#ifndef NDEBUG_DEVEL
        puts("yes");
        printf("\tchecking for topology prefix...\n");
        printf("\t\tself topology prefix: ");
#endif
        int *ip, *jp;
#ifndef NDEBUG_DEVEL
        for (ip = self_topology_prefix->array; ip != self_topology_prefix->last; ++ip) {
            printf("%d", *ip);
        }
        printf("\n");
        printf("\t\tspec topology prefix: ");

        for (jp = topology_prefix->array; jp != topology_prefix->last; ++jp) {
            printf("%d", *jp);
        }   
#endif
        int is_subset = 1;
        for (ip = self_topology_prefix->array, jp = topology_prefix->array;
                ip != self_topology_prefix->last && jp != topology_prefix->last;
                ++ip, ++jp) {
            if (*ip != *jp) {
                is_subset = 0;
                break;
            }
        }
#ifndef NDEBUG_DEVEL
        if (is_subset) printf("\tsubset!\n"); else printf("\tnot subset!\n");
#endif
        if (is_subset) return exit_branches_count + 1;
    } else {
#ifndef NDEBUG_DEVEL
        puts("no");
#endif
    }

    int child_idx;
    for (p = n->child, child_idx = 0; p != NULL; p = p->next, ++child_idx) {
        append_int_array(self_topology_prefix, child_idx);
        exit_branches_count += do_get_exit_branches(p->node, limit, distance + n->dist, topology_prefix, self_topology_prefix, max_dist_from_root);
        // remove topology index once the sub-tree has been walked
        --(self_topology_prefix->last);
    }

    return exit_branches_count;
}

int get_exit_branches(newick_node *n, real_t distance, int_array *topology_prefix, real_t max_dist_from_root) {
    int_array self_topology_prefix;
    init_int_array(&self_topology_prefix);
    int exit_branches_count = do_get_exit_branches(n, distance, 0.f, topology_prefix, &self_topology_prefix, max_dist_from_root);
    free(self_topology_prefix.array);
    return exit_branches_count;
}



int do_get_exit_branches_fixed(newick_node *n, real_t limit, real_t distance, real_t max_dist_from_root, char_ptr_array *species_leaf_descendants, hash_table *gene_leaf_descendants) {
    int exit_branches_count = 0;
    newick_child *p;
    
#ifndef NDEBUG_DEVEL
    printf("\t%.16f <= %.16f? ", max_dist_from_root - distance - n->dist, limit);
#endif
    
    //    if ((max_dist_from_root - distance - n->dist) <= limit) {
    if (real_cmp((max_dist_from_root - distance - n->dist), limit) <= 0) {
        
#ifndef NDEBUG_DEVEL
        puts("yes");
        printf("\tchecking for descendants...\n");
        printf("\t\tgene node descendants: ");
#endif
        char_ptr *iip, *jjp;
        char_ptr_array *leaf_descendants = htab_lookup(gene_leaf_descendants, n, sizeof(n), compar_addr);

#ifndef NDEBUG_DEVEL
        for (iip = leaf_descendants->array; iip != leaf_descendants->last; ++iip) {
            printf("%s ", *iip);
        }
        printf("\n");
        printf("\t\tspecies descendants: ");
        
        for (iip = species_leaf_descendants->array; iip != species_leaf_descendants->last; ++iip) {
            printf("%s ", *iip);
        }
#endif
        int is_subset = 1;
        for (iip = leaf_descendants->array; iip != leaf_descendants->last; ++iip) {
            int found = 0;
            for (jjp = species_leaf_descendants->array; jjp != species_leaf_descendants->last; ++jjp) {
                if (strcmp(*iip, *jjp) == 0) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                is_subset = 0;
                break;
            }
        }
#ifndef NDEBUG_DEVEL
        if (is_subset) printf("\tsubset!\n"); else printf("\tnot subset!\n");
#endif
        if (is_subset) return exit_branches_count + 1;
    }
#ifndef NDEBUG_DEVEL
      else puts("no");
#endif
    
    int child_idx;
    for (p = n->child, child_idx = 0; p != NULL; p = p->next, ++child_idx) {
        exit_branches_count += do_get_exit_branches_fixed(p->node, limit, distance + n->dist, max_dist_from_root, species_leaf_descendants, gene_leaf_descendants);
    }
    
    return exit_branches_count;
}

int get_exit_branches_fixed(newick_node *n, real_t distance, real_t max_dist_from_root, char_ptr_array *species_leaf_descendants, hash_table *gene_leaf_descendants) {
    int_array self_topology_prefix;
    init_int_array(&self_topology_prefix);
    int exit_branches_count = do_get_exit_branches_fixed(n, distance, 0.f, max_dist_from_root, species_leaf_descendants, gene_leaf_descendants);
    free(self_topology_prefix.array);
    return exit_branches_count;
}

void add_nodes_in_interval(newick_node_ptr_array *arr, newick_node *n, real_t start_interval, real_t end_interval, real_t max_root_distance) {
    //TODO: optimize, can call "get_distance_from_root" only once (in driver, move this to do_driver)
    real_t root_dist = get_distance_from_root(n);
    real_t dist = max_root_distance - root_dist;

#ifndef NDEBUG_DEVEL
    printf("\t\tnode@%p with taxon %s: %f <= %f < %f", n, n->taxon, start_interval, dist, end_interval);
#endif

//    if (dist >= start_interval && dist < end_interval) {
    if (real_cmp(dist, start_interval) >= 0 && real_cmp(dist, end_interval) < 0) {

#ifndef NDEBUG_DEVEL
        printf(" y\n");
#endif
        append_newick_node_ptr_array(arr, n);
    }
#ifndef NDEBUG_DEVEL
      else printf(" n\n");
#endif

    newick_child *p;
    for (p = n->child; p != NULL; p = p->next) {
        add_nodes_in_interval(arr, p->node, start_interval, end_interval, max_root_distance);
    }
}

void do_bead_tree(newick_node *t, real_t distance, real_array *speciation_distances, real_t max_dist_from_root) {
    real_t *fp;

    newick_node *bead;
    newick_child *p, *q, *child, *child_head, *new_child;
    
    // if a node is a leaf, we need to extend it so it's distance from root = max_dist_from_root
    real_t root_dist = distance + t->dist;
    assert(real_cmp(root_dist, 0.0) >= 0);
//    if (t->childNum == 0 && root_dist < max_dist_from_root) {
    if (t->childNum == 0 && real_cmp(root_dist, max_dist_from_root) < 0) {
        t->dist += max_dist_from_root - root_dist;
    }

    real_t start_interval, end_interval;

    end_interval = max_dist_from_root - (distance + t->dist);
    for (p = t->child; p != NULL; ) {
        do_bead_tree(p->node, distance + t->dist, speciation_distances, max_dist_from_root);

        // find if there are any speciation times in between this node and it's child
        start_interval = max_dist_from_root - (distance + t->dist + p->node->dist);

        child = p;
        child_head = t->child;

        p = p->next;

        for (fp = (speciation_distances->last - 1); fp != (speciation_distances->array - 1); --fp) {
//            if (*fp > start_interval && *fp < end_interval) {
            if (real_cmp(*fp, start_interval) > 0 && real_cmp(*fp, end_interval) < 0) {

#ifndef NDEBUG_DEVEL
                printf("\t%f in interval <%f, %f>\n", *fp, start_interval, end_interval);
#endif
                bead = (newick_node *) monitored_malloc(sizeof(newick_node));
                bead->child = child;
                bead->parent = child->node->parent;
                child->node->parent = bead;
                bead->childNum = 1;
                bead->taxon = "bead";

                bead->dist = end_interval - *fp;
                child->node->dist = *fp - start_interval; // should've been end_interval - start_interval

                // now attach bead as one of t's children:
                new_child = (newick_child *) monitored_malloc(sizeof(newick_child));
                new_child->node = bead;
                new_child->next = child->next;

                // we need to insert a bead right into the place where the old node was to keep the same topology
                // (treat the previous node correctly)
                if (child_head == child) {
                    child_head = new_child;
                } else {
                    for (q = child_head; q != NULL; q = q->next) {
                        if (q->next == child) break;
                    }
                    assert(q);
                    q->next = new_child;
                }

                bead->child->next = NULL;

                t->child = child_head;
                child = new_child;
                start_interval = *fp;
            }
        }
    }
}

void bead_tree(newick_node *n, real_array *speciation_times, real_t max_dist_from_root) {
    do_bead_tree(n, 0.f, speciation_times, max_dist_from_root);
}

/*
 * nodes array
 * nodes are sorted by distance
 */
node2int_array *get_indexed_array(newick_node *t) {
    node2real *n2f;
    node2real_array *dist_array = (node2real_array *)malloc(sizeof(node2real_array));
    node2int_array *nodes_array = (node2int_array *)malloc(sizeof(node2int_array));

    init_node2real_array(dist_array);
    init_node2int_array(nodes_array);

    do_get_distance_array(t, 0.f, dist_array);
    qsort(dist_array->array, dist_array->last - dist_array->array,
        sizeof(node2real), (int(*)(const void*,const void*))node2real_cmp);

    int c = 0;
    for (n2f = dist_array->array; n2f != dist_array->last; ++n2f) {
        node2int pair;
        n2f->node->id = c;
        pair.node = n2f->node;
        pair.val = c++;

        append_node2int_array(nodes_array, pair);
    }

    return nodes_array;
}

void construct_y_matrix(hash_table *mat_idx_tab, newick_node_ptr_array_array *Y, real_array *spec_dists, newick_node *species_tree, real_t max_dist_from_species_root) {
    int i, j;
    real_t *fp;
    
    for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
        //int_array cur_interval_nodes;
        //init_int_array(&cur_interval_nodes);
        newick_node_ptr_array cur_interval_nodes;
        init_newick_node_ptr_array(&cur_interval_nodes);
        
        // species_tree is a beaded tree in here
#ifndef NDEBUG_DEVEL
        printf("\tNodes in interval [%f, %f)\n", *(fp + 1), *fp);
#endif
        add_nodes_in_interval(&cur_interval_nodes, species_tree, *(fp + 1), *fp, max_dist_from_species_root);
        
        //index each node in a hashtable
        int i = fp - spec_dists->array, j;
        for (j = 0; j < array_size(cur_interval_nodes); ++j) {
            matidx *indices = (matidx *) malloc(sizeof(matidx));
            indices->i = i, indices->j = j;
            htab_insert(mat_idx_tab, cur_interval_nodes.array[j], sizeof(newick_node), indices);
#ifndef NDEBUG_DEVEL
            hash_t h_val = htab_hash(cur_interval_nodes.array[j], sizeof(newick_node));
            printf("\t\tMapped node@%p to <%d,%d> hash: %llu\n", cur_interval_nodes.array[j], i, j, h_val);
#endif
        }
        // get nodes for the current tau
        append_newick_node_ptr_array_array(Y, cur_interval_nodes);
    }
    
    // one more iteration
    
    //int_array cur_interval_nodes;
    //init_int_array(&cur_interval_nodes);
    newick_node_ptr_array cur_interval_nodes;
    init_newick_node_ptr_array(&cur_interval_nodes);
    
    // species_tree is a beaded tree in here
#ifndef NDEBUG_DEVEL
    printf("\tNodes in interval [%f, %f)\n", 0.f, *fp);
#endif
    add_nodes_in_interval(&cur_interval_nodes, species_tree, 0.f, *fp, max_dist_from_species_root);
    
    //index each node in a hashtable
    i = fp - spec_dists->array;
    for (j = 0; j < array_size(cur_interval_nodes); ++j) {
        matidx *indices = (matidx *) malloc(sizeof(matidx));
        indices->i = i, indices->j = j;
        htab_insert(mat_idx_tab, cur_interval_nodes.array[j], sizeof(newick_node), indices);
#ifndef NDEBUG_DEVEL
        hash_t h_val = htab_hash(cur_interval_nodes.array[j], sizeof(newick_node));
        printf("\t\tMapped node@%p to <%d,%d> hash: %llu\n", cur_interval_nodes.array[j], i, j, h_val);
#endif
    }
    // get nodes for the current tau
    append_newick_node_ptr_array_array(Y, cur_interval_nodes);
    
}

int get_species_node_id_for_taxon(newick_node *n, char *t) {
    int id = -1;
    newick_child *p;

    if (t == NULL) return id;

    if (n->taxon != NULL && strcmp(n->taxon, t) == 0) {
        return n->id;
    }

    int cur_id;
    for (p = n->child; p != NULL; p = p->next) {
        cur_id = get_species_node_id_for_taxon(p->node, t);
        id = cur_id != -1 ? cur_id : id;
    }

    return id;
}

void do_get_all_descedants_taxa(newick_node *n, char_ptr_array *descedants_array_taxa) {
    newick_child *p;
    
    if (n->taxon != NULL) append_char_ptr_array(descedants_array_taxa, n->taxon);
    
    for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedants_taxa(p->node, descedants_array_taxa);
    }
}

char_ptr_array *get_all_descedants_taxa(newick_node *n) {
    newick_child *p;
    char_ptr_array *descedants_array_taxa = (char_ptr_array*) malloc(sizeof(char_ptr_array));
    init_char_ptr_array(descedants_array_taxa);

    for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedants_taxa(p->node, descedants_array_taxa);
    }

    return descedants_array_taxa;
}

void do_get_all_descedant_leaves_taxa(newick_node *n, char_ptr_array *descedants_array_taxa) {
    //printf("node %p, its child %p and taxon %s\n", n, n->child, n->taxon);
    if (n->child == NULL) append_char_ptr_array(descedants_array_taxa, n->taxon);
    
    newick_child *p;
    for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedant_leaves_taxa(p->node, descedants_array_taxa);
    }
}

char_ptr_array *get_all_descedant_leaves_taxa(newick_node *n) {
    //newick_child *p;
    char_ptr_array *descedants_array_taxa = (char_ptr_array*) malloc(sizeof(char_ptr_array));
    init_char_ptr_array(descedants_array_taxa);
    
    //for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedant_leaves_taxa(n, descedants_array_taxa);
    //}
    
    return descedants_array_taxa;
}

int get_leaves_count(newick_node *n) {
    int leaves_count = 0;
    newick_child *p;

    if (n->childNum == 0) return 1;

    for (p = n->child; p != NULL; p = p->next) {
        leaves_count += get_leaves_count(p->node);
    }
    return leaves_count;
}

newick_node *tree_from_file(const char *filename) {
    char c;
    newick_node *tree;
    char_array tree_string;

    FILE *f = fopen(filename, "r+");

    init_char_array(&tree_string);
    for (c = fgetc(f); c != EOF && c != '\n'; c = fgetc(f)) {
        append_char_array(&tree_string, c);
    }
    append_char_array(&tree_string, '\0');

    tree = parseTree(tree_string.array);

    free(tree_string.array);
    fclose(f);

    return tree;
}

struct globalArgs_t {
    const char *geneTreeFileName;
    const char *speciesTreeFileName;
    const char *outFileName;
} globalArgs;

void version(char *prog_name) {
    printf("%s, v0.1 developed by Islam Ismailov, with support of NESCENT as part of the Google Summer of Code 2012,\nMentors: J Degnan, T Stadler\n", prog_name);
}

void usage(char *prog_name) {
    printf("Usage: %s -g (--gtree) gene tree (newick tree file) -s (--stree) species tree (newick tree file) -o (--out) outputfile\n", prog_name);
}

static const char *optString = "g:s:o:vh";

static const struct option longOpts[] = {
    { "stree", required_argument, NULL, 's' },
    { "gtree", required_argument, NULL, 'g' },
    { "out", required_argument, NULL, 'o' },
    { "version", no_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 } /* long options with no short equivalent go to "case 0:" */
};

void init_global_args() {
    globalArgs.geneTreeFileName = NULL;
    globalArgs.speciesTreeFileName = NULL;
    globalArgs.outFileName = NULL;
}

void parse_global_args(int argc, char **argv) {
    int longIndex;
    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (opt != -1) {
        switch (opt) {
            case 's':
                globalArgs.speciesTreeFileName = optarg;
                break;
            case 'g':
                globalArgs.geneTreeFileName = optarg;
                break;
            case 'o':
                globalArgs.outFileName = optarg;
                break;
            case 'v':
                version(argv[0]);
                break;
            case 'h':
            case '?':
                usage(argv[0]);
                break;
            default:
                break;
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
}

int check_tree_parent_references(newick_node *n) {
    int status = 1;
    newick_child *p;
    for (p = n->child; p != NULL; p = p->next) {
        status &= (p->node->parent == n && p->node->parent != NULL);
        status &= check_tree_parent_references(p->node);
    }
    return status;
}

int check_leaves_not_null(newick_node *n) {
    int ok = 1;
    newick_child *p;
    
    for (p = n->child; p != NULL; p=p->next) {
        if (p->node->child == NULL) {
            ok &= (p->node->taxon != NULL);
        } else {
            ok &= check_leaves_not_null(p->node);
        }
        if (!ok) break;
    }
    return ok;
}

void get_dfs_leaf_order(newick_node *n, char_ptr_array *out) {
    newick_child *p;

    for (p = n->child; p != NULL; p = p->next) get_dfs_leaf_order(p->node, out);
    
    if (n->child == NULL) append_char_ptr_array(out, n->taxon);
}

int check_leaf_equivalence(char_ptr_array *species_tree_leaves, char_ptr_array *gene_tree_leaves) {
    char_ptr *p, *q;
    
    if (array_size(*species_tree_leaves) != array_size(*gene_tree_leaves)) return 0;
    
    for (p = species_tree_leaves->array; p != species_tree_leaves->last; ++p) {
        int found = 0;
        for (q = gene_tree_leaves->array; q != gene_tree_leaves->last; ++q) {
            if (strcmp(*p, *q) == 0) {
                found = 1;
                break;
            }
        }
        if (!found) return 0;
    }
    return 1;
}

int check_leaves_unique(char_ptr_array *leaves) {
    char_ptr *p, *q;
    for (p = leaves->array; p != leaves->last; ++p) {
        int unique = 1;
        for (q = p + 1; q != leaves->last; ++q) {
            if (strcmp(*p, *q) == 0) {
                unique = 0;
                break;
            }
        }
        if (!unique) return 0;
    }
    return 1;
}

int same_leaf_order(char_ptr_array *species_leaf_order, char_ptr_array *gene_leaf_order) {
    char_ptr *p, *q;
    if (array_size(*species_leaf_order) != array_size(*gene_leaf_order)) return 0;
    for (p = species_leaf_order->array, q = gene_leaf_order->array; p != species_leaf_order->last && q != gene_leaf_order->last; ++p, ++q) {
        if (strcmp(*p, *q) != 0) return 0;
    }
    return 1;
}

void taxon_swap_front(newick_node *n, char_ptr *taxon) {
    newick_child *p, *p_prev;
    for (p = n->child; p != NULL; p_prev = p, p = p->next) {
        if (p->node->taxon) {
            if (p == n->child) return;
            
            p_prev->next = p->next;
            p->next = n->child;
            n->child = p;
            break;
        }
    }
}

char_ptr_array *force_leaf_order_fixed(newick_node *n, char_ptr_array *correct_order) {
    if (n->child == NULL) {
        char_ptr_array *suborder = (char_ptr_array *) malloc(sizeof(char_ptr_array));
        init_char_ptr_array(suborder);
        append_char_ptr_array(suborder, n->taxon);

//printf("returning from leaf: %s\n", n->taxon);

        return suborder;
    }
    
    int failed = 0;
    char_ptr *i, *j;
    newick_child *p, *q, *p_prev, *q_prev;
    
    char_ptr_array *correct_suborder = (char_ptr_array *) malloc(sizeof(char_ptr_array));
    init_char_ptr_array(correct_suborder);
    
    hash_table *xtab = htab_get_new();
    
    for (p = n->child; p != NULL; p = p->next) {
        char_ptr_array *p_suborder = force_leaf_order_fixed(p->node, correct_order);
        if (p_suborder == NULL) {
            failed = 1;
            break;
        }
        htab_insert(xtab, p->node, sizeof(p->node), p_suborder);
    }
    
    if (!failed) {
        int started = 0;
        for (i = correct_order->array; i < correct_order->last; ++i) {
            int found = 0;
            for (p = n->child; p != NULL; p = p->next) {
                char_ptr_array *p_subarray = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);
                for (j = p_subarray->array; j < p_subarray->last; ++j) {

//printf("i:%s j:%s\n", *i, *j);

                    if (strcmp(*i, *j) == 0) {
                        found = 1;
                        started = 1;
                        append_char_ptr_array(correct_suborder, *i);
                        break;
                    }
                }
                if (found) break;
            }
            if (started && !found) {
                failed = 1;
                break;
            }
        }
    }
/*
    printf("BEFOR SUBORDER: ");
    for (p = n->child; p != NULL; p = p->next) {
        char_ptr_array *p_subarray = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);
        char_ptr *i;
        for (i = p_subarray->array; i < p_subarray->last; ++i) {
            printf("%s", *i);
        }
        printf(" ");
    }
    printf("\n");
    
    printf("RIGHT SUBORDER: ");
    for (i = correct_suborder->array; i < correct_suborder->last; ++i) {
        printf("%s ", *i);
    }
    printf("\n");
*/
    if (!failed) {
        p = n->child;
        p_prev = NULL;

        for (i = correct_suborder->array; i < correct_suborder->last;) {
            char_ptr_array *p_suborder = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);

#ifndef NDEBUG_DEVEL
            printf("comparing %s and %s\n", *i, *(p_suborder->array));
#endif
            if (strcmp(*i, *(p_suborder->array)) != 0) {
                // misorder detected: find suborder to swap with p
                char_ptr_array *q_suborder = NULL;
                for (q = p->next, q_prev = p; q != NULL; q_prev = q, q = q->next) {
                    q_suborder = (char_ptr_array *) htab_lookup(xtab, q->node, sizeof(q->node), compar_addr);
                    // we will compare only the first element of the suborder assuming
                    // that these suborders returned from recursive calls are
                    // ordered correctly. later we increment i with suborder size
                    if (strcmp(*i, *(q_suborder->array)) == 0) {
                        break;
                    }
                }
                
                if (q == NULL) {
                    failed = 1;
                    break;
                    //assert(q != NULL);
                } else {
#ifndef NDEBUG_DEVEL
                    printf("swapping %p with %p\n", p, q);
#endif
                    if (p == n->child) {
                        n->child = q;
                        q_prev->next = p;

                        newick_child *p_next = p->next;
                        p->next = q->next;
                        q->next = p_next;
                    } else {
                        p_prev->next = q;
                        q_prev->next = p;

                        newick_child *p_next = p->next;
                        p->next = q->next;
                        q->next = p_next;
                    }
                    p = q->next;
#ifndef NDEBUG_DEVEL
                    printf("inc i with %ld\n", array_size(*q_suborder));
#endif
                    i += array_size(*q_suborder);
                }
            } else {
#ifndef NDEBUG_DEVEL
                printf("inc i with %ld\n", array_size(*p_suborder));
#endif
                i += array_size(*p_suborder);
                p = p->next;
            }
        }
    }
    
    if (failed) {
        for (p = n->child; p != NULL; p = p->next) {
            char_ptr_array *p_subarray = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);
            if (p_subarray != NULL) {
                free(p_subarray->array);
                free(p_subarray);
            }
        }
        htab_free_table(xtab);
        free(correct_suborder->array);
        free(correct_suborder);
        return NULL;
    }

    free(correct_suborder->array);
    free(correct_suborder);

#ifndef NDEBUG_DEVEL
    printf("AFTER SUBORDER: ");
    for (p = n->child; p != NULL; p = p->next) {
        char_ptr_array *p_subarray = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);
        char_ptr *i;
        for (i = p_subarray->array; i < p_subarray->last; ++i) {
            printf("%s", *i);
        }
        printf(" ");
    }
    printf("\n");
#endif
    // everything should be sorted by this point
    // push it in a newly allocated array and return
    char_ptr_array *suborder = (char_ptr_array *) malloc(sizeof(char_ptr_array));
    init_char_ptr_array(suborder);
    
    for (p = n->child; p != NULL; p = p->next) {
        char_ptr_array *p_subarray = (char_ptr_array *) htab_lookup(xtab, p->node, sizeof(p->node), compar_addr);
        for (i = p_subarray->array; i < p_subarray->last; ++i) {
            append_char_ptr_array(suborder, *i);
        }
        free(p_subarray->array);
        free(p_subarray);
    }
    htab_free_table(xtab);

    return suborder;
}

int passes_checks(newick_node *species_tree, newick_node *gene_tree) {
    int ok = 1;

    ok &= check_leaves_not_null(species_tree);
    ok &= check_leaves_not_null(gene_tree);

    if (!ok) return ok;

    char_ptr_array species_leaf_order, gene_leaf_order;
    init_char_ptr_array(&species_leaf_order);
    init_char_ptr_array(&gene_leaf_order);
    
    get_dfs_leaf_order(species_tree, &species_leaf_order);
    get_dfs_leaf_order(gene_tree, &gene_leaf_order);

    ok &= check_leaves_unique(&species_leaf_order);
    ok &= check_leaves_unique(&species_leaf_order);

    if (!ok) return ok;

    ok &= check_leaf_equivalence(&species_leaf_order, &gene_leaf_order);
/*
    // rotation implementation is not powerful enough (yet)
    if (!same_leaf_order(&species_leaf_order, &gene_leaf_order)) {
        char_ptr_array *ret = force_leaf_order_fixed(gene_tree, &species_leaf_order);
        ok &= (ret != NULL);

        if (ret != NULL) {
            free(ret->array);
            free(ret);
        }
    }
*/
    return ok;
}

int main(int argc, char **argv) {
    newick_node *species_tree = NULL;
    newick_node *gene_tree = NULL;

    real_array *spec_dists;
    int_array *gene_lineages;
    int coalescence_count;
    node2real_array *coalescence_array;
    node2int_array *species_indexed_nodes;
    int *ip;
    int i, j, k, n;
    real_t *fp;
    node2real *n2f;

    monitored_memory_init();

    init_global_args();
    parse_global_args(argc, argv);

    if (globalArgs.geneTreeFileName == NULL || globalArgs.speciesTreeFileName == NULL || globalArgs.outFileName == NULL) {
        usage(argv[0]);
        return -1;
    }

    gene_tree = tree_from_file(globalArgs.geneTreeFileName);
    species_tree = tree_from_file(globalArgs.speciesTreeFileName);

    if (check_tree_parent_references(species_tree) == 0 || check_tree_parent_references(gene_tree) == 0) {
        printf("BAD PARENTING\n");
        return EXIT_FAILURE;
    }

#ifndef NDEBUG
    printf("\nGene tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printTree(gene_tree);

    printf("\n\nSpecies tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printTree(species_tree);
#endif
/*
    char_ptr_array species_leaf_order, gene_leaf_order;
    
    init_char_ptr_array(&species_leaf_order);
    init_char_ptr_array(&gene_leaf_order);
    
    get_dfs_leaf_order(species_tree, &species_leaf_order);
    get_dfs_leaf_order(gene_tree, &gene_leaf_order);
    
    printf("\n\nspecies leaf order\n");
    char_ptr *x;
    for (x = species_leaf_order.array; x != species_leaf_order.last; ++x) {
        printf("%s ", *x);
    }
    printf("\n");

    printf("gene leaf order\n");
    for (x = gene_leaf_order.array; x != gene_leaf_order.last; ++x) {
        printf("%s ", *x);
    }
    printf("\n");
    
    force_leaf_order_fixed(gene_tree, &species_leaf_order);
    exit(-1);
*/
    if (!passes_checks(species_tree, gene_tree)) {
        puts("BAD INPUT");
        exit(-1);
    }

    //max_dist_from_species_root
    real_t max_dist_from_species_root = max_dist_from_root(species_tree);
    real_t max_dist_from_gene_root = max_dist_from_root(gene_tree);

#ifndef NDEBUG
    printf("\n\nMax distance from root (height) for species tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printf("\t%f\n", max_dist_from_species_root);
    printf("\n\nMax distance from root (height) for gene tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printf("\t%f\n", max_dist_from_gene_root);
#endif

    // find speciation distances (0 is max distance from species root)
    spec_dists = get_speciation_distances(species_tree, max_dist_from_species_root);

#ifndef NDEBUG
    printf("\n\nSpeciation distances:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (fp = spec_dists->array; fp != spec_dists->last; ++fp) {
        printf("\t%f\n", *fp);
    }
#endif

#ifndef NDEBUG
    printf("\n\nGene Lineages:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

    gene_lineages = get_gene_lineages(spec_dists, gene_tree, max_dist_from_gene_root);

#ifndef NDEBUG
    for (ip = gene_lineages->array; ip != gene_lineages->last; ++ip) {
        printf("\t%d\n", *ip);
    }
#endif

    coalescence_count = get_tree_coalescence_count(gene_tree);

    // get coalescence array of the gene tree (find all nodes that have children)
    coalescence_array = get_coalescence_array(gene_tree);

#ifndef NDEBUG_DEVEL
    printf("\n\nIndexed Coalescence array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (n2f = coalescence_array->array; n2f != coalescence_array->last; ++n2f) {
        printf("\tnode@%p with dist from root %f(%ld array idx) childnum %d taxon %s\n", n2f->node, n2f->val, (n2f - coalescence_array->array), n2f->node->childNum, n2f->node->taxon);
    }
#endif
    
    // assign each node of species tree integer id (sorted by distance from the root)
    // to fast search for our lca implementation
    species_indexed_nodes = get_indexed_array(species_tree);

#ifndef NDEBUG_DEVEL
    printf("\n\nIndexed species array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    node2int *n2i;
    for (n2i = species_indexed_nodes->array; n2i != species_indexed_nodes->last; ++n2i) {
        printf("\tnode@%p ith node id %d and array index %ld taxon %s\n", n2i->node, n2i->node->id, (n2i - species_indexed_nodes->array), n2i->node->taxon);
    }
#endif

#ifndef NDEBUG_DEVEL
    printf("\n\nLCA preprocessing:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif
    lca_init(species_indexed_nodes->last - species_indexed_nodes->array, species_indexed_nodes);

    // for each node in a gene tree we need to get all descendants' taxa, so we can find
    // the lca of the nodes in species tree with the equivalent taxa
    int_array *equivalent_node_ids = (int_array *) malloc(sizeof(int_array));
    init_int_array(equivalent_node_ids);

    int lca_idx;

    int_array g;
    init_int_array(&g);
    
    int species_leaves_count = get_leaves_count(species_tree);
    int gene_leaves_count = get_leaves_count(gene_tree);
    
    assert(species_leaves_count == gene_leaves_count);

    n = species_leaves_count;
    for (i = 0; i < n; ++i) {
        int min_lineages = 0;

        for (j = i + 1; j < n - 1; ++j) {
            real_t product = 1.0;
            for (k = j; k < n - 1; ++k) {
                //printf("i:%d, j:%d, k:%d\n", i, j, k);

                n2f = (coalescence_array->array + k);

#ifndef NDEBUG_DEVEL
                printf("\n\nAll descedants taxa and their equivalent ids for node@%p:\n---- ---- ---- ---- ---- ---- ---- ----\n", n2f->node);
#endif

                char_ptr_array *taxa = get_all_descedants_taxa(n2f->node);

                clear_int_array(equivalent_node_ids);

                int id;
                char **taxon;
                for (taxon = taxa->array; taxon != taxa->last; ++taxon) {
                    id = get_species_node_id_for_taxon(species_tree, *taxon);
                    if (id != -1) append_int_array(equivalent_node_ids, id);

#ifndef NDEBUG_DEVEL
                    printf("\ttaxon:%s id:%d\n", *taxon == NULL ? "NULL" : *taxon, id);
#endif
                }

#ifndef NDEBUG_DEVEL
                printf("\n\tLCA search:\n");
#endif

                // now we need to find lowest common ancestor for these nodes
                lca_idx = *(equivalent_node_ids->array);
                for (ip = (equivalent_node_ids->array + 1); ip != equivalent_node_ids->last; ++ip) {
#ifndef NDEBUG_DEVEL
                    printf("\t\tLCA(%d, %d): ", lca_idx, *ip);
#endif
                    lca_idx = lca(lca_idx, *ip);
#ifndef NDEBUG_DEVEL
                    printf("%d\n", lca_idx);
#endif
                }

                // let's find tau interval for a given lowest common ancestor
                // and to do that we need to calculate it's distance from the max distance from species root
                real_t lca_dist = max_dist_from_species_root - get_distance_from_root((species_indexed_nodes->array + lca_idx)->node);
                for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
//                    if (lca_dist >= *(fp + 1) && lca_dist < *fp) {
                    if (real_cmp(lca_dist, *(fp + 1)) >= 0 && real_cmp(lca_dist, *fp) < 0) {
                        break;
                    }
                }

                int tau_idx = fp - spec_dists->array;

#ifndef NDEBUG_DEVEL
                printf("\n\tLCA Overall: %d with Tau index %d\n", lca_idx, tau_idx);
#endif

                product *= (tau_idx > i) ? 1 : 0;
            }
            min_lineages += product;
        }
        append_int_array(&g, n - min_lineages);
    }

#ifndef NDEBUG
    printf("\n\ng array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (ip = g.array; ip != g.last; ++ip) {
        printf("\tg[%ld]:%d\n", ip - g.array, *ip);
    }
#endif

#ifndef NDEBUG
    printf("\n\nBeaded species tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif
    bead_tree(species_tree, spec_dists, max_dist_from_species_root);
    printTree(species_tree);
    
#ifndef NDEBUG_DEVEL
    printf("\n\nY matrix construction:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif
    
    hash_table *mat_idx_tab = htab_get_new(); // node -> <i,j> mapping for Y array
    newick_node_ptr_array_array Y;
    init_newick_node_ptr_array_array(&Y);
    construct_y_matrix(mat_idx_tab, &Y, spec_dists, species_tree, max_dist_from_species_root);

#ifndef NDEBUG
    printf("\n\nY matrix:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (i = 0; i < array_size(Y); ++i) {
        for (j = 0; j < array_size(Y.array[i]); ++j) {
            printf("\tY[%d][%d]: node@%p taxon:%s\n", i, j, Y.array[i].array[j], Y.array[i].array[j]->taxon);
        }
    }
#endif

#ifndef NDEBUG
    printf("\n\nm array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif
    
    // m[i] arrays contains coalescenses that occur in tau[i] intervals
    node2real_array_array m;
    init_node2real_array_array(&m);

    for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
        int coalescence_events_count = 0;
        node2real_array coalescences;
        init_node2real_array(&coalescences);
        for (n2f = coalescence_array->array; n2f != coalescence_array->last; ++n2f) {
            real_t dist = max_dist_from_gene_root - n2f->val;
#ifndef NDEBUG_DEVEL
            printf("\t%f =< %f < %f?", *(fp + 1), dist, *fp);
#endif
//            if (dist >= *(fp + 1) && dist < *fp) {
            if (real_cmp(dist, *(fp + 1)) >= 0 && real_cmp(dist, *fp) < 0) {
#ifndef NDEBUG_DEVEL
                printf(" y\n");
#endif
                node2real coalescence;
                coalescence.node = n2f->node;
                coalescence.val = dist;
#ifndef NDEBUG_DEVEL
                printf("m[%d][%d]=%f\n", m.last - m.array, coalescences.last - coalescences.array, dist);
#endif
                append_node2real_array(&coalescences, coalescence);
                ++coalescence_events_count;
            }
#ifndef NDEBUG_DEVEL
            else {
                printf(" n\n");
            }
#endif
        }
        append_node2real_array_array(&m, coalescences);
#ifndef NDEBUG
        printf("\t%d coalescence events in interval m[%ld]: [%f, %f)\n", coalescence_events_count, m.last - m.array - 1, *(fp + 1), *fp);
#endif
    }
    
    node2int_array *gene_indexed_nodes = get_indexed_array(gene_tree);
    node2int_array *beaded_indexed_nodes = get_indexed_array(species_tree);

#ifndef NDEBUG
    printf("\n\nGene descendants:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

    hash_table *gene_leaf_descendants = htab_get_new();
    for (i = 0; i < array_size(*gene_indexed_nodes); ++i) {
        char_ptr_array *descendants = get_all_descedant_leaves_taxa(gene_indexed_nodes->array[i].node);

#ifndef NDEBUG
        printf("\t%p: ", gene_indexed_nodes->array[i].node);
        char_ptr *x;
        for (x = descendants->array; x != descendants->last; ++x) {
            printf("%s ", *x);
        }
        printf("\n");
#endif
        
        htab_insert(gene_leaf_descendants, gene_indexed_nodes->array[i].node, sizeof(gene_indexed_nodes->array[i].node), descendants);
    }

#ifndef NDEBUG
    printf("\n\nSpecies descendants:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

    hash_table *species_leaf_descendants = htab_get_new();
    for (i = 0; i < array_size(*beaded_indexed_nodes); ++i) {
        char_ptr_array *descendants = get_all_descedant_leaves_taxa(beaded_indexed_nodes->array[i].node);

#ifndef NDEBUG
        printf("\t%p: ", beaded_indexed_nodes->array[i].node);
        char_ptr *x;
        for (x = descendants->array; x != descendants->last; ++x) {
            printf("%s ", *x);
        }
        printf("\n");
#endif

        htab_insert(species_leaf_descendants, beaded_indexed_nodes->array[i].node, sizeof(beaded_indexed_nodes->array[i].node), descendants);
    }
    
#ifndef NDEBUG_DEVEL
    printf("\n\nk array constuction:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

    // k[][][] array
    int speciation_count = spec_dists->last - spec_dists->array; // this is i dimension
    coalescence_count = get_tree_coalescence_count(gene_tree);   // this is j dimension

    int z;
    int_array_array_array K;
    init_int_array_array_array(&K);
    for (i = 0; i < speciation_count; ++i) {
        int_array_array mtx_to_add;
        init_int_array_array(&mtx_to_add);
        for (j = 0; j < coalescence_count; ++j) {
            int_array to_add;
            init_int_array(&to_add);
            for (z = 0; z < array_size(Y.array[i]); ++z) {
                int int_to_add = 0;
                append_int_array(&to_add, int_to_add);
            }
            append_int_array_array(&mtx_to_add, to_add);
        }
        append_int_array_array_array(&K, mtx_to_add);
    }

#ifndef NDEBUG_DEVEL
    printf("speciation count: %d\n", speciation_count);
#endif

    for (i = 0; i < speciation_count; ++i) {

        int arsize = Y.array[i].last - Y.array[i].array;
        
        assert(arsize == array_size(Y.array[i]));
//        printf("y[%d] has %ld elements\n", i, array_size((Y.array)[i]));
        
        for (z = 0; z < array_size(Y.array[i]); ++z) {
#ifndef NDEBUG_DEVEL
            printf("k[%d][0][%d]:\n", i, z);
#endif
            //int_array topology_prefix = get_topology_prefix(species_tree, Y.array[i].array[z]);
#ifndef NDEBUG_DEVEL
            printf("run for spec_dists[%d]=%f\n", i, (spec_dists->array)[i]);
#endif
//            K.array[i].array[0].array[z] = get_exit_branches(gene_tree, (spec_dists->array)[i - 1], &topology_prefix, max_dist_from_gene_root);
            char_ptr_array *descendants = htab_lookup(species_leaf_descendants, Y.array[i].array[z], sizeof(Y.array[i].array[z]), compar_addr);

#ifndef NDEBUG_DEVEL
            printf("looking up %p\n", Y.array[i].array[z]);
#endif

            assert(descendants != NULL);
            K.array[i].array[0].array[z] = get_exit_branches_fixed(gene_tree, spec_dists->array[i], max_dist_from_gene_root, descendants, gene_leaf_descendants);
        }
    }

    for (i = 0; i < array_size(Y); ++i) {
        for (z = 0; z < array_size(Y.array[i]); ++z) {
            // start j = 1 because j = 0 already calculated
            for (j = 1; j <= array_size(m.array[i]); ++j) {
                real_t dist = m.array[i].array[j - 1].val;
                real_t EPS = 2 * fmax(REAL_T_EPS, REAL_T_EPS * fabs(dist));

#ifndef NDEBUG_DEVEL
                printf("k[%d][%d][%d]: with distance %f\n", i, j, z, dist);
#endif

                K.array[i].array[j].array[z] = get_exit_branches_fixed(gene_tree, dist - EPS, max_dist_from_gene_root, htab_lookup(species_leaf_descendants, Y.array[i].array[z], sizeof(Y.array[i].array[z]), compar_addr), gene_leaf_descendants);
            }
        }
    }

#ifndef NDEBUG
    printf("\n\nk array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    
    for (i = 0; i < array_size(Y); ++i) {
        for (z = 0; z < array_size(Y.array[i]); ++z) {
            // start j = 1 because j = 0 already calculated
            for (j = 0; j <= array_size(m.array[i]); ++j) {
                printf("\tK[%d][%d][%d] = %d\n", i, j, z, K.array[i].array[j].array[z]);
            }
        }
    }
#endif

    int_array_array lambda;
    init_int_array_array(&lambda);
    
    for (i = 0; i < array_size(Y); ++i) {
        int_array lambda_i;
        init_int_array(&lambda_i);
        for (j = 0; j <= array_size(m.array[i]); ++j) {
            int sum = 0;
            for (z = 0; z <= i; ++z) {
                sum += combinations(K.array[i].array[j].array[z], 2);
            }
            append_int_array(&lambda_i, sum);
        }
        append_int_array_array(&lambda, lambda_i);
    }

#ifndef NDEBUG
    printf("\n\nlambda array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (i = 0; i < array_size(Y); ++i) {
        for (j = 0; j <= array_size(m.array[i]); ++j) {
            printf("\tlambda[%d][%d] = %d\n", i, j, lambda.array[i].array[j]); 
        }
    }
#endif
    
    printf("\n\nfinal steps:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    
    int l_1;
    n = speciation_count + 1;
    real_t P = 0;
    for (l_1 = g.array[0]; l_1 < n; ++l_1) {
        P += theorem2_func(0, l_1, n, &g, &lambda, spec_dists, &m) / H(l_1);
    }
    
    printf("\tfinal answer: %f\n", P);
    
    lca_end();
    monitored_memory_end();
    return 0;
}
