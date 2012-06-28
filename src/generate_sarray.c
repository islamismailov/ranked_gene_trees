#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "newick_tree.h"

#include "utils.h"
#include "getopt.h"
#include "hash_table.h"
#include "monitored_memory.h"

float max(float x, float y) {
    return x > y? x : y;
}

typedef struct node2float {
    newick_node *node;
    float val;
} node2float;

typedef struct node2int {
    newick_node *node;
    int val;
} node2int;

int node2float_cmp(node2float *a, node2float *b) {
    float diff = a->val - b->val;

    if (diff > 1e-6)        return  1;
    else if (diff < -1e-6)  return -1;

    return 0;
}

typedef char * char_ptr;
typedef int * int_ptr;
typedef int ** int_ptr_ptr;
typedef int *** int_ptr_ptr_ptr;
typedef newick_node * newick_node_ptr;

DEF_ARRAY(int);
DEF_ARRAY(char);
DEF_ARRAY(float);

DEF_ARRAY(newick_node_ptr);
DEF_ARRAY(newick_node_ptr_array);

DEF_ARRAY(node2float);
DEF_ARRAY(node2int);
DEF_ARRAY(char_ptr);
DEF_ARRAY(int_array);
DEF_ARRAY(int_array_array);

float get_distance_from_root(newick_node *n) {
    float dist = n->dist;
    if (n->parent != NULL) dist += get_distance_from_root(n->parent);
    return dist;
}

float max_dist_from_root(newick_node *t) {
    newick_child *p;
    float dist = 0.f;
    for (p = t->child; p != NULL; p = p->next) {
        dist = max(dist, max_dist_from_root(p->node));
    }
    return dist + t->dist;
}

void do_get_speciation_distances(newick_node *t, float distance, float_array *speciation_distances, float max_dist_from_root) {
    newick_child *p;
    if (t->childNum > 0) // speciation happened
        append_float_array(speciation_distances, max_dist_from_root - (distance + t->dist));
    for (p = t->child; p != NULL; p = p->next) {
        do_get_speciation_distances(p->node, distance + t->dist, speciation_distances, max_dist_from_root);
    }
}

float_array *get_speciation_distances(newick_node *t, float max_dist_from_root) {
    float_array *speciation_distances = (float_array *)malloc(sizeof(float_array));
    init_float_array(speciation_distances);
    do_get_speciation_distances(t, 0.0f, speciation_distances, max_dist_from_root);
    append_float_array(speciation_distances, FLT_MAX);

    qsort(speciation_distances->array, speciation_distances->last - speciation_distances->array,
        sizeof(float), (int(*)(const void*,const void*))flt_cmp_desc);
    return speciation_distances;
}

int do_get_gene_lineages(newick_node *t, float limit, float distance, float max_dist_from_root) {
    int lineages = 0;
    newick_child *p;
    printf("\t%f <= %f? ", max_dist_from_root - distance - t->dist, limit);
    if ((max_dist_from_root - distance - t->dist) <= limit) {
        puts("yes");
        return lineages + 1;
    } else puts("no");

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist, max_dist_from_root);
    }

    return lineages;
}

int_array *get_gene_lineages(float_array *speciation_distances, newick_node *t, float max_dist_from_root) {
    float *p;
    int_array *lineages = (int_array *)malloc(sizeof(int_array));

    init_int_array(lineages);
    for (p = speciation_distances->array; p != speciation_distances->last; ++p) {
        append_int_array(lineages, do_get_gene_lineages(t, *p, 0.f, max_dist_from_root));
    }
    return lineages;
}

int do_get_gene_lineages_for_k(newick_node *t, float limit, float distance, float max_dist_from_root, int_array *species_turns, int_array *gene_turns, int turn) {
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

    printf("\t%f <= %f? ", max_dist_from_root - distance - t->dist, limit);

    if ((max_dist_from_root - distance - t->dist) <= limit) {
        puts("yes");
        return lineages + 1;
    } else puts("no");

    int child_index = 0;
    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages_for_k(p->node, limit, distance + t->dist, max_dist_from_root, species_turns, gene_turns, child_index++);
    }

    // remove added turn
    if (gene_turns->last - gene_turns->array > 0) --(gene_turns->last);

    return lineages;
}

int get_gene_lineages_for_k(float *speciation_distance, newick_node *t, float max_dist_from_root, int_array *species_turns) {
    int lineages = 0;
    int_array *gene_turns = (int_array *) malloc(sizeof(int_array));
    init_int_array(gene_turns);
    lineages = do_get_gene_lineages_for_k(t, *speciation_distance, 0.f, max_dist_from_root, species_turns, gene_turns, -1);
    return lineages;
}

/*
 * distances from the root for each node
 */
void do_get_distance_array(newick_node *t, float distance, node2float_array *dist_array) {
    newick_child *p;
    node2float pair;
    pair.val = distance + t->dist;
    pair.node = t;
    append_node2float_array(dist_array, pair);
    for (p = t->child; p != NULL; p = p->next) {
        do_get_distance_array(p->node, distance + t->dist, dist_array);
    }
}

// TODO: remove later, used for debugging
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
 * nodes are sorted by distance
 */
node2int_array *get_coalescence_array(newick_node *t) {
    node2float *n2f;
    node2float_array *dist_array = (node2float_array *)malloc(sizeof(node2float_array));
    node2int_array *coalescence_array = (node2int_array *)malloc(sizeof(node2int_array));

    init_node2float_array(dist_array);
    init_node2int_array(coalescence_array);

    do_get_distance_array(t, 0.f, dist_array);
    qsort(dist_array->array, dist_array->last - dist_array->array,
        sizeof(node2float), (int(*)(const void*,const void*))node2float_cmp);

    int c = 0;
    for (n2f = dist_array->array; n2f != dist_array->last; ++n2f) {
        if (n2f->node->childNum > 0) {
            node2int pair;
            pair.node = n2f->node;
            pair.val = c++;

            append_node2int_array(coalescence_array, pair);
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
    do_get_topology_prefix(&topology_prefix, n, target);
    return topology_prefix;
}
/*
 *
int do_get_gene_lineages(newick_node *t, float limit, float distance, float max_dist_from_root) {
    int lineages = 0;
    newick_child *p;
    printf("\t%f <= %f? ", max_dist_from_root - distance - t->dist, limit);
    if ((max_dist_from_root - distance - t->dist) <= limit) {
        puts("yes");
        return lineages + 1;
    } else puts("no");

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist, max_dist_from_root);
    }

    return lineages;
*/

int do_get_exit_branches(newick_node *n, float limit, float distance, int_array *topology_prefix, int_array *self_topology_prefix, float max_dist_from_root) {
    int exit_branches_count = 0;
    newick_child *p;

#ifndef NDEBUG
    printf("\t%f <= %f? ", max_dist_from_root - distance - n->dist, limit);
#endif

    if ((max_dist_from_root - distance - n->dist) <= limit) {

#ifndef NDEBUG
        puts("yes");
        printf("\tchecking for topology prefix...\n");
        printf("\t\tself topology prefix:\n\t\t\t");
#endif
        int *ip, *jp;
#ifndef NDEBUG
        for (ip = self_topology_prefix->array; ip != self_topology_prefix->last; ++ip) {
            printf("%d", *ip);
        }
        printf("\n");
        printf("\t\tspecies topology prefix:\n\t\t\t");

        for (jp = topology_prefix->array; jp != topology_prefix->last; ++jp) {
            printf("%d", *jp);
        }
        printf("\n");
#endif
        int is_subset = 0;
        for (ip = self_topology_prefix->array, jp = topology_prefix->array;
                ip != self_topology_prefix->last && jp != topology_prefix->last;
                ++ip, ++jp) {
            if (*ip != *jp) {
                is_subset = 0;
                break;
            }
        }
#ifndef NDEBUG
        if (is_subset) printf("subset!\n"); else printf("not subset!\n");
#endif
        if (is_subset) return exit_branches_count + 1;
    } else {
#ifndef NDEBUG
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

int get_exit_branches(newick_node *n, float distance, int_array *topology_prefix, float max_dist_from_root) {
    int_array self_topology_prefix;
    init_int_array(&self_topology_prefix);
    int exit_branches_count = do_get_exit_branches(n, distance, 0.f, topology_prefix, &self_topology_prefix, max_dist_from_root);
    free(self_topology_prefix.array);
    return exit_branches_count;
}

void add_nodes_in_interval(newick_node_ptr_array *arr, newick_node *n, float start_interval, float end_interval, float max_root_distance) {
    //TODO: optimize, can call "get_distance_from_root" only once (in driver, move this to do_driver)
    float dist = max_root_distance - get_distance_from_root(n);
    if (dist >= start_interval && dist < end_interval) {
        append_newick_node_ptr_array(arr, n);
    }

    newick_child *p;
    for (p = n->child; p != NULL; p = p->next) {
        add_nodes_in_interval(arr, p->node, start_interval, end_interval, max_root_distance);
    }
}

void do_bead_tree(newick_node *t, float distance, float_array *speciation_distances, float max_dist_from_root) {
    float *fp;

    newick_node *bead;
    newick_child *p, *q, *child, *child_head, *new_child;

    float start_interval, end_interval;

    start_interval = max_dist_from_root - (distance + t->dist);
    for (p = t->child; p != NULL; ) {
        do_bead_tree(p->node, distance + t->dist, speciation_distances, max_dist_from_root); //do_get_speciation_distances(p->node, distance + t->dist, speciation_distances, max_dist_from_root);

        // find if there are any speciation times in between this node and it's child
        end_interval = max_dist_from_root - (distance + t->dist + p->node->dist);

        child = p;
        child_head = t->child;

        p = p->next;

        for (fp = (speciation_distances->last - 1); fp != (speciation_distances->array - 1); --fp) {
            if (*fp > start_interval && *fp < end_interval) {
                // delete son from the old node
                if (child_head->next == p) {
                    child_head->next = child->next;
                } else {
                    for (q = child_head; q != NULL; q = q->next) {
                        if (q->next == p) break;
                    }
                    q->next = child->next;
                }

                bead = (newick_node *) monitored_malloc(sizeof(newick_node));
                bead->dist = *fp;
                bead->child = child;
                bead->child->node->parent = bead;
                bead->childNum = 1;
                bead->child->next = NULL;

                // now attach bead as one of t's children:
                new_child = (newick_child *) monitored_malloc(sizeof(newick_child));
                new_child->node = bead;
                new_child->next = t->child;
                t->child = new_child;
                bead->parent = t;

                // parents number of children stayed same

                child = new_child;
                child_head = new_child;
                end_interval = *fp;
            }
        }
    }
}

void bead_tree(newick_node *n, float_array *speciation_times, float max_dist_from_root) {
    do_bead_tree(n, 0.f, speciation_times, max_dist_from_root);
}

/*
 * nodes array
 * nodes are sorted by distance
 */
node2int_array *get_indexed_array(newick_node *t) {
    node2float *n2f;
    node2float_array *dist_array = (node2float_array *)malloc(sizeof(node2float_array));
    node2int_array *nodes_array = (node2int_array *)malloc(sizeof(node2int_array));

    init_node2float_array(dist_array);
    init_node2int_array(nodes_array);

    do_get_distance_array(t, 0.f, dist_array);
    qsort(dist_array->array, dist_array->last - dist_array->array,
        sizeof(node2float), (int(*)(const void*,const void*))node2float_cmp);

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










int timer, n, l;
int *tin, *tout;
int **up;

void lca_preprocess(node2int_array *coalescence_array, int v, int p) {
    int i, to;
    node2int *cp;
    newick_child *np;

    tin[v] = ++timer;
    up[v][0] = p;
    for (i = 1; i <= l; ++i)
        up[v][i] = up[up[v][i - 1]][i - 1];

    // iterate through children of coalescence u[v]:

    for (np = coalescence_array->array[v].node->child; np != NULL; np = np->next) {
        printf("\titerating thru children of u[%d]@%p:\n", v, coalescence_array->array[v].node);
        // find child's index (to)
        for (cp = coalescence_array->array, to = 0; cp != coalescence_array->last; ++to, ++cp) {
            //if (cp->node == coalescence_array->array[v].node)
            if (cp->node == np->node)
                break;
        }

        printf("\tu[%d](u[%d])@%p is a child? ", to, cp->val, np->node);
        if (to != p && to != coalescence_array->last - coalescence_array->array) {
            puts("y");
            lca_preprocess(coalescence_array, cp->val, v);
        } else puts("n");
    }

    tout[v] = ++timer;
}

int upper(int a, int b) {
    return tin[a] <= tin[b] && tout[a] >= tout[b];
}

int lca(int a, int b) {
    int i;
    if (upper(a, b)) return a;
    if (upper(b, a)) return b;

    for (i = l; i >= 0; --i) {
        if (!upper(up[a][i], b)) {
            a = up[a][i];
        }
    }

    return up[a][0];
}

void do_get_all_descedants_taxa(newick_node *n, char_ptr_array *descedants_array_taxa) {
    newick_child *p;

    if (n->taxon != NULL) append_char_ptr_array(descedants_array_taxa, n->taxon);

    for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedants_taxa(p->node, descedants_array_taxa);
    }
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



char_ptr_array *get_all_descedants_taxa(newick_node *n) {
    newick_child *p;
    char_ptr_array *descedants_array_taxa = (char_ptr_array*) malloc(sizeof(char_ptr_array));
    init_char_ptr_array(descedants_array_taxa);

    for (p = n->child; p != NULL; p = p->next) {
        do_get_all_descedants_taxa(p->node, descedants_array_taxa);
    }

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

void lca_init(int n, node2int_array *coalescence_array) {
    int i;

    tin = (int *) malloc(n * sizeof(int));
    tout = (int *) malloc(n * sizeof(int));
    up = (int **) malloc(n * sizeof(int *));

    l = 1;
    while ((1 << l) <= n)  ++l;
    for (i = 0; i < n; ++i)  up[i] = (int *) malloc((l + 1) * sizeof(int));

    lca_preprocess(coalescence_array, 0, 0);
}

void lca_end() {
    int i;
    free(tin);
    free(tout);
    for (i = 0; i < n; ++i) free(up[i]);
    free(up);
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

int main(int argc, char **argv) {
    newick_node *species_tree = NULL;
    newick_node *gene_tree = NULL;

    float_array *spec_dists;
    int_array *gene_lineages;
    int coalescence_count;
    node2int_array *coalescence_array;
    node2int_array *species_indexed_nodes;
    int *ip;
    int i, j, k, n;
    float *fp;
    node2int *n2i;
    float farthest_leaf_dist = 0.f;

    monitored_memory_init();

    init_global_args();
    parse_global_args(argc, argv);

    if (globalArgs.geneTreeFileName == NULL || globalArgs.speciesTreeFileName == NULL || globalArgs.outFileName == NULL) {
        usage(argv[0]);
        return -1;
    }

    gene_tree = tree_from_file(globalArgs.geneTreeFileName);
    species_tree = tree_from_file(globalArgs.speciesTreeFileName);

#ifndef NDEBUG
    printf("\n\nGene tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printTree(gene_tree);

    printf("\n\nSpecies tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printTree(species_tree);
#endif

    farthest_leaf_dist = max_dist_from_root(species_tree);

#ifndef NDEBUG
    printf("\n\nMax distance from root for species tree:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    printf("\t%f\n", farthest_leaf_dist);
#endif

    // find speciation distances (0 is farthest leaf from the root)
    spec_dists = get_speciation_distances(species_tree, farthest_leaf_dist);

#ifndef NDEBUG
    printf("\n\nSpeciation distances:\n---- ---- ---- ---- ---- ---- ---- ----\n");
    for (fp = spec_dists->array; fp != spec_dists->last; ++fp) {
        printf("\t%f\n", *fp);
    }
#endif

#ifndef NDEBUG
    printf("\n\nGene Lineages:\n");
#endif

    gene_lineages = get_gene_lineages(spec_dists, gene_tree, farthest_leaf_dist);

#ifndef NDEBUG
    for (ip = gene_lineages->array; ip != gene_lineages->last; ++ip) {
        printf("\t%d\n", *ip);
    }
#endif

    coalescence_count = get_tree_coalescence_count(gene_tree);

    // get coalescence array of the gene tree (find all nodes that have children)
    coalescence_array = get_coalescence_array(gene_tree);

#ifndef NDEBUG
    printf("\n\nCoalescence array (size %d):\n---- ---- ---- ---- ---- ---- ---- ----\n", coalescence_count);
    for (n2i = coalescence_array->array; n2i != coalescence_array->last; ++n2i) {
        printf("\tval:%d node@%p childnum:%d\n", n2i->val, n2i->node, n2i->node->childNum);
    }
#endif

    //lca_init(coalescence_count, root, coalescence_array);
    //printf("lca of %d and %d is %d\n", 2, 3, lca (2, 3));
    //lca_end();

#ifndef NDEBUG
    printf("\n\nIndexed species array:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

    // assign each node of species tree integer id (sorted by distance from the root)
    // to fast search for our lca implementation
    species_indexed_nodes = get_indexed_array(species_tree);

#ifndef NDEBUG
    for (n2i = species_indexed_nodes->array; n2i != species_indexed_nodes->last; ++n2i) {
        printf("\tnode@%p with node id %d and array index %d\n", n2i->node, n2i->node->id, (n2i - species_indexed_nodes->array));
    }
#endif

    // this array holds number of coalescence events in interval tau[i]
    // TODO: add m array

#ifndef NDEBUG
    printf("\n\nLCA:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif
    lca_init(species_indexed_nodes->last - species_indexed_nodes->array, species_indexed_nodes);

    // for each node in a gene tree we need to get all descendants' taxa, so we can find
    // the lca of the nodes in species tree with the equivalent taxa
    int_array *equivalent_node_ids = (int_array *) malloc(sizeof(int_array));
    init_int_array(equivalent_node_ids);

    int lca_idx;
/*
    for (n2i = coalescence_array->array; n2i != coalescence_array->last; ++n2i) {

#ifndef NDEBUG
        printf("\n\nAll descedants taxa and their equivalent ids for node@%p:\n---- ---- ---- ---- ---- ---- ---- ----\n", n2i->node);
#endif
        char_ptr_array *taxa = get_all_descedants_taxa(n2i->node);

        clear_int_array(equivalent_node_ids);

        int id;
        char **taxon;
        for (taxon = taxa->array; taxon != taxa->last; ++taxon) {
            id = get_species_node_id_for_taxon(species_tree, *taxon);
            if (id != -1) append_int_array(equivalent_node_ids, id);

#ifndef NDEBUG
            printf("\ttaxon:%s id:%d\n", *taxon == NULL ? "NULL" : *taxon, id);
#endif
        }

#ifndef NDEBUG
        printf("\n\nLCA search:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

        // now we need to find lowest common ancestor for these nodes
        lca_idx = *(equivalent_node_ids->array);
        for (ip = (equivalent_node_ids->array + 1); ip != equivalent_node_ids->last; ++ip) {
            lca_idx = lca(lca_idx, *ip);

#ifndef NDEBUG
            printf("\tLCA(%d, %d): %d\n", lca_idx, *ip, lca_idx);
#endif
        }

        // let's find tau interval for a given lowest common ancestor
        // so that we need to calculate it's distance from the farthest leaf from the root
        float lca_dist = farthest_leaf_dist - get_distance_from_root((species_indexed_nodes->array + lca_idx)->node);
        for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
            if (lca_dist >= *(fp + 1) && lca_dist < *fp) {
                break;
            }
        }

        int tau_idx = fp - spec_dists->array;

#ifndef NDEBUG
        printf("---- ---- ---- ---- ---- ---- ---- ----\n\tOverall LCA id: %d with Tau index %d\n", lca_idx, tau_idx);
#endif
    }
*/
    int_array g;
    init_int_array(&g);

    n = get_leaves_count(species_tree); // should be same for gene_tree
    for (i = 0; i < n; ++i) {
        int min_lineages = 0;

        for (j = i + 1; j < n - 1; ++j) {
            double product = 1.0;
            for (k = j; k < n - 1; ++k) {
                printf("i:%d, j:%d, k:%d\n", i, j, k);

                n2i = (coalescence_array->array + k);

#ifndef NDEBUG
                printf("\n\nAll descedants taxa and their equivalent ids for node@%p:\n---- ---- ---- ---- ---- ---- ---- ----\n", n2i->node);
#endif

                char_ptr_array *taxa = get_all_descedants_taxa(n2i->node);

                clear_int_array(equivalent_node_ids);

                int id;
                char **taxon;
                for (taxon = taxa->array; taxon != taxa->last; ++taxon) {
                    id = get_species_node_id_for_taxon(species_tree, *taxon);
                    if (id != -1) append_int_array(equivalent_node_ids, id);

#ifndef NDEBUG
                    printf("\ttaxon:%s id:%d\n", *taxon == NULL ? "NULL" : *taxon, id);
#endif
                }

#ifndef NDEBUG
                printf("\n\nLCA search:\n---- ---- ---- ---- ---- ---- ---- ----\n");
#endif

                // now we need to find lowest common ancestor for these nodes
                lca_idx = *(equivalent_node_ids->array);
                for (ip = (equivalent_node_ids->array + 1); ip != equivalent_node_ids->last; ++ip) {
                    lca_idx = lca(lca_idx, *ip);

#ifndef NDEBUG
                    printf("\tLCA(%d, %d): %d\n", lca_idx, *ip, lca_idx);
#endif
                }

                // let's find tau interval for a given lowest common ancestor
                // and to do that we need to calculate it's distance from the farthest leaf from the root
                float lca_dist = farthest_leaf_dist - get_distance_from_root((species_indexed_nodes->array + lca_idx)->node);
                for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
                    if (lca_dist >= *(fp + 1) && lca_dist < *fp) {
                        break;
                    }
                }

                int tau_idx = fp - spec_dists->array;

#ifndef NDEBUG
                printf("---- ---- ---- ---- ---- ---- ---- ----\n\tOverall LCA id: %d with Tau index %d\n", lca_idx, tau_idx);
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
        printf("\tg[%d]:%d\n", ip - g.array, *ip);
    }
#endif

    // alrighty-roo, lets' bead species tree!
    bead_tree(species_tree, spec_dists, farthest_leaf_dist);

    newick_node_ptr_array_array Y;
    init_newick_node_ptr_array_array(&Y);
    for (fp = spec_dists->array; fp != (spec_dists->last - 1); ++fp) {
        //int_array cur_interval_nodes;
        //init_int_array(&cur_interval_nodes);
        newick_node_ptr_array cur_interval_nodes;
        init_newick_node_ptr_array(&cur_interval_nodes);

        // species_tree is a beaded tree in here
        printf("add nodes in interval [%f, %f)", *(fp + 1), *fp); puts("");
        add_nodes_in_interval(&cur_interval_nodes, species_tree, *(fp + 1), *fp, farthest_leaf_dist);

        // get nodes for the current tau
        append_newick_node_ptr_array_array(&Y, cur_interval_nodes);
    }
/*
    // let's allocate k[][][] array
    int k_z_dim = 5, k_i_dim = 5, k_j_dim = 5;
    int ***k_arr;
    k_arr = (int ***) malloc(sizeof(int**)*k_i_dim);
    for (i = k_i_dim; i < k_i_dim; ++i) {
        *(k_arr + i) = (int**)malloc(sizeof(int*)*k_j_dim);
        for (j = 0; j < k_j_dim; ++j) {
            *(*(k_arr + i) + j) = (int *)malloc(sizeof(int)*k_z_dim);
        }
    }
*/

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

    for (i = 2; i < speciation_count; ++i) {
        for (z = 0; z < array_size(Y.array[i]); ++z) {
            int_array topology_prefix = get_topology_prefix(species_tree, Y.array[i].array[z]);
            K.array[i].array[0].array[z] = get_exit_branches(gene_tree, (spec_dists->array)[i - 1], &topology_prefix, farthest_leaf_dist);
        }
    }

    lca_end();
    monitored_memory_end();
    return 0;
}
