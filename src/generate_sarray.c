#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>

#include "Newickform.h"

#include "utils.h"
#include "getopt.h"

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

DEF_ARRAY(node2float);
DEF_ARRAY(node2int);
DEF_ARRAY(char);

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
    printf("%f <= %f? ", max_dist_from_root - distance - t->dist, limit);
    if ((max_dist_from_root - distance - t->dist) <= limit) {
        puts("yes");
        return lineages + 1;
    } else puts("no");

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist, max_dist_from_root);
    }

    return lineages;
}
/*
int do_get_gene_lineages(newick_node *t, float limit, float distance) {
    int lineages = 0;
    newick_child *p;

    if (distance >= limit) return lineages + 1;

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist);
    }

    return lineages;
}

float do_get_gene_lineages_bottom(newick_node *t, float limit, int *lineages, float distance_from_the_root, float max_dist_from_root) {
    newick_child *p;

    float dist = 0.f, max_dist = 0.f;
    for (p = t->child; p != NULL; p = p->next) {
        //dist = max(dist, do_get_speciation_distances(p->node, distance + t->dist, speciation_distances));
        dist = do_get_gene_lineages_bottom(p->node, limit, lineages, distance_from_the_root + t->dist, max_dist_from_root);
        max_dist = max(dist, max_dist);
        //printf("checking between %f and %f limit: %f\n", dist, dist + t->dist, limit);
        //if (dist < limit && dist + t->dist > limit) {
        //    ++(*lineages);
        //}
    }

    printf("%f in [%f .. %f)? ", limit, max_dist, max_dist + t->dist);
    if (limit >= max_dist  && limit < max_dist + t->dist) {
        ++(*lineages);
        puts("yes!");
    } else {
        puts("no!");
    }

    return max_dist + t->dist;
}
*/
int_array *get_gene_lineages(float_array *speciation_distances, newick_node *t, float max_dist_from_root) {
    float *p;
    int_array *lineages = (int_array *)malloc(sizeof(int_array));

    init_int_array(lineages);
    for (p = speciation_distances->array; p != speciation_distances->last; ++p) {
        append_int_array(lineages, do_get_gene_lineages(t, *p, 0.f, max_dist_from_root));
    }
    return lineages;
}

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
        printf("iterating thru childs of u[%d]@%d:\n", v, coalescence_array->array[v].node);
        // find child's index (to)
        for (cp = coalescence_array->array, to = 0; cp != coalescence_array->last; ++to, ++cp) {
            //if (cp->node == coalescence_array->array[v].node)
            if (cp->node == np->node)
                break;
        }

        printf("u[%d](u[%d])@%d is a child? ", to, cp->val, np->node);
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

void lca_init(int n, newick_node *t, node2int_array *coalescence_array) {
    int i;

    tin = (int *) malloc(n * sizeof(int));
    tout = (int *) malloc(n * sizeof(int));
    up = (int **) malloc(n * sizeof(int *));

    l = 1;
    while ((1 << l) <= n)  ++l;
    for (i = 0; i < n; ++i)  up[i] = (int *) malloc((l + 1) * sizeof(int));

    lca_preprocess(coalescence_array, 0, 0);

    // queries:
    //for (;;) {
    //    int a, b; for u[a], u[b]
    //    int res = lca (a, b);
    //}
}

void lca_end() {
    int i;
    free(tin);
    free(tout);
    for (i = 0; i < n; ++i) free(up[i]);
    free(up);
}

struct globalArgs_t {
    const char *inFileName;
    const char *outFileName;
} globalArgs;

static const char *optString = "i:o:vh";

static const struct option longOpts[] = {
    { "in", required_argument, NULL, 'i' },
    { "out", required_argument, NULL, 'o' },
    { "version", no_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 } /* long options with no short equivalent go to "case 0:" */
};

void init_global_args() {
    globalArgs.inFileName = NULL;
    globalArgs.outFileName = NULL;
}

void version(char *prog_name) {
    printf("%s, v0.1 developed by Islam Ismailov, with support of NESCENT as part of the Google Summer of Code 2012,\nMentors: J Degnan, T Stadler\n", prog_name);
}

void usage(char *prog_name) {
    printf("Usage: %s -i (--in) (newick tree file) -o (--out) outputfile\n", prog_name);
}

int main(int argc, char **argv) {
    newick_node *root = NULL;
    char_array *tree_string = (char_array*) malloc(sizeof(char_array));
    float_array *spec_dists;
    int_array *gene_lineages;
    int coalescence_count;
    node2int_array *coalescence_array;
    int *ip;
    float *fp;
    node2int *n2i;
    float farthest_leaf_dist = 0.f;

    FILE *f;
    int i;

    char c;
    int longIndex;
    init_global_args();

    int opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    while (opt != -1) {
        switch (opt) {
            case 'i':
                globalArgs.inFileName = optarg;
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

    if (globalArgs.inFileName == NULL || globalArgs.outFileName == NULL) {
        usage(argv[0]);
        return -1;
    }

    // Open tree file
    f = fopen(globalArgs.inFileName, "r+");

    // Read in tree string
    init_char_array(tree_string);
    printf("\n\n\n");
    for (c = fgetc(f); c != EOF && c != '\n'; c = fgetc(f)) {
        putchar(c);
        append_char_array(tree_string, c);
    }
    printf("\n\n\n");
    append_char_array(tree_string, '\0');
    fclose(f);

    printf("tree string:\n%s\n\n", tree_string->array );

    // Parse tree string
    root = parseTree(tree_string->array);
    printTree(root);
    printf("\n\n");

    farthest_leaf_dist = max_dist_from_root(root);

    printf("max dist from root: %f\n", farthest_leaf_dist);

    // let's do smth interesting now
    spec_dists = get_speciation_distances(root, farthest_leaf_dist);
    printf("Speciation distances:\n----\n");
    for (fp = spec_dists->array; fp != spec_dists->last; ++fp) {
        printf("%f\n", *fp);
    }
    printf("----\n");

    gene_lineages = get_gene_lineages(spec_dists, root, farthest_leaf_dist);
    printf("Gene Lineages:\n");
    for (ip = gene_lineages->array; ip != gene_lineages->last; ++ip) {
        printf("%d\n", *ip);
    }

    coalescence_count = get_tree_coalescence_count(root);
    coalescence_array = get_coalescence_array(root);
    printf("Coalescence array (size %d):\n---- ---- ---- ----\n", coalescence_count);
    for (n2i = coalescence_array->array; n2i != coalescence_array->last; ++n2i) {
        printf("val:%d node:%d childnum:%d\n", n2i->val, n2i->node, n2i->node->childNum);
    }
    printf("---- ---- ---- ---\n");
    lca_init(coalescence_count,root, coalescence_array);

    printf("lca of %d and %d is %d\n", 2, 3, lca (2, 3));

    lca_end();

    free(tree_string->array);

    return 0;
}
