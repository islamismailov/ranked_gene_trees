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

int node2float_cmp(node2float a, node2float b) {
    float diff = a.val - b.val;

    if (diff > 1e-6)        return  1;
    else if (diff < -1e-6)  return -1;

    return 0;
}

DEF_ARRAY(node2float);
DEF_ARRAY(node2int);
DEF_ARRAY(char);

void do_get_speciation_distances(newick_node *t, float distance, float_array *speciation_distances) {
    newick_child *p;
    if (t->childNum > 0) // speciation happened
        append_float_array(speciation_distances, distance + t->dist);
    for (p = t->child; p != NULL; p = p->next) {
        do_get_speciation_distances(p->node, distance + t->dist, speciation_distances);
    }
}

float do_get_speciation_distances_bottom(newick_node *t, float distance, float_array *speciation_distances) {
    newick_child *p;

    //if (t->childNum == 0) return t->dist;

    float dist = 0.f;
    for (p = t->child; p != NULL; p = p->next) {
        dist = max(dist, do_get_speciation_distances_bottom(p->node, distance + t->dist, speciation_distances));
    }

    //dist += t->dist;

    if (t->childNum > 0) // speciation happened
        append_float_array(speciation_distances, dist);

    return dist + t->dist;
}

float_array *get_speciation_distances(newick_node *t) {
    float_array *speciation_distances = (float_array *)malloc(sizeof(float_array));
    init_float_array(speciation_distances);
    //do_get_speciation_distances(t, 0.0f, speciation_distances);
    do_get_speciation_distances_bottom(t, 0.0f, speciation_distances);
    append_float_array(speciation_distances, FLT_MAX);

    qsort(speciation_distances->array, speciation_distances->last - speciation_distances->array,
        sizeof(float), (int(*)(const void*,const void*))flt_cmp_desc);
    return speciation_distances;
}

int do_get_gene_lineages(newick_node *t, float limit, float distance) {
    int lineages = 0;
    newick_child *p;

    if (distance >= distance) return lineages + 1;

    for (p = t->child; p != NULL; p = p->next) {
        lineages += do_get_gene_lineages(p->node, limit, distance + t->dist);
    }

    return lineages;
}

int_array *get_gene_lineages(float_array *speciation_distances, newick_node *t) {
    float *p;
    int_array *lineages = (int_array *)malloc(sizeof(int_array));

    init_int_array(lineages);
    for (p = speciation_distances->array; p != speciation_distances->last; ++p) {
        append_int_array(lineages, do_get_gene_lineages(t, *p, 0.f));
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
    int i;
    node2float_array *dist_array = (node2float_array *)malloc(sizeof(node2float_array));
    node2int_array *coalescence_array = (node2int_array *)malloc(sizeof(node2int_array));

    do_get_distance_array(t, 0.f, dist_array);
    qsort(dist_array->array, dist_array->last - dist_array->array,
        sizeof(node2float), (int(*)(const void*,const void*))node2float_cmp);

    int coalescence_count = get_tree_coalescence_count(t);

    for (i = 0; i < coalescence_count; ++i) {
        node2int pair;
        pair.node = (dist_array->array + i)->node;
        pair.val = i;

        append_node2int_array(coalescence_array, pair);
    }

    return coalescence_array;
}












int timer, n, l;
int *tin, *tout;
int **up;

void lca_preprocess(node2int_array *coalescence_array, int v, int p) {
    int i;
    node2int *cp;
    newick_child *np;

    tin[v] = ++timer;
    up[v][0] = p;
    for (i = 1; i <= l; ++i)
        up[v][i] = up[up[v][i - 1]][i - 1];

    // iterate through children of coalescence u[v]:
    for (np = coalescence_array->array[v].node->child; np != NULL; np = np->next) {
        int to = 0;
        for (cp = coalescence_array->array; cp != coalescence_array->last; ++cp) {
            if (cp->node == coalescence_array->array[v].node)
                break;
            ++to;
        }
        if (to != p)
            lca_preprocess(coalescence_array, to, v);
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
    float *fp;

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
    for (c = fgetc(f); c != EOF && c != '\n'; c = fgetc(f)) {
        append_char_array(tree_string, c);
    }
    append_char_array(tree_string, '\0');
    fclose(f);

    // Parse tree string
    root = parseTree(tree_string->array);
    printTree(root);
    printf("\n");

    // let's do smth interesting now
    spec_dists = get_speciation_distances(root);
    printf("Speciation distances:\n");
    for (fp = spec_dists->array; fp != spec_dists->last; ++fp) {
        printf("dist: %f\n", *fp);
    }
    //gene_lineages = get_gene_lineages(spec_dists, root);
    //coalescence_count = get_tree_coalescence_count(root);
    //coalescence_array = get_coalescence_array(root);

    free(tree_string->array);

    return 0;
}