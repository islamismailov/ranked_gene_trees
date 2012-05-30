#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include "Newickform.h"

#include "utils.h"

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

void do_get_speciation_distances(newick_node *t, float distance, float_array *speciation_distances) {
    newick_child *p;
    append_float_array(speciation_distances, distance + t->dist);
    for (p = t->child; p != NULL; p = p->next) {
        do_get_speciation_distances(p->node, distance + t->dist, speciation_distances);
    }
}

float_array *get_speciation_distances(newick_node *t) {
    float_array *speciation_distances = (float_array *)malloc(sizeof(float_array));
    init_float_array(speciation_distances);
    do_get_speciation_distances(t, 0.0f, speciation_distances);
    append_float_array(speciation_distances, FLT_MAX);
    qsort(speciation_distances->array, speciation_distances->last - speciation_distances->array,
        sizeof(float), (int(*)(const void*,const void*))flt_cmp);
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
    append_tree_hash_dist_array(dist_array, pair);
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
        count += get_tree_node_count(p->node);
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

/*
    for (size_t i = 0; i < g[v].size(); ++i) {
        int to = g[v][i];
        if (to != p)
            lca_preprocess (to, v);
    }
*/
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
