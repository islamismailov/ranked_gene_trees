#include "Newickform.h"

typedef struct float_array {
    float *array;
    float *end;
    float *last;
} float_array;

init_float_array(float_array *a) {
    size_t size = 16;
    a->array = (float *) malloc(size * sizeof(float));
    a->last = a->array;
    a->end = a->array + size;
}

void append_float_array(float_array *a, float f) {
    if (a->last != a->end) {
        *(a->last) = f;
        ++(a->last);
    } else {
        size_t sz = a->end - a->array;
        a->array = (float *)realloc(a->array, (sz * 2) * sizeof(float));
        a->end = a->array + (sz * 2);
        a->last = a->array + sz;
    }
}

int flt_cmp(float a, float b) {
    float diff = a - b;

    if (diff > 1e-6)        return  1;
    else if (diff < -1e-6)  return -1;

    return 0;
}

void do_dfs(newick_node *t, float distance, float_array *speciation_distances) {
    append_float_array(speciation_distances, distance + t->dist);
    for (newick_child *p = t->child; p != NULL; p = p->next) {
        dfs(p, distance + t->dist, speciation_distances);
    }
}

float_array *get_speciation_distances(newick_node *t) {
    float_array *speciation_distances = (float_array *)malloc(sizeof(float_array));
    init_float_array(speciation_distances);
    do_dfs(t, 0.0f, speciation_distances);
    qsort(speciation_distances->array, speciation_distances->last - speciation_distances->array, sizeof(float), flt_cmp);
    return speciation_distances;
}
