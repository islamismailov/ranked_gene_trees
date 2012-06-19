#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

int flt_cmp(float *a, float *b) {
    float diff = *a - *b;

    if (diff > 1e-6)        return  1;
    else if (diff < -1e-6)  return -1;

    return 0;
}

int flt_cmp_desc(float *a, float *b) {
    float diff = *a - *b;

    if (diff > 1e-6)        return -1;
    else if (diff < -1e-6)  return  1;

    return 0;
}

int init_float_array(float_array *a) {
    size_t size = 16;
    a->array = (float *) malloc(size * sizeof(float));
    a->last = a->array;
    a->end = a->array + size;
    return a->array != NULL;
}

int append_float_array(float_array *a, float val) {
    if (a->last == a->end) {
        size_t sz = a->end - a->array;
        float *realloc_array = (float *)realloc(a->array, (sz * 2) * sizeof(float));
        if (a->array != realloc_array) {
            free(a->array);
            a->array = realloc_array;
        }
        a->last = a->array + sz;
        a->end = a->array + (sz * 2);
    }
    *(a->last) = val;
    ++(a->last);
    return a->array != NULL;
}

void clear_float_array(float_array *a) {
    a->last = a->array;
}

int init_int_array(int_array *a) {
    size_t size = 16;
    a->array = (int *) malloc(size * sizeof(int));
    a->last = a->array;
    a->end = a->array + size;
    return a->array != NULL;
}

int append_int_array(int_array *a, int val) {
    if (a->last == a->end) {
        size_t sz = a->end - a->array;
        int *realloc_array = (int *)realloc(a->array, (sz * 2) * sizeof(int));
        if (a->array != realloc_array) {
            free(a->array);
            a->array = realloc_array;
        }
        a->last = a->array + sz;
        a->end = a->array + (sz * 2);
    }
    *(a->last) = val;
    ++(a->last);
    return a->array != NULL;
}

void clear_int_array(int_array *a) {
    a->last = a->array;
}
