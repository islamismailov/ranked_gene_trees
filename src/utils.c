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

DEF_ARRAY_IMPL(int);
DEF_ARRAY_IMPL(char);
DEF_ARRAY_IMPL(float);

DEF_ARRAY_IMPL(newick_node_ptr);
DEF_ARRAY_IMPL(newick_node_ptr_array);

DEF_ARRAY_IMPL(node2float);
DEF_ARRAY_IMPL(node2float_array);

DEF_ARRAY_IMPL(node2int);
DEF_ARRAY_IMPL(char_ptr);
DEF_ARRAY_IMPL(int_array);
DEF_ARRAY_IMPL(int_array_array);

/* mpfr array stuff */
int init_mpfr_t_array(mpfr_t_array *a) {
    size_t size = 16;
    a->array = (mpfr_t *) malloc(size * sizeof(mpfr_t));
    a->last = a->array;
    a->end = a->array + size;
    return a->array != NULL;
}

int append_mpfr_t_array(mpfr_t_array *a, mpfr_t val) {
    if (a->last == a->end) {
        size_t sz = a->end - a->array;
        mpfr_t *realloc_array = (mpfr_t *)realloc(a->array, (sz * 2) * sizeof(mpfr_t));
        if (a->array != NULL) {
            a->array = realloc_array;
        }
        a->last = a->array + sz;
        a->end = a->array + (sz * 2);
    }
    //*(a->last) = val;
    mpfr_init2(*(a->last), MPFR_RNDN);
    mpfr_set(*(a->last), val, MPFR_RNDN);
    
    ++(a->last);
    return a->array != NULL;
}

int init_node2mpfr_t_array(node2mpfr_t_array *a) {
    size_t size = 16;
    a->array = (node2mpfr_t *) malloc(size * sizeof(node2mpfr_t));
    a->last = a->array;
    a->end = a->array + size;
    return a->array != NULL;
}

int append_node2mpfr_t_array(node2mpfr_t_array *a, node2mpfr_t val) {
    if (a->last == a->end) {
        size_t sz = a->end - a->array;
        node2mpfr_t *realloc_array = (node2mpfr_t *)realloc(a->array, (sz * 2) * sizeof(node2mpfr_t));
        if (a->array != NULL) {
            a->array = realloc_array;
        }
        a->last = a->array + sz;
        a->end = a->array + (sz * 2);
    }
    //*(a->last) = val;
    (*(a->last)).node = val.node;
    mpfr_init2((*(a->last)).val, MPFR_RNDN);
    mpfr_set((*(a->last)).val, val.val, MPFR_RNDN);
    
    ++(a->last);
    return a->array != NULL;
}
/* mpfr array stuff */




/*
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
        if (a->array != NULL) {
            a->array = realloc_array;
        }
        a->last = a->array + sz;
        a->end = a->array + (sz * 2);
    }
    *(a->last) = val;
    ++(a->last);
    return a->array != NULL;
}
*/
