#ifndef __RANKED_GENE_TREES_H__
#define __RANKED_GENE_TREES_H__

#include <assert.h>
#include <mpfr.h>

#include "newick_tree.h"

#define MPFR_PRECISION 1024

#define DEF_ARRAY_DECL(ITEM_TYPE)\
typedef struct { \
    ITEM_TYPE* array; \
    ITEM_TYPE* last; \
    ITEM_TYPE* end; \
} ITEM_TYPE##_array; \
int init_##ITEM_TYPE##_array(ITEM_TYPE##_array *); \
int append_##ITEM_TYPE##_array(ITEM_TYPE##_array*, ITEM_TYPE); \
void clear_##ITEM_TYPE##_array(ITEM_TYPE##_array*); \

#define DEF_ARRAY_IMPL(ITEM_TYPE) \
int init_##ITEM_TYPE##_array(ITEM_TYPE##_array *a) { \
    size_t size = 16; \
    a->array = (ITEM_TYPE *) malloc(size * sizeof(ITEM_TYPE)); \
    assert(a->array != NULL); \
    a->last = a->array; \
    a->end = a->array + size; \
    return a->array != NULL; \
} \
int append_##ITEM_TYPE##_array(ITEM_TYPE##_array* a, ITEM_TYPE val) { \
    /* printf("%12p %12p %12p\n", a->array, a->last, a->end); */ \
    if (a->last == a->end) { \
        size_t sz = a->end - a->array; \
        ITEM_TYPE * realloc_array = (ITEM_TYPE *)realloc(a->array, (sz * 2) * sizeof(ITEM_TYPE)); \
        if (realloc_array != NULL) { \
            a->array = realloc_array; \
        } \
        assert(realloc_array != NULL); \
        a->last = a->array + sz; \
        a->end = a->array + (sz * 2); \
    } \
    *(a->last) = val; \
    ++(a->last); \
    return a->array != NULL; \
} \
void clear_##ITEM_TYPE##_array(ITEM_TYPE##_array* a) { \
    a->last = a->array; \
} \

#define DEF_ARRAY(ITEM_TYPE) \
DEF_ARRAY_DECL(ITEM_TYPE) \
DEF_ARRAY_IMPL(ITEM_TYPE) \

#define clear_array(X) (X).last = (X).array;
#define array_size(X) ((X).last - (X).array)

typedef struct matidx {
    int i, j;
} matidx;

typedef struct node2float {
    newick_node *node;
    float val;
} node2float;

typedef struct node2mpfr_t {
    newick_node *node;
    mpfr_t val;
} node2mpfr_t;

typedef struct node2int {
    newick_node *node;
    int val;
} node2int;

typedef char * char_ptr;
typedef int * int_ptr;
typedef int ** int_ptr_ptr;
typedef int *** int_ptr_ptr_ptr;
typedef newick_node * newick_node_ptr;

DEF_ARRAY_DECL(int);
DEF_ARRAY_DECL(char);
DEF_ARRAY_DECL(float);
DEF_ARRAY_DECL(mpfr_t);

DEF_ARRAY_DECL(newick_node_ptr);
DEF_ARRAY_DECL(newick_node_ptr_array);

DEF_ARRAY_DECL(node2float);
DEF_ARRAY_DECL(node2mpfr_t);

DEF_ARRAY_DECL(node2float_array);
DEF_ARRAY_DECL(node2int);
DEF_ARRAY_DECL(char_ptr);
DEF_ARRAY_DECL(int_array);
DEF_ARRAY_DECL(int_array_array);

int flt_cmp(float *, float *);
int flt_cmp_desc(float *, float *);

#endif
