#ifndef __RANKED_GENE_TREES_H__
#define __RANKED_GENE_TREES_H__

#include <assert.h>

typedef struct float_array {
    float *array;
    float *end;
    float *last;
} float_array;

int init_float_array(float_array *);
int append_float_array(float_array *, float);
void clear_float_array(float_array *);

typedef struct int_array {
    int *array;
    int *end;
    int *last;
} int_array;

int init_int_array(int_array *);
int append_int_array(int_array *, int);
void clear_int_array(int_array *);

int flt_cmp(float *, float *);
int flt_cmp_desc(float *, float *);

#define DEF_ARRAY(ITEM_TYPE) \
   typedef struct { \
      ITEM_TYPE* array; \
      ITEM_TYPE* last; \
      ITEM_TYPE* end; \
   } ITEM_TYPE##_array; \
   \
   int init_##ITEM_TYPE##_array(ITEM_TYPE##_array *a) { \
       size_t size = 16; \
       a->array = (ITEM_TYPE *) malloc(size * sizeof(ITEM_TYPE)); \
       assert(a->array != NULL); \
       a->last = a->array; \
       a->end = a->array + size; \
       return a->array != NULL; \
   } \
   int append_##ITEM_TYPE##_array(ITEM_TYPE##_array* a, ITEM_TYPE val) { \
       /* printf("%12u %12u %12u\n", a->array, a->last, a->end); */ \
       if (a->last == a->end) { \
           size_t sz = a->end - a->array; \
           ITEM_TYPE *realloc_array = (ITEM_TYPE *)realloc(a->array, (sz * 2) * sizeof(ITEM_TYPE)); \
           if (a->array != realloc_array) { \
               free(a->array); \
               a->array = realloc_array; \
           } \
           assert(a->array != NULL); \
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

#define clear_array(X) X.last = X.array;

#endif
