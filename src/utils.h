#ifndef __RANKED_GENE_TREES_H__
#define __RANKED_GENE_TREES_H__

#include <assert.h>

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

#define clear_array(X) (X).last = (X).array;
#define array_size(X) ((X).last - (X).array)

#endif
