#ifndef __RANKED_GENE_TREES_H__
#define __RANKED_GENE_TREES_H__

typedef struct float_array {
    float *array;
    float *end;
    float *last;
} float_array;

int init_float_array(float_array *);

int append_float_array(float_array *, float);

typedef struct int_array {
    int *array;
    int *end;
    int *last;
} int_array;

int init_int_array(int_array *);

int append_int_array(int_array *, int);

int flt_cmp(float *, float *);

int flt_cmp_desc(float *, float *);

#define DEF_ARRAY(ITEM_TYPE) \
   typedef struct { \
      ITEM_TYPE* array; \
      ITEM_TYPE* last; \
      ITEM_TYPE* end; \
   } ITEM_TYPE##_array; \
   int init_##ITEM_TYPE##_array(ITEM_TYPE##_array *a) { \
       size_t size = 16; \
       a->array = (ITEM_TYPE *) malloc(size * sizeof(ITEM_TYPE)); \
       a->last = a->array; \
       a->end = a->array + size; \
       return a->array != NULL; \
   } \
   int append_##ITEM_TYPE##_array(ITEM_TYPE##_array* a, ITEM_TYPE val) \
   { \
       if (a->last != a->end) { \
           *(a->last) = val; \
           ++(a->last); \
       } else { \
           size_t sz = a->end - a->array; \
           a->array = (ITEM_TYPE *)realloc(a->array, (sz * 2) * sizeof(ITEM_TYPE)); \
           a->last = a->end; \
           a->end = a->array + (sz * 2); \
       } \
       return a->array != NULL; \
   } \

#endif
