#ifndef __RANKED_GENE_TREES_HASH_TABLE_H__
#define __RANKED_GENE_TREES_HASH_TABLE_H__

#include <stdlib.h>

int compar_int(int *a, int *b);
int compar_addr(const void *a, const void *b);

typedef unsigned long long hash_t;

#ifndef __COMPAR_FN_T
#define __COMPAR_FN_T
typedef int(*   __compar_fn_t )(const void *, const void *);
#endif

extern hash_t HASH_PRIME;

extern size_t HASH_PRIMES[];

typedef struct hash_tab_element {
    void *key;
    void *val;
    hash_t hash_key;
    struct hash_tab_element *next;
} hash_tab_element;

typedef hash_tab_element* hash_tab_element_ptr;

typedef struct hash_table {
    hash_tab_element_ptr *elements;
    int capacity_idx;
    int elements_inserted;
} hash_table;

unsigned long long htab_hash(void *p, size_t len_bytes);

hash_table *htab_get_new(/*size_t capacity*/);

void htab_do_insert(hash_table *table, void *key, hash_t hash_key, void *val);
void htab_insert(hash_table *table, void *key, size_t len_bytes, void *val);

void *htab_do_lookup(hash_table *table, void *key, hash_t hash_key, __compar_fn_t __compar);
void *htab_lookup(hash_table *table, void *key, size_t len_bytes, __compar_fn_t __compar);

void htab_do_remove(hash_table *table, void *key, hash_t hash_key, __compar_fn_t __compar);
void htab_remove(hash_table *table, void *key, size_t len_bytes, __compar_fn_t __compar);

void htab_free_table(hash_table *hash);

#endif
