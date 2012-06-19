#ifndef __RANKED_GENE_TREES_HASH_TABLE_H__
#define __RANKED_GENE_TREES_HASH_TABLE_H__

#include <stdlib.h>

int compar_int(int *a, int *b);

typedef unsigned long long hash_t;

extern unsigned long long HASH_PRIME;

extern size_t HASH_PRIMES[];

typedef struct hash_tab_element {
    void *ptr; // void *val;
    void *key;
    unsigned long long h_val;
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
void htab_insert(hash_table *hash, void *p, size_t len_bytes);
void htab_insert_ex(hash_table *hash, void *p, void *k, size_t k_len_bytes);
void htab_do_insert(hash_table *table, void *p, unsigned long long h_val);
void *htab_lookup(hash_table *table, void *p, size_t len_bytes, __compar_fn_t __compar);
void htab_do_remove(hash_table *table, void *p, hash_t h_val, __compar_fn_t __compar);
void htab_remove(hash_table *table, void *p, size_t len_bytes, __compar_fn_t __compar);
void htab_free_table(hash_table *hash);

#endif
