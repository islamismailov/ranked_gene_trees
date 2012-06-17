#ifndef __RANKED_GENE_TREES_HASH_TABLE_H__
#define __RANKED_GENE_TREES_HASH_TABLE_H__

#include <stdlib.h>

extern unsigned long long HASH_PRIME;

extern int HASH_PRIMES[];

typedef struct hash_tab_element {
    void *ptr;
    unsigned long long h_val;
    struct hash_tab_element *next;
} hash_tab_element;

typedef hash_tab_element* hash_tab_element_ptr;

typedef struct hash_table {
    hash_tab_element_ptr *elements;
    int capacity_idx;
    int elements_inserted;
} hash_table;

unsigned long long hash(void *p, size_t len_bytes);
hash_table *get_new_hash_table(/*size_t capacity*/);
void insert(hash_table *hash, void *p, size_t len_bytes);
void do_insert(hash_table *table, void *p, unsigned long long h_val);
void *lookup(hash_table *table, void *p, size_t len_bytes, __compar_fn_t __compar);

#endif
