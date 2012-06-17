#include <stdio.h>
#include "hash_table.h"

unsigned long long HASH_PRIME = 1181783497276652981LL; // if c == 0      x_n = (h*x_n-1 +c) % 2^64
//unsigned long long HASH_PRIME = 2862933555777941757LL; // if c is odd in x_n = (h*x_n-1 +c) % 2^64
int HASH_PRIMES[] = { 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739 };

unsigned long long hash(void *p, size_t len_bytes) {
    char *byte_ptr;
    unsigned long long h = 0;

    for (byte_ptr = (char *)p; byte_ptr < (char *)p + len_bytes; ++byte_ptr) {
        h *= HASH_PRIME;
        h += *byte_ptr;
    }

    return h;
}

hash_table *get_new_hash_table(/*size_t capacity*/) {
    hash_table *new_table = (hash_table *) malloc(sizeof(hash_table));
    new_table->capacity_idx = 0;

    // calloc nulls array for us
    new_table->elements = (hash_tab_element_ptr *) calloc(HASH_PRIMES[new_table->capacity_idx], sizeof(hash_tab_element_ptr));
    //new_table->tab = (hash_tab_element *) malloc(HASH_PRIMES[new_table->capacity_idx] * sizeof(hash_tab_element));

    return new_table;
}

void grow(hash_table *table) {
    hash_tab_element_ptr *old_elements = table->elements;
    ++table->capacity_idx;
    table->elements = (hash_tab_element_ptr *) calloc(HASH_PRIMES[table->capacity_idx], sizeof(hash_tab_element_ptr));

    int i;
    hash_tab_element *iter, *iter_next;
    for (i = 0; i < HASH_PRIMES[table->capacity_idx - 1]; ++i) {
        iter = old_elements[i];
        while (iter != NULL) {
            do_insert(table, iter->ptr, iter->h_val);
            iter_next = iter->next;
            free(iter);
            iter = iter_next;
        }
    }
    free(old_elements);
    old_elements = NULL;
}

void do_insert(hash_table *table, void *p, unsigned long long h_val) {
    int pos = h_val % HASH_PRIMES[table->capacity_idx];
    hash_tab_element *new_element = (hash_tab_element *) malloc(sizeof(hash_tab_element));
    new_element->ptr = p;
    new_element->h_val = h_val;

    if (table->elements[pos] == NULL) {
        new_element->next = NULL;
        table->elements[pos] = new_element;
    } else {
        new_element->next = table->elements[pos];
        table->elements[pos] = new_element;
    }
    ++table->elements_inserted;
}

void insert(hash_table *table, void *p, size_t len_bytes) {
    if (table->elements_inserted >= HASH_PRIMES[table->capacity_idx] * 0.95) {
        grow(table);
    }
    do_insert(table, p, hash(p, len_bytes));
}

void *lookup(hash_table *table, void *p, size_t len_bytes, __compar_fn_t __compar) {
    hash_tab_element *iter;
    unsigned long long h_val = hash(p, len_bytes);
    int pos = h_val % HASH_PRIMES[table->capacity_idx];
    hash_tab_element *chain_head = table->elements[pos];

    for (iter = chain_head; iter != NULL; iter = iter->next) {
        if (__compar(p, iter->ptr) == 0) {
            return iter->ptr;
        }
    }
    return NULL;
}
