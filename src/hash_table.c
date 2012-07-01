#include <stdio.h>
#include "hash_table.h"

int compar_int(int *a, int *b) {
    return (*a > *b) - (*a < *b);
}

int compar_addr(const void *a, const void *b) {
    return a - b;
}

const hash_t HASH_PRIME = 1181783497276652981LL; // if c == 0      x_n = (h*x_n-1 +c) % 2^64
//const hash_t HASH_PRIME = 2862933555777941757LL; // if c is odd in x_n = (h*x_n-1 +c) % 2^64

const size_t HASH_PRIMES[] = {
         23,       47,       97,      193,       389,       769,      1543,      3079,       6151,
      12289,    24593,    49157,    98317,    196613,    393241,    786433,   1572869,    3145739,
    6291487, 12582929, 25165843, 50331683, 100663291, 201326557, 402653213, 805306339, 1610612747
};

hash_t htab_hash(void *p, size_t len_bytes) {
    char *byte_ptr;
    hash_t h = 0;

    for (byte_ptr = (char *)p; byte_ptr < (char *)p + len_bytes; ++byte_ptr) {
        h *= HASH_PRIME;
        h += *byte_ptr;
    }
    return h;
}

hash_table *htab_get_new(/*size_t capacity*/) {
    hash_table *new_table = (hash_table *) malloc(sizeof(hash_table));
    new_table->capacity_idx = 0;

    // calloc nulls array for us
    new_table->elements = (hash_tab_element_ptr *) calloc(HASH_PRIMES[new_table->capacity_idx], sizeof(hash_tab_element_ptr));
    return new_table;
}

void htab_grow(hash_table *table) {
    hash_tab_element_ptr *old_elements = table->elements;
    ++table->capacity_idx;
    // calloc nulls array for us
    table->elements = (hash_tab_element_ptr *) calloc(HASH_PRIMES[table->capacity_idx], sizeof(hash_tab_element_ptr));

    size_t i;
    hash_tab_element *iter, *iter_next;
    for (i = 0; i < HASH_PRIMES[table->capacity_idx - 1]; ++i) {
        iter = old_elements[i];
        while (iter != NULL) {
            htab_do_insert(table, iter->key, iter->hash_key, iter->val);
            iter_next = iter->next;
            free(iter);
            iter = iter_next;
        }
        old_elements[i] = NULL;
    }
    free(old_elements);
    old_elements = NULL;
}

void htab_do_insert(hash_table *table, void *key, hash_t hash_key, void *val) {
    size_t pos = hash_key % HASH_PRIMES[table->capacity_idx];
    hash_tab_element *new_element = (hash_tab_element *) malloc(sizeof(hash_tab_element));
    new_element->key = key;
    new_element->val = val;
    new_element->hash_key = hash_key;

    if (table->elements[pos] == NULL) {
        new_element->next = NULL;
        table->elements[pos] = new_element;
    } else {
        new_element->next = table->elements[pos];
        table->elements[pos] = new_element;
    }
    ++table->elements_inserted;
}

void htab_insert(hash_table *table, void *key, size_t len_bytes, void *val) {
    if (table->elements_inserted >= HASH_PRIMES[table->capacity_idx] * 0.95) {
        htab_grow(table);
    }
    htab_do_insert(table, key, htab_hash(key, len_bytes), val);
}

void *htab_do_lookup(hash_table *table, void *key, hash_t hash_key, __compar_fn_t __compar) {
    hash_tab_element *iter;
    size_t pos = hash_key % HASH_PRIMES[table->capacity_idx];
    
    for (iter = table->elements[pos]; iter != NULL; iter = iter->next) {
        if (__compar(key, iter->key) == 0) {
            return iter->val;
        }
    }
    return NULL;    
}

void *htab_lookup(hash_table *table, void *key, size_t len_bytes, __compar_fn_t __compar) {
    return htab_do_lookup(table, key, htab_hash(key, len_bytes), __compar);
}

void htab_do_remove(hash_table *table, void *key, hash_t hash_key, __compar_fn_t __compar) {
    size_t pos = hash_key % HASH_PRIMES[table->capacity_idx];

    hash_tab_element *next_iter = NULL, *prev_iter = NULL, *iter = table->elements[pos];
    while (iter != NULL) {
        if (__compar(key, iter->key) == 0) {
            if (prev_iter == NULL) {
                table->elements[pos] = NULL;
            } else {
                prev_iter->next = iter->next;
            }
            next_iter = iter->next;
            free(iter);
            iter = next_iter;
        } else {
            prev_iter = iter;
            iter = iter->next;
        }
    }
}

void htab_remove(hash_table *table, void *key, size_t len_bytes, __compar_fn_t __compar) {
    htab_do_remove(table, key, htab_hash(key, len_bytes), __compar);
}

void htab_free_table(hash_table *table) {
    size_t i;
    hash_tab_element *iter, *iter_next;
    for (i = 0; i < HASH_PRIMES[table->capacity_idx]; ++i) {
        iter = table->elements[i];
        while (iter != NULL) {
            //htab_do_insert(table, iter->ptr, iter->h_val);
            iter_next = iter->next;
            free(iter);
            iter = iter_next;
        }
        table->elements[i] = NULL;
    }
    free(table->elements);
    table->capacity_idx = 0;
    table->elements_inserted = 0;
}
