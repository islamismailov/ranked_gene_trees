#include <stdlib.h>

#include "monitored_memory.h"
#include "hash_table.h"

hash_table *mem_tab;

void monitored_memory_init() {
    mem_tab = htab_get_new();
}

void *monitored_malloc(size_t size) {
    void *p = malloc(size);
    htab_do_insert(mem_tab, p, (size_t)p, p); // use address as a hash value
    return p;
}

void *monitored_calloc(size_t n, size_t size) {
    void *p = calloc(n, size);
    htab_do_insert(mem_tab, p, (size_t)p, p); // use address as a hash value
    return p;
}

void *monitored_realloc(void *p, size_t size) {
    void *realloc_p = realloc (p, size);
    if (realloc_p != p) {
        htab_do_remove(mem_tab, p, (size_t) p, (__compar_fn_t) compar_addr);
        htab_do_insert(mem_tab, realloc_p, (size_t)realloc_p, p);  // use address as a hash value
    }
    return realloc_p;
}

void monitored_free(void *p) {
    htab_do_remove(mem_tab, p, (size_t)p, (__compar_fn_t)compar_int);  // use address as a hash value
}

void monitored_free_all() { // free memory for each element in hash table, and remove that element (mems + free hash table)
    size_t i;
    hash_tab_element *iter, *iter_next;
    for (i = 0; i < HASH_PRIMES[mem_tab->capacity_idx]; ++i) {
        iter = mem_tab->elements[i];
        while (iter != NULL) {
            htab_do_insert(mem_tab, iter->key, iter->hash_key, iter->val);
            iter_next = iter->next;
            free(iter->val);
            free(iter);
            iter = iter_next;
        }
        mem_tab->elements[i] = NULL;
    }
    mem_tab->capacity_idx = 0;
    mem_tab->elements_inserted = 0;
}

void monitored_memory_end() {
    monitored_free_all();
    htab_free_table(mem_tab);
}
