#ifndef __MONITORED_MEMORY_H__
#define __MONITORED_MEMORY_H__

#include <stdlib.h>

void monitored_memory_init();
void *monitored_malloc(size_t size);
void *monitored_calloc(size_t n, size_t size);
void *monitored_realloc(void *p, size_t size);
void monitored_free(void *p);
void monitored_free_all();
void monitored_memory_end();

#endif
