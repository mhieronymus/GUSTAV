#ifndef _SPLINE_TABLE_H
#define _SPLINE_TABLE_H

#include <stdlib.h>
#include "../types.hpp"

#ifdef __cplusplus
extern "C" {
#endif

struct splinetable_buffer {
	void *data;
	size_t size;
	void *(*mem_alloc)(size_t newsize);
	void *(*mem_realloc)(void *p, size_t newsize);
};

typedef enum {
	SPLINETABLE_index_t,
	SPLINETABLE_value_t
} splinetable_dtype;

index_t readsplinefitstable(const char *path, struct splinetable *table);
index_t readsplinefitstable_mem(struct splinetable_buffer *buffer,
    struct splinetable *table);
index_t writesplinefitstable(const char *path, const struct splinetable *table);
index_t writesplinefitstable_mem(struct splinetable_buffer *buffer,
    const struct splinetable *table);
void splinetable_free(struct splinetable *table);
void splinetable_permute(struct splinetable *table, index_t *permutation);
char * splinetable_get_key(const struct splinetable *table, const char *key);
index_t splinetable_read_key(const struct splinetable *table, splinetable_dtype type,
    const char *key, void *result);

#ifdef __cplusplus
}
#endif

#endif /* _SPLINE_TABLE_H */

