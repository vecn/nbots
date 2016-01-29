#ifndef __VCN_TEST_H__
#define __VCN_TEST_H__

#include <stdbool.h>

typedef struct {
	void* (*init)();
	bool (*run)();
	void (*clear)(void*);
	char *label;
} test_t;

test_t* test_create(void);
double test_get_avg_seconds(const test_t *test, void *init_data);

#endif
