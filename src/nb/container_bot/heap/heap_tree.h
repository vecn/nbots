#ifndef __NB_NB_CONTAINER_BOT_HEAP_HEAP_TREE_H__
#define __NB_NB_CONTAINER_BOT_HEAP_HEAP_TREE_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/memory_bot.h"

typedef struct htree_s htree_t;

struct htree_s {
	/* Binary-Tree */
	bool branch_is_left;
	htree_t *left;
	htree_t *right;
	void *val;
};

uint32_t htree_get_memsize(void);
htree_t* htree_create(nb_membank_t *membank);
htree_t* htree_clone(nb_membank_t *membank, const htree_t *const tree,
		     void* (*clone)(const void *const));
void htree_destroy(nb_membank_t *membank, htree_t* tree);
void htree_destroy_values_recursively(htree_t *tree,
				      void (*destroy)(void*));
htree_t* htree_link(htree_t *tree1, htree_t *tree2,
		    uint32_t (*key)(const void *const));
htree_t* htree_delete_and_get_new_root(nb_membank_t *membank, htree_t *root,
				       uint32_t (*key)(const void *const));
htree_t* htree_containing_val(htree_t *tree, const void *const val,
			      uint32_t (*key)(const void *const),
			      int8_t (*compare)(const void*, const void*));
void htree_delete(nb_membank_t *membank, htree_t *tree,
		  uint32_t (*key)(const void *const));

#endif
