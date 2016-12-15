#ifndef __NB_CONTAINER_BOT_SORTED_AVL_TREE_H__
#define __NB_CONTAINER_BOT_SORTED_AVL_TREE_H__

#include <stdbool.h>
#include <stdint.h>

#include "nb/memory_bot.h"

typedef struct tree_s tree_t;

struct tree_s {
	/* Binary-Tree */
	int32_t height;
	tree_t *left;
	tree_t *right;
	void *val;
};

uint32_t tree_get_memsize(void);
tree_t* tree_create(nb_membank_t *membank);
void tree_destroy(nb_membank_t *membank, tree_t* tree);
void tree_destroy_values_recursively(tree_t* tree, void (*destroy)(void*));
tree_t* tree_clone(nb_membank_t *membank,
		   const tree_t *const tree,
		   void* (*clone)(const void *const));
void* tree_exist(const tree_t *const tree, const void *val,
		 int8_t (*compare)(const void*, const void*));
bool tree_is_leaf(const tree_t *const tree);
uint32_t tree_get_height(const tree_t *const tree);
bool tree_insert(nb_membank_t *membank,
		 tree_t *tree, const void *const val,
		 int8_t (*compare)(const void*, const void*));
tree_t* tree_create_leaf(nb_membank_t *membank, const void *const val);
tree_t* tree_unlink_most_left(tree_t* tree);
tree_t* tree_unlink_most_right(tree_t* tree);
void tree_replace_root(nb_membank_t *membank, tree_t* tree);
void* tree_delete(nb_membank_t *membank,
		  tree_t *tree, const void *val,
		  int8_t (*compare)(const void*, const void*));
#endif
