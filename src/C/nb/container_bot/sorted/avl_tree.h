#ifndef __NB_CONTAINER_BOT_SORTED_AVL_TREE_H__
#define __NB_CONTAINER_BOT_SORTED_AVL_TREE_H__

#include <stdbool.h>
#include <stdint.h>

typedef struct tree_s tree_t;

struct tree_s {
	/* Binary-Tree */
	int32_t height;
	tree_t *left;
	tree_t *right;
	void *val;
};

tree_t* tree_create(void);
void tree_destroy(tree_t* tree, void (*destroy)(void*));
void tree_destroy_recursively(tree_t* tree, void (*destroy)(void*));
tree_t* tree_clone(const tree_t *const tree,
		   void* (*clone)(const void *const));
void* tree_exist(const tree_t *const tree, const void *val,
		 int8_t (*compare)(const void*, const void*));
bool tree_is_leaf(const tree_t *const tree);
uint32_t tree_get_height(const tree_t *const tree);
bool tree_insert(tree_t *tree, const void *const val,
		 int8_t (*compare)(const void*, const void*));
tree_t* tree_create_leaf(const void *const val);
tree_t* tree_unlink_most_left(tree_t* tree);
tree_t* tree_unlink_most_right(tree_t* tree);
void tree_replace_root(tree_t* tree);
void* tree_delete(tree_t *tree, const void *val,
		  int8_t (*compare)(const void*, const void*));
#endif
