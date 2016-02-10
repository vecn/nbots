/******************************************************************************
 *   AVL iterator                                                             *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "avl_tree.h"
#include "avl_dst.h"
#include "avl_iterator.h"

typedef struct {
	/* Inorder iterator (non-recursive) */
	const tree_t *root;
	int32_t parent_id;
	int32_t N_parent;
	tree_t** parent;
	tree_t* tree;
} iter_t;

static void jump_to_most_left(iter_t *iter);
static void jump_to_parent(iter_t *iter);
static void jump_to_next(iter_t *iter);

inline void* avl_iter_create(void)
{
	return calloc(1, sizeof(iter_t));
}

void avl_iter_set_dst(void *iter_ptr, const void *const avl_ptr)
{
  	iter_t *iter = iter_ptr;
	iter->root = avl_get_iterator_start(avl_ptr);
	iter->tree = (tree_t*) iter->root;
	if (NULL != iter->root) {
		iter->N_parent = tree_get_height(iter->root);
		iter->parent = calloc(iter->N_parent, sizeof(tree_t*));
		jump_to_most_left(iter);
	}
}

static void jump_to_most_left(iter_t *iter)
{
	while (NULL != iter->tree->left) {
		iter->parent[iter->parent_id] = iter->tree;
		iter->parent_id += 1;
		iter->tree = iter->tree->left;
	}
}

void* avl_iter_clone(const void *const iter_ptr)
{
  	const iter_t *const restrict iter = iter_ptr;
	iter_t* clone = calloc(1, sizeof(iter_t));
	clone->root = iter->root;
	clone->tree = iter->tree;
	clone->parent_id = iter->parent_id;
	clone->N_parent = iter->N_parent;
	if (clone->N_parent > 0) {
		clone->parent = calloc(clone->N_parent, sizeof(tree_t*));
		memcpy(clone->parent, iter->parent, 
		       clone->N_parent * sizeof(tree_t*));
	}
	return clone;
}

void avl_iter_destroy(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	if (iter->N_parent > 0)
       		free(iter->parent);
	free(iter);
}

void avl_iter_restart(void *iter_ptr)
{
  	iter_t *iter = iter_ptr;
	iter->parent_id = 0;
	iter->tree = (tree_t*) iter->root;
	if (NULL != iter->root) {
		memset(iter->parent, 0, iter->N_parent * sizeof(tree_t*));
		jump_to_most_left(iter);
	}
}


const void* avl_iter_get_next(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	if (NULL == iter->tree)
		jump_to_parent(iter);

	void* val = NULL;
	if (NULL != iter->tree) {
		val = iter->tree->val;
		jump_to_next(iter);
	}
	return val;
}

static void jump_to_parent(iter_t *iter)
{
	iter->parent_id -= 1;
	if (iter->parent_id >= 0)
		iter->tree = (tree_t*)iter->parent[iter->parent_id];
	else
		iter->tree = NULL;
}

static void jump_to_next(iter_t *iter)
{
	if (iter->tree->right == NULL) {
		iter->tree = NULL;
	} else {
		iter->tree = iter->tree->right;
		jump_to_most_left(iter);
	}
}

inline bool avl_iter_has_more(const void *const iter_ptr)
{
	const iter_t *const restrict iter = iter_ptr;
	return (iter->parent_id > 0 || iter->tree != NULL);
}
