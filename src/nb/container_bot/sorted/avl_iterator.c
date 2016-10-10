#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "avl_tree.h"
#include "avl_struct.h"
#include "avl_iterator.h"

typedef struct {
	/* Inorder iterator (non-recursive) */
	const tree_t *root;
	int32_t parent_id;
	int32_t N_parent;
	tree_t** parent;
	tree_t* tree;
} iter_t;

static void* allocate_iter(void);
static void jump_to_most_left(iter_t *iter);
static void jump_to_parent(iter_t *iter);
static void jump_to_next(iter_t *iter);

inline uint16_t avl_iter_get_memsize(void)
{
	return sizeof(iter_t);
}

inline void avl_iter_init(void *iter_ptr)
{
	memset(iter_ptr, 0, avl_iter_get_memsize());
}

void avl_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	iter_t *iter = iter_ptr;
	const iter_t *src_iter = src_iter_ptr;
	iter->root = src_iter->root;
	iter->tree = src_iter->tree;
	iter->parent_id = src_iter->parent_id;
	iter->N_parent = src_iter->N_parent;
	if (0 < iter->N_parent) {
		iter->parent = nb_allocate_mem(iter->N_parent *
					       sizeof(*(iter->parent)));
		memcpy(iter->parent, src_iter->parent, 
		       iter->N_parent * sizeof(*(iter->parent)));
	} else {
		iter->parent = NULL;
	}
}

inline void avl_iter_finish(void *iter_ptr)
{
	avl_iter_clear(iter_ptr);
}

inline void* avl_iter_create(void)
{
	void *iter = allocate_iter();
	avl_iter_init(iter);
	return iter;
}

static inline void* allocate_iter(void)
{
	uint16_t size = avl_iter_get_memsize();
	return nb_allocate_mem(size);
}

inline void* avl_iter_clone(const void *iter_ptr)
{
	void *iter = allocate_iter();
	avl_iter_copy(iter, iter_ptr);
	return iter;
}

inline void avl_iter_destroy(void *iter_ptr)
{
	avl_iter_finish(iter_ptr);
	nb_free_mem(iter_ptr);
}

inline void avl_iter_clear(void *iter_ptr)
{
	iter_t *iter = iter_ptr;
	if (iter->N_parent > 0)
       		nb_free_mem(iter->parent);
}

void avl_iter_set_dst(void *iter_ptr, const void *avl_ptr)
{
  	iter_t *iter = iter_ptr;
	const avl_t *avl = avl_ptr;
	iter->root = avl->root;
	iter->tree = (tree_t*) iter->root;
	if (NULL != iter->root) {
		iter->N_parent = tree_get_height(iter->root);
		iter->parent = nb_allocate_zero_mem(iter->N_parent *
						    sizeof(tree_t*));
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
