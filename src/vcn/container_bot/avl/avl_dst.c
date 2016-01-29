/******************************************************************************
 *   AVL DST: AVL Tree, a self-balanced binary search tree                    *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "avl_tree.h"
#include "avl_dst.h"

typedef struct {
	uint32_t length;
	tree_t* root;
} avl_t;

static bool is_not_empty(const avl_t *const avl);
static void null_destroy(void *val);
static void* delete_root(avl_t *avl);

inline void* avl_create(void)
{
	return calloc(1, sizeof(avl_t));
}

void* avl_clone(const void *const avl_ptr,
		void* (*clone)(const void *const))
{
	const avl_t *const restrict avl = avl_ptr;
	avl_t *avl_clone = avl_create();
	if (is_not_empty(avl)) {
		avl_clone->root = tree_clone(avl->root, clone);
		avl_clone->length = avl->length;
	}
	return avl_clone;
}

void avl_merge(void *avl1_ptr, void *avl2_ptr,
	       uint32_t (*key)(const void *const))
{
	avl_t *avl1 = avl1_ptr;
	avl_t *avl2 = avl2_ptr;
	while (is_not_empty(avl2)) {
		void *val = avl_delete_first(avl2, key);
		avl_insert(avl1, val, key);
	}
}

static inline bool is_not_empty(const avl_t *const restrict avl)
{
	return (NULL != avl->root);
}

inline void avl_destroy(void *avl_ptr,
			void (*destroy)(void*))
{
	avl_clear(avl_ptr, destroy);
	free(avl_ptr);
}

inline void avl_clear(void *avl_ptr,
		      void (*destroy)(void*))
{
	avl_t *avl = avl_ptr;
	if (is_not_empty(avl)) {
		tree_destroy_recursively(avl->root, destroy);
		avl->root = NULL;
		avl->length = 0;
	}
}

bool avl_insert(void *avl_ptr, const void *const restrict val,
		uint32_t (*key)(const void *const))
{
	avl_t *avl = avl_ptr;
	bool success;
	if (is_not_empty(avl)) {
		success = tree_insert(avl->root, val, key);
	} else {
		avl->root = tree_create_leaf(val);
		success = true;
	}
	if (success)
		avl->length += 1;
	return success;
}

void* avl_get_first(const void *const avl_ptr)
{
	const avl_t *const restrict avl = avl_ptr;
	void *val = NULL;
	if (is_not_empty(avl)) {
		tree_t* iter = avl->root;
		while (NULL != iter->left)
			iter = iter->left;
		val = iter->val;
	}
	return val;
}

void* avl_delete_first(void *avl_ptr, uint32_t (*key)(const void *const))
{
	avl_t *avl = avl_ptr;
	void *val = NULL;
	if (is_not_empty(avl)) {
		tree_t *most_left = tree_unlink_most_left(avl->root);
		if (most_left == avl->root)
			avl->root = most_left->right;
		val = most_left->val;
		tree_destroy(most_left, null_destroy);
		avl->length -= 1;
	}
	return val;
}

void* avl_delete_last(void *avl_ptr, uint32_t (*key)(const void *const))
{
	avl_t *avl = avl_ptr;
	void *val = NULL;
	if (is_not_empty(avl)) {
		tree_t *most_right = tree_unlink_most_right(avl->root);
		if (most_right == avl->root)
			avl->root = most_right->left;
		val = most_right->val;
		tree_destroy(most_right, null_destroy);
		avl->length -= 1;
	}
	return val;
}

static inline void null_destroy(void *val)
{
  ; /* Null statement */
}

void* avl_exist(const void *const avl_ptr, const void *const val,
		uint32_t (*key)(const void *const),
		bool (*are_equal)(const void *const, const void *const))
{
	const avl_t *const restrict avl = avl_ptr;
	void *existing_val = NULL;
	if (is_not_empty(avl))
		existing_val = tree_exist(avl->root, val, key, are_equal);
	return existing_val;
}

void* avl_delete(void *avl_ptr, const void *const val,
		 uint32_t (*key)(const void *const),
		 bool (*are_equal)(const void *const, const void *const))
{
	avl_t *avl = avl_ptr;
	void *deleted_val = NULL;
	if(is_not_empty(avl)) {
		if(are_equal(avl->root->val, val))
			deleted_val = delete_root(avl);
		else
			deleted_val = tree_delete(avl->root, val, 
						  key, are_equal);
	}
	if (NULL != deleted_val)
		avl->length -= 1;
	return deleted_val;
}

static void* delete_root(avl_t *avl)
{
	void *val = avl->root->val;
	if (tree_is_leaf(avl->root)) {
		tree_destroy(avl->root, null_destroy);
		avl->root = NULL;
	} else {
		tree_replace_root(avl->root);
	}
	return val;
}

inline uint32_t avl_get_length(const void *const avl_ptr)
{
	const avl_t *const restrict avl = avl_ptr;
	return avl->length;
}

inline bool avl_is_empty(const void *const avl_ptr)
{
	const avl_t *const restrict avl = avl_ptr;
	return (NULL == avl->root);
}

inline bool avl_is_not_empty(const void *const restrict avl_ptr)
{
	return is_not_empty(avl_ptr);
}

inline const void* avl_get_iterator_start(const void *const avl_ptr)
{
	const avl_t *const restrict avl = avl_ptr;
	return avl->root;
}
