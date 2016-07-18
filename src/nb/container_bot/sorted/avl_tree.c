#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "avl_tree.h"

#include "nb/memory_bot.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static int8_t compare_by_ptr_if_equal(const void *ptr1, const void *ptr2,
				      int8_t (*compare)(const void*,
							const void *));
static int8_t compare_by_ptr(const void *ptr1, const void *ptr2);
static bool insert_in_left(nb_membank_t *membank,
			   tree_t *tree, const void *const val,
			   int8_t (*compare)(const void*, const void*));
static bool insert_in_right(nb_membank_t *membank,
			    tree_t *tree, const void *const val,
			    int8_t (*compare)(const void*, const void*));
static void update_height(tree_t *tree);
static void self_balance(tree_t *tree);
static int get_balance(const tree_t *const tree);
static void double_right_rotation(tree_t *tree);
static void simple_right_rotation(tree_t *tree);
static void update_height_after_rotation(tree_t *tree);
static void double_left_rotation(tree_t *tree);
static void simple_left_rotation(tree_t *tree);
static void replace_with_most_right_of_left(nb_membank_t *membank,
					    tree_t* tree);
static void replace_with_left(tree_t *tree);
static tree_t* unlink_most_right(tree_t* tree);
static void replace_with_most_left_of_right(nb_membank_t *membank,
					    tree_t* tree);
static void replace_with_right(tree_t *tree);
static void* delete_from_left(nb_membank_t *membank,
			      tree_t *tree, const void *val,
			      int8_t (*compare)(const void*, const void*));
static void* delete_from_right(nb_membank_t *membank,
			       tree_t *tree, const void *val,
			       int8_t (*compare)(const void*, const void*));
static void* delete_left(nb_membank_t *membank, tree_t *tree);
static void* delete_right(nb_membank_t *membank, tree_t *tree);

uint32_t tree_get_memsize(void)
{
	return sizeof(tree_t);
}

tree_t* tree_create(nb_membank_t *membank)
{
	return nb_membank_data_calloc(membank);
}

void tree_destroy(nb_membank_t *membank, tree_t* tree)
{
	nb_membank_data_free(membank, tree);
}

void tree_destroy_values_recursively(tree_t* tree, void (*destroy)(void*))
{
	if (NULL != tree->left)
		tree_destroy_values_recursively(tree->left, destroy);
	if (NULL != tree->right)
		tree_destroy_values_recursively(tree->right, destroy);

	destroy(tree->val);
}

tree_t* tree_clone(nb_membank_t *membank,
		   const tree_t *const restrict tree,
		   void* (*clone)(const void *const))
{
	tree_t *tc = tree_create(membank);
	tc->height = tree->height;
	tc->val = clone(tree->val);
	if (NULL != tree->left)
		tc->left = tree_clone(membank, tree->left, clone);
	if (NULL != tree->right)
		tc->right = tree_clone(membank, tree->right, clone);
	return tc;
}


void* tree_exist(const tree_t *const restrict tree,  const void *val,
		 int8_t (*compare)(const void*, const void*))
{
	void *existing_val = NULL;
	int8_t comparison = compare(val, tree->val);
	if (0 == comparison) {
		existing_val = tree->val;
	} else {
		if (comparison < 0) {
			if (NULL != tree->left)
				existing_val = tree_exist(tree->left, val,
							  compare);
		} else {
			if (NULL != tree->right)
				existing_val = tree_exist(tree->right, val,
							  compare);
		}
	}
	return existing_val;
}

static int8_t compare_by_ptr_if_equal(const void *ptr1, const void *ptr2,
				      int8_t (*compare)(const void*,
							const void *))
{
	int8_t comparison = compare(ptr1, ptr2);
	if (0 == comparison)
		comparison = compare_by_ptr(ptr1, ptr2);
	return comparison;
}

static int8_t compare_by_ptr(const void *ptr1, const void *ptr2)
{
	int8_t out;
	if (ptr1 < ptr2)
		out = -1;
	else if (ptr1 > ptr2)
		out = 1;
	else
		out = 0;
	return out;
}

inline bool tree_is_leaf(const tree_t *const restrict tree)
{
	return (NULL == tree->left) && (NULL == tree->right);
}

inline uint32_t tree_get_height(const tree_t *const restrict tree)
{
	return tree->height;
}

bool tree_insert(nb_membank_t *membank,
		 tree_t *restrict tree, const void *const restrict val,
		 int8_t (*compare)(const void*, const void*))
{
	int8_t comparison = compare_by_ptr_if_equal(val, tree->val, compare);
  	bool success = false;
	if (0 != comparison) {
		if (0 > comparison)
			success = insert_in_left(membank, tree, val, compare);
		else /* if (0 < comparison) */
			success = insert_in_right(membank, tree, val, compare);
		update_height(tree);
		self_balance(tree);
	}
	return success;
}

static bool insert_in_left(nb_membank_t *membank,
			   tree_t *tree, const void *const restrict val,
			   int8_t (*compare)(const void*, const void*))
{
	bool success;
	if (NULL != tree->left) {
		success = tree_insert(membank, tree->left, val, compare);
	} else {
		tree->left = tree_create_leaf(membank, val);
		success = true;
	}
	return success;
}

tree_t* tree_create_leaf(nb_membank_t *membank, const void *const val)
{
	tree_t *leaf = tree_create(membank);
	leaf->val = (void*) val;
	leaf->height = 1;
	return leaf;
}

static bool insert_in_right(nb_membank_t *membank,
			    tree_t *tree, const void *const restrict val,
			    int8_t (*compare)(const void*, const void*))
{
	bool success;
	if (NULL != tree->right) {
		success = tree_insert(membank, tree->right, val, compare);
	} else {
		tree->right = tree_create_leaf(membank, val);
		success = true;
	}
	return success;
}

static void update_height(tree_t *restrict tree)
{
	int lheight = 0;
	if (NULL != tree->left)
		lheight = tree->left->height;
	int rheight = 0;
	if (NULL != tree->right)
		rheight = tree->right->height;

	tree->height = 1 + MAX(lheight, rheight);
}

static void self_balance(tree_t *restrict tree)
{
	if (get_balance(tree) < -1) {
		if (get_balance(tree->left) > 0)
			double_right_rotation(tree);
		else
			simple_right_rotation(tree);
		update_height_after_rotation(tree);
	} else if (get_balance(tree) > 1) {
		if(get_balance(tree->right) < 0)
			double_left_rotation(tree);
		else
			simple_left_rotation(tree);
		update_height_after_rotation(tree);
	}
}

static int get_balance(const tree_t *const restrict tree)
{
	int lheight = 0;
	if (tree->left != NULL)
		lheight = tree->left->height;
	int rheight = 0;
	if (tree->right != NULL)
		rheight = tree->right->height;
	return rheight - lheight;
}

static void double_right_rotation(tree_t *restrict tree)
{
	void* aux = tree->right;
	tree->right = tree->left->right;
	tree->left->right = tree->right->left;
	tree->right->left = tree->right->right;
	tree->right->right = aux;
	aux = tree->val;
	tree->val = tree->right->val;
	tree->right->val = aux;
}

static void simple_right_rotation(tree_t *restrict tree)
{
	void* aux = tree->right;
	tree->right = tree->left;
	tree->left = tree->right->left;
	tree->right->left = tree->right->right;
	tree->right->right = aux;
	aux = tree->val;
	tree->val = tree->right->val;
	tree->right->val = aux;
}

static inline void update_height_after_rotation(tree_t *restrict tree)
{
	update_height(tree->left);
	update_height(tree->right);
	update_height(tree);
}

static void double_left_rotation(tree_t *restrict tree)
{
	void* aux = tree->left;
	tree->left = tree->right->left;
	tree->right->left = tree->left->right;
	tree->left->right = tree->left->left;
	tree->left->left = aux;
	aux = tree->val;
	tree->val = tree->left->val;
	tree->left->val = aux;
}

static void simple_left_rotation(tree_t *restrict tree)
{
	void* aux = tree->left;
	tree->left = tree->right;
	tree->right = tree->left->right;
	tree->left->right = tree->left->left;
	tree->left->left = aux;
	aux = tree->val;
	tree->val = tree->left->val;
	tree->left->val = aux;
}

tree_t* tree_unlink_most_left(tree_t* tree)
{  
	tree_t *most_left = tree;
	if (NULL != tree->left) {
		most_left = tree_unlink_most_left(tree->left);

		if (most_left == tree->left) {
			tree->left = most_left->right;
			most_left->right = NULL;
		}

		update_height(tree);
		self_balance(tree);
	}
	return most_left;
}

tree_t* tree_unlink_most_right(tree_t* tree)
{  
	tree_t *most_right = tree;
	if (NULL != tree->right) {
		most_right = tree_unlink_most_right(tree->right);

		if (most_right == tree->right) {
			tree->right = most_right->left;
			most_right->left = NULL;
		}

		update_height(tree);
		self_balance(tree);
	}
	return most_right;
}

void tree_replace_root(nb_membank_t *membank, tree_t* tree)
{
	if (tree->left != NULL)
		replace_with_most_right_of_left(membank, tree);
	else
		replace_with_most_left_of_right(membank, tree);
}

static void replace_with_most_right_of_left(nb_membank_t *membank,
					    tree_t* tree)
{
	tree_t *replace = unlink_most_right(tree->left);
	tree->val = replace->val;
	if (tree->left == replace)
		replace_with_left(tree);
	tree_destroy(membank, replace);
}

static void replace_with_left(tree_t *tree)
{
	tree->left = tree->left->left;
	update_height(tree);
	self_balance(tree);
}

static tree_t* unlink_most_right(tree_t* tree)
{  
	tree_t *most_right = tree;
	if (NULL != tree->right) {
		most_right = unlink_most_right(tree->right);

		if (most_right == tree->right) {
			tree->right = most_right->left;
			most_right->left = NULL;
		}

		update_height(tree);
		self_balance(tree);
	}
	return most_right;
}

static void replace_with_most_left_of_right(nb_membank_t *membank,
					    tree_t* tree)
{
	tree_t *replace = tree_unlink_most_left(tree->right);
	tree->val = replace->val;
	if (tree->right == replace)
		replace_with_right(tree);
	tree_destroy(membank, replace);
}

static void replace_with_right(tree_t *tree)
{
	tree->right = tree->right->right;
	update_height(tree);
	self_balance(tree);
}

void* tree_delete(nb_membank_t *membank,
		  tree_t *tree, const void *val,
		  int8_t (*compare)(const void*, const void*))
{
	int8_t comparison = compare_by_ptr_if_equal(val, tree->val, compare);
	void *deleted_val = NULL;
	if (0 > comparison)
		deleted_val = delete_from_left(membank, tree, val, compare);
	else
		deleted_val = delete_from_right(membank, tree, val, compare);

	if (NULL != deleted_val) {
		update_height(tree);
		self_balance(tree);
	}
	return deleted_val;
}

static void* delete_from_left(nb_membank_t *membank,
			      tree_t *tree, const void *const val,
			      int8_t (*compare)(const void*, const void*))
{
	void *deleted_val = NULL;
	if (NULL != tree->left) {
		int8_t comparison = compare(tree->left->val, val);
		if (0 == comparison)
			deleted_val = delete_left(membank, tree);
		else
			deleted_val = tree_delete(membank,
						  tree->left, val,
						  compare);
	}
	return deleted_val;
}

static void* delete_from_right(nb_membank_t *membank,
			       tree_t *tree, const void *const val,
			       int8_t (*compare)(const void*, const void*))
{
	void *deleted_val = NULL;
	if (NULL != tree->right) {
		int8_t comparison = compare(tree->right->val, val);
		if (0 == comparison)
		        deleted_val = delete_right(membank, tree);
		else
			deleted_val = tree_delete(membank,
						  tree->right, val,
						  compare);
	}
	return deleted_val;
}

static void* delete_left(nb_membank_t *membank, tree_t *tree)
{
	void *val = tree->left->val;
	if (tree_is_leaf(tree->left)) {
		tree_destroy(membank, tree->left);
		tree->left = NULL;
	} else {
		tree_replace_root(membank, tree->left);
	}
	return val;
}

static void* delete_right(nb_membank_t *membank, tree_t *tree)
{
	void *val = tree->right->val;
	if (tree_is_leaf(tree->right)) {
		tree_destroy(membank, tree->right);
		tree->right = NULL;
	} else {
		tree_replace_root(membank, tree->right);
	}
	return val;
}
