#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "avl_tree.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static void get_keys(const void *const val1, const void *const val2,
		     uint32_t (*key)(const void *const), uint32_t keys[2]);
static uint32_t key_ptr(const void *const ptr);
static bool insert_in_left(tree_t *tree, const void *const val,
			   uint32_t (*key)(const void *const));
static bool insert_in_right(tree_t *tree, const void *const val,
			    uint32_t (*key)(const void *const));
static void update_height(tree_t *tree);
static void self_balance(tree_t *tree);
static int get_balance(const tree_t *const tree);
static void double_right_rotation(tree_t *tree);
static void simple_right_rotation(tree_t *tree);
static void update_height_after_rotation(tree_t *tree);
static void double_left_rotation(tree_t *tree);
static void simple_left_rotation(tree_t *tree);
static void replace_with_most_right_of_left(tree_t* tree);
static void null_destroy(void *val);
static void replace_with_left(tree_t *tree);
static tree_t* unlink_most_right(tree_t* tree);
static void replace_with_most_left_of_right(tree_t* tree);
static void replace_with_right(tree_t *tree);
static void* delete_from_left(tree_t *tree, const void *const val,
			      uint32_t (*key)(const void *const),
			      bool (*are_equal)(const void *const,
						const void *const));
static void* delete_from_right(tree_t *tree, const void *const val,
			       uint32_t (*key)(const void *const),
			       bool (*are_equal)(const void *const,
						 const void *const));
static void* delete_left(tree_t *tree);
static void* delete_right(tree_t *tree);

inline tree_t* tree_create(void)
{
	return calloc(1, sizeof(tree_t));
}

inline void tree_destroy(tree_t* tree, void (*destroy)(void*))
{
	destroy(tree->val);
	free(tree);
}


void tree_destroy_recursively(tree_t* tree, void (*destroy)(void*))
{
	if (NULL != tree->left)
		tree_destroy_recursively(tree->left, destroy);
	if (NULL != tree->right)
		tree_destroy_recursively(tree->right, destroy);
	tree_destroy(tree, destroy);
}

tree_t* tree_clone(const tree_t *const restrict tree,
		   void* (*clone)(const void *const))
{
	tree_t *tc = tree_create();
	tc->height = tree->height;
	tc->val = clone(tree->val);
	if (NULL != tree->left)
		tc->left = tree_clone(tree->left, clone);
	if (NULL != tree->right)
		tc->right = tree_clone(tree->right, clone);
	return tc;
}


void* tree_exist(const tree_t *const restrict tree,
		 const void *const restrict val,
		 uint32_t (*key)(const void *const),
		 bool (*are_equal)(const void *const, const void *const))
{
	void *existing_val = NULL;
	if (are_equal(tree->val, val)) {
		existing_val = tree->val;
	} else {
		uint32_t keys[2];
		get_keys(val, tree->val, key, keys);
		if (keys[0] < keys[1]) {
			if (NULL != tree->left)
				existing_val = tree_exist(tree->left, val,
							  key, are_equal);
		} else if (keys[0] > keys[1]) {
			if (NULL != tree->right)
				existing_val = tree_exist(tree->right, val,
							  key, are_equal);
		}
	}
	return existing_val;
}

static void get_keys(const void *const val1, const void *const val2,
		     uint32_t (*key)(const void *const), uint32_t keys[2])
{
	keys[0] = key(val1);
	keys[1] = key(val2);

	if (keys[0] == keys[1]) {
		/* Key by pointer to support repeated elements */
		keys[0] = key_ptr(val1);
		keys[1] = key_ptr(val2);
	}
}

static inline uint32_t key_ptr(const void *const restrict ptr)
{
	return (uint32_t)((uintptr_t)ptr);
}

inline bool tree_is_leaf(const tree_t *const restrict tree)
{
	return (NULL == tree->left) && (NULL == tree->right);
}

inline uint32_t tree_get_height(const tree_t *const restrict tree)
{
	return tree->height;
}

bool tree_insert(tree_t *restrict tree, const void *const restrict val,
		 uint32_t (*key)(const void *const))
{
	uint32_t keys[2];
	get_keys(val, tree->val, key, keys);
  
	bool success = false;
	if (keys[0] != keys[1]) {
		if (keys[0] < keys[1])
			success = insert_in_left(tree, val, key);
		else /* if (keys[0] > keys[1]) */
			success = insert_in_right(tree, val, key);
		update_height(tree);
		self_balance(tree);
	}
	return success;
}

static bool insert_in_left(tree_t *tree, const void *const restrict val,
			   uint32_t (*key)(const void *const))
{
	bool success;
	if (NULL != tree->left) {
		success = tree_insert(tree->left, val, key);
	} else {
		tree->left = tree_create_leaf(val);
		success = true;
	}
	return success;
}

tree_t* tree_create_leaf(const void *const val)
{
	tree_t *leaf = tree_create();
	leaf->val = (void*) val;
	leaf->height = 1;
	return leaf;
}

static bool insert_in_right(tree_t *tree, const void *const restrict val,
			    uint32_t (*key)(const void *const))
{
	bool success;
	if (NULL != tree->right) {
		success = tree_insert(tree->right, val, key);
	} else {
		tree->right = tree_create_leaf(val);
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

void tree_replace_root(tree_t* tree)
{
	if (tree->left != NULL)
		replace_with_most_right_of_left(tree);
	else
		replace_with_most_left_of_right(tree);
}

static void replace_with_most_right_of_left(tree_t* tree)
{
	tree_t *replace = unlink_most_right(tree->left);
	tree->val = replace->val;
	if (tree->left == replace)
		replace_with_left(tree);
	tree_destroy(replace, null_destroy);
}

static inline void null_destroy(void *val)
{
	; /* Null statement */
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

static void replace_with_most_left_of_right(tree_t* tree)
{
	tree_t *replace = tree_unlink_most_left(tree->right);
	tree->val = replace->val;
	if (tree->right == replace)
		replace_with_right(tree);
	tree_destroy(replace, null_destroy);
}

static void replace_with_right(tree_t *tree)
{
	tree->right = tree->right->right;
	update_height(tree);
	self_balance(tree);
}

void* tree_delete(tree_t *tree, const void *const val,
		  uint32_t (*key)(const void *const),
		  bool (*are_equal)(const void *const, const void *const))
{
	uint32_t keys[2];
	get_keys(val, tree->val, key, keys);

	void *deleted_val = NULL;
	if (keys[0] < keys[1])
		deleted_val = delete_from_left(tree, val, key, are_equal);
	else
		deleted_val = delete_from_right(tree, val, key, are_equal);

	if (NULL != deleted_val) {
		update_height(tree);
		self_balance(tree);
	}
	return deleted_val;
}

static void* delete_from_left(tree_t *tree, const void *const val,
			      uint32_t (*key)(const void *const),
			      bool (*are_equal)(const void *const,
						const void *const))
{
	void *deleted_val = NULL;
	if (NULL != tree->left) {
		if (are_equal(tree->left->val, val))
			deleted_val = delete_left(tree);
		else
			deleted_val = tree_delete(tree->left, val,
						  key, are_equal);
	}
	return deleted_val;
}

static void* delete_from_right(tree_t *tree, const void *const val,
			       uint32_t (*key)(const void *const),
			       bool (*are_equal)(const void *const,
						 const void *const))
{
	void *deleted_val = NULL;
	if (NULL != tree->right) {
		if (are_equal(tree->right->val, val))
		        deleted_val = delete_right(tree);
		else
			deleted_val = tree_delete(tree->right, val,
						  key, are_equal);
	}
	return deleted_val;
}

static void* delete_left(tree_t *tree)
{
	void *val = tree->left->val;
	if (tree_is_leaf(tree->left)) {
		free(tree->left);
		tree->left = NULL;
	} else {
		tree_replace_root(tree->left);
	}
	return val;
}

static void* delete_right(tree_t *tree)
{
	void *val = tree->right->val;
	if (tree_is_leaf(tree->right)) {
		free(tree->right);
		tree->right = NULL;
	} else {
		tree_replace_root(tree->right);
	}
	return val;
}
