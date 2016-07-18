#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "heap_tree.h"

static void set_left_tree(htree_t *tree, htree_t *branch);
static void set_right_tree(htree_t *tree, htree_t *branch);
static htree_t* get_left_tree(const htree_t *const tree);
static htree_t* get_right_tree(const htree_t *const tree);
static void add_link(htree_t *main_tree, htree_t *linked_tree);
static htree_t* pairing_left_to_right(htree_t *tree,
				      uint32_t (*key)(const void *const));
static void make_root(htree_t *tree);
static bool is_not_root(const htree_t *const tree);
static htree_t* get_most_right(htree_t *tree);
static bool is_unique_branch(const htree_t *const tree);
static bool is_left_branch(const htree_t *const tree);
static bool is_right_branch(const htree_t *const tree);
static htree_t* link_right_to_left(htree_t *tree,
				   uint32_t (*key)(const void *const));
static htree_t* get_parent(const htree_t *const restrict tree);

inline htree_t* htree_create(void)
{
	return calloc(1, sizeof(htree_t));
}


htree_t* htree_clone(const htree_t *const tree,
		     void* (*clone)(const void *const))
{
	htree_t *tree_cloned = htree_create();
	tree_cloned->val = clone(tree->val);
	tree_cloned->branch_is_left = tree->branch_is_left;

	htree_t *left_tree = get_left_tree(tree);
	if (NULL != left_tree) {
		htree_t *cloned_left = htree_clone(left_tree, clone);
		set_left_tree(tree_cloned, cloned_left);
	}
	
	htree_t *right_tree = get_right_tree(tree);
	if (NULL != right_tree) {
		htree_t *cloned_right = htree_clone(right_tree, clone);
		set_right_tree(tree_cloned, cloned_right);
	}
	return tree_cloned;
}


static void set_left_tree(htree_t *tree, htree_t *branch)
{
	if (NULL != branch) {
		if (NULL == tree->left) {
			branch->right = tree;
		} else {
			if (tree->branch_is_left)
				branch->right = tree;
			else
				branch->right = tree->left;
		}
	}
	tree->left = branch;
	tree->branch_is_left = true;
}

static void set_right_tree(htree_t *tree, htree_t *branch)
{
	if (NULL != branch)
		branch->right = tree;
	if (NULL == tree->left) {
		tree->branch_is_left = false;
		tree->left = branch;
	} else {
		if (tree->branch_is_left)
			tree->left->right = branch;
		else
			tree->left = branch;
	}
}

void htree_destroy(htree_t* tree)
{
	free(tree);
}

void htree_destroy_recursively(htree_t *tree,
			       void (*destroy)(void*))
{
	if (NULL != tree->left) {
		htree_t *right = tree->left->right;
		htree_destroy_recursively(tree->left, destroy);
		if (right != tree && right != NULL)
			htree_destroy_recursively(right, destroy);
	}
	if (NULL != destroy)
		destroy(tree->val);
	htree_destroy(tree);
}

static htree_t* get_left_tree(const htree_t *const restrict tree)
{
	htree_t *left = NULL;
	if (NULL != tree->left) {
		htree_t *branch = tree->left;
		if (branch->right == tree) {
			if (tree->branch_is_left)
				left = branch;
		} else {
			left = branch;
		}
	}
	return left;
}

static htree_t* get_right_tree(const htree_t *const restrict tree)
{
	htree_t *right = NULL;
	if (NULL != tree->left) {
		htree_t *branch = tree->left;
		if (branch->right == tree) {
			bool branch_is_right = !tree->branch_is_left;
			if (branch_is_right)
				right = branch;
		} else {
			right = branch->right;
		}
	}
	return right;
}

htree_t* htree_link(htree_t *tree1, htree_t *tree2,
		    uint32_t (*key)(const void *const))
{
	uint32_t key1 = key(tree1->val);
	uint32_t key2 = key(tree2->val);
	htree_t *root;
	if (key1 < key2) {
		add_link(tree1, tree2);
		root = tree1;
	} else {
		add_link(tree2, tree1);
		root = tree2;
	}
	return root;
}

static inline void add_link(htree_t *main_tree, htree_t *linked_tree)
{
	htree_t* first_son = get_left_tree(main_tree);
	set_right_tree(linked_tree, first_son);
	set_left_tree(main_tree, linked_tree);
}

htree_t* htree_delete_and_get_new_root(htree_t *root,
				       uint32_t (*key)(const void *const))
{
	htree_t *first_son = get_left_tree(root);
	htree_destroy(root);

	htree_t *new_root = NULL;
	if (NULL != first_son) {
		first_son = pairing_left_to_right(first_son, key);
		make_root(first_son);
		htree_t *last_son = get_most_right(first_son);
		new_root = link_right_to_left(last_son, key);
	}
	return new_root;
}

static htree_t* pairing_left_to_right(htree_t *tree,
				      uint32_t (*key)(const void *const))
{
	htree_t *next_bro = get_right_tree(tree);
	if (NULL != next_bro) {    
		htree_t *next_pair = get_right_tree(next_bro);
		tree = htree_link(tree, next_bro, key);
		if (NULL != next_pair)
			next_pair = pairing_left_to_right(next_pair, key);
		set_right_tree(tree, next_pair);
	}
	return tree;
}

static void make_root(htree_t *restrict tree)
{
	if (is_not_root(tree)) {
		if (is_unique_branch(tree))
			tree->right->left = NULL;
		else if (is_left_branch(tree))
			tree->right->right->left = tree->right;
		else if (is_right_branch(tree))
			tree->right->left->right = tree->right;
		tree->right = NULL;    
	}
}

static inline bool is_not_root(const htree_t *const restrict tree)
{
	return (NULL != tree->right);
}

static inline htree_t* get_most_right(htree_t *tree)
{
	while (NULL != get_right_tree(tree))
		tree = get_right_tree(tree);
	return tree;
}

static inline bool is_unique_branch(const htree_t *const restrict tree)
{
	return (tree->right->left == tree);
}

static inline bool is_left_branch(const htree_t *const restrict tree)
{
	bool is_left = false;
	if (NULL != tree->right->right) {
		if (tree->right->right->left == tree)
			is_left = true;
	}
	return is_left;
}

static inline bool is_right_branch(const htree_t *const restrict tree)
{
	bool is_right = false;
	if (NULL != tree->right->left) {
		if (tree->right->left->right == tree)
			is_right = true;
	}
	return is_right;
}

static htree_t* link_right_to_left(htree_t *tree,
				  uint32_t (*key)(const void *const))
{
	htree_t *parent = get_parent(tree);
	if (NULL != parent) {
		htree_t *grand_parent = get_parent(parent);
		tree = htree_link(tree, parent, key);
		if (NULL != grand_parent) {
			set_right_tree(grand_parent, tree);
			tree = link_right_to_left(tree, key);
		}
	}
	return tree;
}

static htree_t* get_parent(const htree_t *const restrict tree)
{
	htree_t *parent = NULL;
	if (is_not_root(tree)) {
		if (is_left_branch(tree))
			parent = tree->right->right;
		else /* if (tree is the unique or the right branch) */
			parent = tree->right;
	}
	return parent;
}

htree_t* htree_containing_val(htree_t *tree, const void *const val,
			      uint32_t (*key)(const void *const),
			      int8_t (*compare)(const void*, const void*))
{
	htree_t *tree_with_val = NULL;
	if (0 == compare(tree->val, val)) {
		tree_with_val = tree;
	} else {
		uint32_t key_val = key(val);
		uint32_t key_tree = key(tree->val);
		if (key_val >= key_tree /* Half ordered tree */) {
			htree_t *left_tree = get_left_tree(tree);
			if (NULL != left_tree) {
				tree_with_val = htree_containing_val(left_tree,
								     val, key,
								     compare);
			}
		}
		if (NULL == tree_with_val) {
			htree_t *right_tree = get_right_tree(tree);
			if (NULL != right_tree) {
				tree_with_val = htree_containing_val(right_tree,
								     val, key,
								     compare);
			}
		}
	}
	return tree_with_val;
}

void htree_delete(htree_t *tree,
		  uint32_t (*key)(const void *const))
{
	htree_t *parent = get_parent(tree);
	bool is_unique = is_unique_branch(tree);
	bool is_left = is_left_branch(tree);
	tree = htree_delete_and_get_new_root(tree, key);
	if (is_unique || is_left)
		set_left_tree(parent, tree);
	else
		set_right_tree(parent, tree);
}
