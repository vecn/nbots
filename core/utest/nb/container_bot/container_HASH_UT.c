#define CONTAINER_ID NB_HASH
#define N_ITEMS 200
#include "container_TEMPLATE_UT.c"

static bool check_realloc_if_max_load(void);
static bool check_get_size(void);
static bool check_get_N_collisions(void);
static bool check_get_collisions(void);
static bool check_exist_if_not_with_collision(void);

static uint32_t small_hash_key(const void *const val_ptr);
static bool check_collisions(nb_container_t *collisions);

inline int vcn_test_get_driver_id(void)
{
	return TEMPLATE_get_driver_id();
}

void vcn_test_load_tests(void *tests_ptr)
{
	TEMPLATE_load_tests(tests_ptr);
	vcn_test_add(tests_ptr, check_realloc_if_max_load,
		     "Check realloc() if max load");
	vcn_test_add(tests_ptr, check_get_size,
		     "Check do('get_size')");
	vcn_test_add(tests_ptr, check_get_N_collisions,
		     "Check do('get_N_collisions')");
	vcn_test_add(tests_ptr, check_get_collisions,
		     "Check do('get_collisions')");
	vcn_test_add(tests_ptr, check_exist_if_not_with_collision,
		     "Check exist() if not exist and the key collides");
}

static bool check_realloc_if_max_load(void)
{
	int N = 3000;
	nb_container_t *cnt = get_container(N);
	uint32_t length = nb_container_get_length(cnt);
	nb_container_destroy(cnt);
	return (N == length);  
}

static bool check_get_size(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	int8_t status;
	uint32_t *size_ptr = nb_container_do(cnt, "get_size", NULL, &status);
	uint32_t size = 0;
	if (0 == status) {
		size = *size_ptr;
		free(size_ptr);
	}
	nb_container_destroy(cnt);
	return (0 < size);
}

static bool check_get_N_collisions(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	nb_container_set_key_generator(cnt, small_hash_key);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < 100; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		nb_container_insert(cnt, val);
	}
	int8_t status;
	uint32_t *N_collisions_ptr =
		nb_container_do(cnt, "get_N_collisions", NULL, &status);
	uint32_t N_collisions = 0;
	if (0 == status) {
		N_collisions = *N_collisions_ptr;
		free(N_collisions_ptr);
	}
	return (3 == N_collisions);
}

static bool check_get_collisions(void)
{
	nb_container_t *cnt = nb_container_create(NB_HASH);
	nb_container_set_key_generator(cnt, small_hash_key);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < 100; i++) {
		int32_t *val = malloc(sizeof(*val));
		*val = i;
		nb_container_insert(cnt, val);
	}
	int8_t status;
	nb_container_t *collisions =
		nb_container_do(cnt, "get_collisions", NULL, &status);
	bool is_ok = false;
	if (0 == status) {
		is_ok = check_collisions(collisions);
		nb_container_destroy(collisions);
	}
	return is_ok;
}


static bool check_exist_if_not_with_collision(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = get_container(N);
	nb_container_set_comparer(cnt, compare);
	int8_t status;
	uint32_t *size = nb_container_do(cnt, "get_size", NULL, &status);
	int32_t to_find = *size + N_ITEMS; /* Not exist and key collides */
	free(size);
	int32_t *val = nb_container_exist(cnt, &to_find);
	bool is_ok = (NULL == val);
	nb_container_destroy(cnt);
	return is_ok;
}

static uint32_t small_hash_key(const void *const val_ptr)
{
	const int32_t *const val = val_ptr;
	return *val % 3;
}

static bool check_collisions(nb_container_t *collisions)
{
	bool is_ok = true;
	while (nb_container_is_not_empty(collisions)) {
		nb_container_t *items = 
			nb_container_delete_first(collisions);
		
		int32_t *val = nb_container_delete_first(items);
		uint32_t key = small_hash_key(val);
		while (nb_container_is_not_empty(items)) {
			val = nb_container_delete_first(items);
			if (small_hash_key(val) != key) {
				is_ok = false;
				break;
			}
		}
		nb_container_destroy(items);
	}
	return is_ok;
}
