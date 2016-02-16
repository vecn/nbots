#define CONTAINER_ID NB_SORTED
#define N_ITEMS 200
#include "container_TEMPLATE_UT.c"

static bool check_delete_last(void);

inline int vcn_test_get_driver_id(void)
{
	return TEMPLATE_get_driver_id();
}

void vcn_test_load_tests(void *tests_ptr)
{
	TEMPLATE_load_tests(tests_ptr);
	vcn_test_add(tests_ptr, check_delete_last,
		     "Check do('delete_last')");
}

static bool check_delete_last(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	nb_container_set_key_generator(cnt, keygen);
	nb_container_set_destroyer(cnt, free);
	for (int i = 0; i < N; i++) {
		int *value = malloc(sizeof(*value));
		*value = i;
		nb_container_insert(cnt, value);
	}
	int8_t status;
	int32_t *val = nb_container_do(cnt, "delete_last", NULL, &status);
	bool is_ok = ((N - 1) == *val);
	free(val);
	bool length_ok = ((N - 1) == nb_container_get_length(cnt));
	nb_container_destroy(cnt);
	return is_ok && length_ok;
}
