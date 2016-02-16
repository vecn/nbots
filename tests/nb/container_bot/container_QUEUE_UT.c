#define CONTAINER_ID NB_QUEUE
#define N_ITEMS 200
#include "container_TEMPLATE_UT.c"

static bool check_insert_first(void);

inline int vcn_test_get_driver_id(void)
{
	return TEMPLATE_get_driver_id();
}

void vcn_test_load_tests(void *tests_ptr)
{
	TEMPLATE_load_tests(tests_ptr);
	vcn_test_add(tests_ptr, check_insert_first,
		     "Check do('insert_first')");
}


static bool check_insert_first(void)
{
	int N = N_ITEMS;
	nb_container_t *cnt = nb_container_create(CONTAINER_ID);
	bool is_ok = true;
	for (int i = 0; i < N; i++) {
		int *value = malloc(sizeof(*value));
		*value = i;
		int8_t status;
		nb_container_do(cnt, "insert_first", value, &status);
		is_ok = is_ok && (0 == status);
	}
	int32_t length = nb_container_get_length(cnt);
	int *first_ptr = nb_container_get_first(cnt);
	int first_val = *first_ptr;
	nb_container_set_destroyer(cnt, free);
	nb_container_destroy(cnt);
	return is_ok && (N == length) && ((N - 1) == first_val);
}
