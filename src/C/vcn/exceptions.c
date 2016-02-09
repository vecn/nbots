#include <stdlib.h>

#include "vcn/exceptions.h"
#include "vcn/container_bot/container.h"

struct  vcn_exception_s {
	jmp_buf buf; /* Used in the try/catch statements */
	int id;      /* Used to identify the exception */
	char* info;  /* Info about the exception */
	vcn_container_t *memory_handler; /* Used in the 'finally' statement.
					  * Memory allocated in the try block */
};

static void destroy_allocation(void *allocation);

inline void* vcn_exception_try_USE_MACRO_INSTEAD_OF_THIS
	(vcn_exception_t *exception)
{
  return exception->buf;
}

vcn_exception_t* vcn_exception_create()
{
	return calloc(1, sizeof(vcn_exception_t));
}

void vcn_exception_destroy(vcn_exception_t *exception)
{
	vcn_exception_clear_alloc(exception);
	free(exception);
}

void vcn_exception_throw(vcn_exception_t *exception, int id, void *info)
{
	exception->id = id;
	exception->info = info;
	longjmp(exception->buf, 1);    
}

int vcn_exception_get_id(vcn_exception_t *exception)
{
	return exception->id;
}

char* vcn_exception_get_info(vcn_exception_t *exception)
{
  return exception->info;
}

void vcn_exception_set_alloc(vcn_exception_t *exception,
			     void* ptr, void (*free)(void*))
{
	vcn_container_t *list = exception->memory_handler;
	if (NULL == list) {
		list = vcn_container_create(VCN_CONTAINER_QUEUE);
		exception->memory_handler = list;
	}
	void** alloc_and_free = malloc(2 * sizeof(*alloc_and_free));
	alloc_and_free[0] = ptr;
	alloc_and_free[1] = free;
	vcn_container_insert(list, alloc_and_free);
}

void vcn_exception_clear_alloc(vcn_exception_t *exception)
{
	vcn_container_t *list = exception->memory_handler;
	if (NULL != list) {
		vcn_container_set_destroyer(list, destroy_allocation);
		vcn_container_destroy(list);
	}
	exception->memory_handler = NULL;
}
  
static void destroy_allocation(void *allocation)
{
	void **alloc_and_free = allocation;
	void *ptr = alloc_and_free[0];
	void (*free_ptr)(void*) = alloc_and_free[1];
	free_ptr(ptr);
	free(alloc_and_free);
}
