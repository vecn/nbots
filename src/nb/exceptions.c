#include <stdlib.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"

#include "nb/exceptions.h"

struct  vcn_exception_s {
	jmp_buf buf; /* Used in the try/catch statements */
	int id;      /* Used to identify the exception */
	char* info;  /* Info about the exception */
	nb_container_t *memory_handler; /* Used in the 'finally' statement.
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
	return nb_allocate_zero_mem(sizeof(vcn_exception_t));
}

void vcn_exception_destroy(vcn_exception_t *exception)
{
	vcn_exception_clear_alloc(exception);
	nb_free_mem(exception);
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
			     void* ptr, void (*nb_free)(void*))
{
	nb_container_t *list = exception->memory_handler;
	if (NULL == list) {
		list = nb_container_create(NB_QUEUE);
		exception->memory_handler = list;
	}
	void** alloc_and_free = nb_allocate_mem(2 * sizeof(*alloc_and_free));
	alloc_and_free[0] = ptr;
	alloc_and_free[1] = nb_free;
	nb_container_insert(list, alloc_and_free);
}

void vcn_exception_clear_alloc(vcn_exception_t *exception)
{
	nb_container_t *list = exception->memory_handler;
	if (NULL != list) {
		nb_container_set_destroyer(list, destroy_allocation);
		nb_container_destroy(list);
	}
	exception->memory_handler = NULL;
}
  
static void destroy_allocation(void *allocation)
{
	void **alloc_and_free = allocation;
	void *ptr = alloc_and_free[0];
	void (*free_ptr)(void*) = alloc_and_free[1];
	free_ptr(ptr);
	nb_free_mem(alloc_and_free);
}
