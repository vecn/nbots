#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/pde_bot.h"

#include "bc_atom.h"
#include "bc_get_container.h"

struct nb_bcond_iter_s {
	const bc_atom_t *bc;
	nb_iterator_t *internal;
};

static void* allocate_bcond_iter(void);

uint16_t nb_bcond_iter_get_memsize(void)
{
	return nb_iterator_get_memsize() + sizeof(nb_bcond_iter_t);
}

void nb_bcond_iter_init(void *iter_ptr)
{
	char *memblock = iter_ptr;
	nb_bcond_iter_t *iter = iter_ptr;
	iter->bc = NULL;
	iter->internal = (void*)(memblock + sizeof(nb_bcond_iter_t));
	nb_iterator_init(iter->internal);
}

void nb_bcond_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	nb_bcond_iter_t *iter = iter_ptr;
	const nb_bcond_iter_t *src_iter = src_iter_ptr;
	iter->bc = src_iter->bc;
	nb_iterator_copy(iter->internal, src_iter->internal);
}

void nb_bcond_iter_finish(void *iter_ptr)
{
	nb_bcond_iter_clear(iter_ptr);
}

void* nb_bcond_iter_create(void)
{
	void *iter = allocate_bcond_iter();
	nb_bcond_iter_init(iter);
	return iter;
}

static void* allocate_bcond_iter(void)
{
	uint16_t size = nb_bcond_iter_get_memsize();
	return nb_allocate_mem(size);
}

void* nb_bcond_iter_clone(const void *iter_ptr)
{
	void *iter = allocate_bcond_iter();
	nb_bcond_iter_copy(iter, iter_ptr);
	return iter;
}

void nb_bcond_iter_destroy(void *iter_ptr)
{
	nb_bcond_iter_clear(iter_ptr);
	nb_free_mem(iter_ptr);
}

void nb_bcond_iter_clear(void *iter_ptr)
{
	nb_bcond_iter_t *iter = iter_ptr;
	iter->bc = NULL;
	nb_iterator_clear(iter->internal);
}

void nb_bcond_iter_set_conditions(nb_bcond_iter_t *iter,
				  const nb_bcond_t *bcond,
				  nb_bcond_id type_id,
				  nb_bcond_where type_elem)
{
	nb_container_t *container = 
		nb_bcond_get_container(bcond, type_id, type_elem);
	nb_iterator_set_container(iter->internal, container);
}

bool nb_bcond_iter_has_more(const nb_bcond_iter_t *iter)
{
	return nb_iterator_has_more(iter->internal);
}

void nb_bcond_iter_go_next(nb_bcond_iter_t *iter)
{
	iter->bc = nb_iterator_get_next(iter->internal);
}

uint32_t nb_bcond_iter_get_id(const nb_bcond_iter_t *iter)
{
	return iter->bc->id;
}

bool nb_bcond_iter_get_mask(const nb_bcond_iter_t *iter,
			    uint8_t dof_id)
{
	return iter->bc->mask[dof_id];
}

bool nb_bcond_iter_val_is_function(const nb_bcond_iter_t *iter)
{
	return NULL != iter->bc->fval;
}

void nb_bcond_iter_get_val(const nb_bcond_iter_t *iter, uint8_t N_dof,
			   double *x, double t, double val[])
{
	if (nb_bcond_iter_val_is_function(iter))
		iter->bc->fval(x, t, val);
	else
		memcpy(val, iter->bc->val, N_dof * sizeof(double));
}
