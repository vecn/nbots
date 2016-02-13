#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/container_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_iter.h"

#include "bc_atom.h"
#include "bc_get_container.h"

struct nb_bcond_iter_s {
	bc_atom_t *bc;
	vcn_iterator_t *internal;
};

static void* malloc_bcond_iter(void);

inline uint16_t nb_bcond_iter_get_memsize(void)
{
	return vcn_iterator_get_memsize() + sizeof(nb_bcond_iter_t);
}

inline void nb_bcond_iter_init(void *iter_ptr)
{
	char *memblock = iter_ptr;
	nb_bcond_iter_t *iter = iter_ptr;
	iter->bc = NULL;
	iter->internal = memblock + sizeof(nb_bcond_iter_t);
	vcn_iterator_init(iter->internal);
}

inline void nb_bcond_iter_copy(void *iter_ptr, const void *src_iter_ptr)
{
	nb_bcond_iter_t *iter = iter_ptr;
	nb_bcond_iter_t *src_iter = src_iter_ptr;
	iter->bc = src_iter->bc;
	vcn_iterator_copy(iter->internal, src_iter->internal);
}

inline void nb_bcond_iter_clear(void *iter_ptr)
{
	nb_bcond_iter_t *iter = iter_ptr;
	iter->bc = NULL;
	vcn_iterator_clear(iter->internal);
}

inline void* nb_bcond_iter_create(void)
{
	void *iter = malloc_bcond_iter();
	nb_bcond_iter_init(iter);
}

static void* malloc_bcond_iter(void)
{
	uint16_t size = nb_bcond_iter_get_memsize();
	return malloc(size);
}

inline void* nb_bcond_iter_clone(const void *iter_ptr)
{
	void *iter = malloc_bcond_iter();
	nb_bcond_iter_copy(iter, iter_ptr);
	return iter;
}

inline void nb_bcond_iter_destroy(void *iter_ptr)
{
	nb_bcond_iter_clear(iter_ptr);
	free(iter_ptr);
}

inline void nb_bcond_iter_set_conditions(nb_bcond_iter_t *iter,
					 nb_bcond_t *bcond,
					 nb_bcond_id type_id,
					 nb_bcond_where type_elem)
{
	vcn_container_t *container = 
		nb_bcond_get_container(bcond, type_id, type_elem);
	vcn_iterator_set_container(iter->internal, container);
}

inline bool nb_bcond_iter_has_more(const nb_bcond_iter_t *iter)
{
	vcn_iterator_has_more(iter->internal);
}

inline void nb_bcond_iter_go_next(nb_bcond_iter_t *iter)
{
	iter->bc = vcn_iterator_get_next(iter->internal);
}

inline uint32_t nb_bcond_iter_get_id(const nb_bcond_iter_t *iter)
{
	return iter->bc->id;
}

inline bool nb_bcond_iter_get_mask(const nb_bcond_iter_t *iter,
				   uint8_t dof_id)
{
	return iter->bc->mask[dof_id];
}

inline double nb_bcond_iter_get_val(const nb_bcond_iter_t *iter,
				    uint8_t dof_id)
{
	return iter->bc->value[dof_id];
}
