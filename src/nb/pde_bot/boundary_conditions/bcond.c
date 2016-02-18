#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <alloca.h>

#include "nb/container_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"

#include "bc_atom.h"
#include "bcond_struct.h"
#include "bc_get_container.h"

#define CONTAINER_ID NB_QUEUE

static void init_containers(nb_bcond_t *bcond);
static void copy_containers(nb_bcond_t *bcond,
			    const nb_bcond_t *const bcond_src);
static void clone_container_elements(nb_container_t *cnt,
				     const nb_container_t *const cnt_src,
				     uint8_t N_dof);
static void* malloc_bcond(void);


uint16_t nb_bcond_get_memsize(uint8_t N_dof)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = nb_container_get_memsize(CONTAINER_ID);
	return size + 4 * size_cnt;
}

inline void nb_bcond_init(void *bcond_ptr, uint8_t N_dof)
{
	nb_bcond_t *bcond = bcond_ptr;
	bcond->N_dof = N_dof;
	init_containers(bcond);	
}

static void init_containers(nb_bcond_t *bcond)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = nb_container_get_memsize(CONTAINER_ID);
	char *memblock = (void*)bcond;
	bcond->dirichlet_vtx = (void*)(memblock + size);
	bcond->neumann_vtx = (void*)(memblock + size + size_cnt);
	bcond->dirichlet_sgm = (void*)(memblock + size + 2 * size_cnt);
	bcond->neumann_sgm = (void*)(memblock + size + 3 * size_cnt);

	nb_container_init(bcond->dirichlet_vtx, CONTAINER_ID);
	nb_container_set_destroyer(bcond->dirichlet_vtx, bc_atom_destroy);

	nb_container_init(bcond->neumann_vtx, CONTAINER_ID);
	nb_container_set_destroyer(bcond->neumann_vtx, bc_atom_destroy);

	nb_container_init(bcond->dirichlet_sgm, CONTAINER_ID);
	nb_container_set_destroyer(bcond->dirichlet_sgm, bc_atom_destroy);

	nb_container_init(bcond->neumann_sgm, CONTAINER_ID);
	nb_container_set_destroyer(bcond->neumann_sgm, bc_atom_destroy);
}

void nb_bcond_copy(void *bcond_ptr, const void *const src_bcond_ptr)
{
	nb_bcond_t *bcond = bcond_ptr;
	const nb_bcond_t *src_bcond = src_bcond_ptr;
	bcond->N_dof = src_bcond->N_dof;
	init_containers(bcond);	
	copy_containers(bcond, src_bcond);	
}


static void copy_containers(nb_bcond_t *bcond, const nb_bcond_t *src_bcond)
{
	clone_container_elements(bcond->dirichlet_vtx,
				 src_bcond->dirichlet_vtx, bcond->N_dof);
	clone_container_elements(bcond->neumann_vtx,
				 src_bcond->neumann_vtx, bcond->N_dof);
	clone_container_elements(bcond->dirichlet_sgm,
				 src_bcond->dirichlet_sgm, bcond->N_dof);
	clone_container_elements(bcond->neumann_sgm,
				 src_bcond->neumann_sgm, bcond->N_dof);
}

static void clone_container_elements(nb_container_t *cnt,
				     const nb_container_t *const cnt_src,
				     uint8_t N_dof)
{
	uint16_t size = nb_iterator_get_memsize();
	nb_iterator_t *iter = alloca(size);
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, cnt_src);
	while (nb_iterator_has_more(iter)) {
		const bc_atom_t *bc = nb_iterator_get_next(iter);
		bc_atom_t *bc_clone = bc_atom_clone(bc, N_dof);
		nb_container_insert(cnt, bc_clone);
	}
	nb_iterator_finish(iter);
}

inline void nb_bcond_finish(void *bcond_ptr)
{
	nb_bcond_clear(bcond_ptr);
}

inline void* nb_bcond_create(uint8_t N_dof)
{
	void *bcond = malloc_bcond();
	nb_bcond_init(bcond, N_dof);
	return bcond;
}

static void* malloc_bcond(void)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = nb_container_get_memsize(CONTAINER_ID);
	return malloc(size + 4 * size_cnt);
}

void* nb_bcond_clone(const void *const bcond)
{
	void *bcond_clone = malloc_bcond();
	nb_bcond_copy(bcond_clone, bcond);
	return bcond_clone;
}


void nb_bcond_destroy(void *bcond)
{
	nb_bcond_finish(bcond);
	free(bcond);
}

void nb_bcond_clear(void *bcond_ptr)
{
	nb_bcond_t *bcond = bcond_ptr;
	nb_container_clear(bcond->dirichlet_vtx);
	nb_container_clear(bcond->neumann_vtx);
	nb_container_clear(bcond->dirichlet_sgm);
	nb_container_clear(bcond->neumann_sgm);
}

inline uint8_t nb_bcond_get_N_dof(const nb_bcond_t *const bcond)
{
	return bcond->N_dof;
}

void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   const bool dof_mask[], const double value[])
{
	nb_container_t *container = 
		nb_bcond_get_container(bcond, type_id, type_elem);
	if (NULL != container) {
		bc_atom_t *bc = bc_atom_create(bcond->N_dof);
		bc_atom_set_data(bc, elem_id, dof_mask, value, bcond->N_dof);
		nb_container_insert(container, bc);
	}
}
