#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "nb/cfreader_cat.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/pde_bot/boundary_conditions.h"

#include "bc_atom.h"

#define CONTAINER_ID NB_CONTAINER_QUEUE

struct nb_bcond_s {
	uint8_t N_dof; /* Degrees of freedom */
	vcn_container_t *dirichlet_vtx;
	vcn_container_t *neumann_vtx;
	vcn_container_t *dirichlet_sgm;
	vcn_container_t *neumann_sgm;
};

static void init_containers(nb_bcond_t *bcond);
static void copy_containers(nb_bcond_t *bcond,
			    const nb_bcond_t *const bcond_src);
static void clone_container_elements(vcn_container_t *cnt,
				     const vcn_container_t *const cnt_src,
				     uint8_t N_dof);
static void* malloc_bcond(void);
static int read_conditions(vcn_container_t *cnt, uint8_t N_dof,
			   vcn_cfreader_t *cfreader);

static void push_dirichlet(nb_bcond_t *bcond,
			   nb_bcond_where type_elem, uint32_t elem_id,
			   const bool dof_mask[], const double value[]);
static void push_neumann(nb_bcond_t *bcond,
			 nb_bcond_where type_elem, uint32_t elem_id,
			 const bool dof_mask[], const double value[]);

uint16_t nb_bcond_get_memsize(uint8_t N_dof)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = vcn_container_get_memsize(CONTAINER_ID);
	return size + 4 * size_cnt;
}

inline void nb_bcond_init(void *bcond_ptr, uint8_t N_dof)
{
	uint16_t size = nb_bcond_get_memsize(N_dof);
	memset(bcond_ptr, 0, size);
	nb_bcond_t *bcond = bcond_ptr;
	bcond->N_dof = N_dof;
	init_containers(bcond);	
}

static void init_containers(nb_bcond_t *bcond)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = vcn_container_get_memsize(CONTAINER_ID);
	char *memblock = bcond;
	bcond->dirichlet_vtx = memblock + size;
	bcond->neumann_vtx = memblock + size + size_cnt;
	bcond->dirichlet_sgm = memblock + size + 2 * size_cnt;
	bcond->neumann_sgm = memblock + size + 3 * size_cnt;

	vcn_container_init(bcond->dirichlet_vtx, CONTAINER_ID);
	vcn_container_set_destroyer(bcond->dirichlet_vtx, bc_atom_destroy);

	vcn_container_init(bcond->neumann_vtx, CONTAINER_ID);
	vcn_container_set_destroyer(bcond->neumann_vtx, bc_atom_destroy);

	vcn_container_init(bcond->dirichlet_sgm, CONTAINER_ID);
	vcn_container_set_destroyer(bcond->dirichlet_sgm, bc_atom_destroy);

	vcn_container_init(bcond->neumann_sgm, CONTAINER_ID);
	vcn_container_set_destroyer(bcond->neumann_sgm, bc_atom_destroy);
}

void nb_bcond_copy(void *bcond_ptr, const void *const src_bcond_ptr)
{
	nb_bcond_t *bcond = bcond_ptr;
	bcond_clone->N_dof = bcond->N_dof;
	copy_containers(bcond_clone, bcond);	
}


static void copy_containers(nb_bcond_t *bcond,
			    const nb_bcond_t *const bcond_src)
{
	clone_container_elements(bcond->dirichlet_vtx,
				 bcond_src->dirichlet_vtx, bcond->N_dof);
	clone_container_elements(bcond->neumann_vtx,
				 bcond_src->neumann_vtx, bcond->N_dof);
	clone_container_elements(bcond->dirichlet_sgm,
				 bcond_src->dirichlet_sgm, bcond->N_dof);
	clone_container_elements(bcond->neumann_sgm,
				 bcond_src->neumann_sgm, bcond->N_dof);
}

static void clone_container_elements(vcn_container_t *cnt,
				     const vcn_container_t *const cnt_src,
				     uint8_t N_dof)
{
	uint16_t size = vcn_iterator_get_memsize();
	vcn_iterator_t *iter = alloca(size);
	vcn_iterator_init(iter);
	vcn_iterator_set_container(iter, cnt_src);
	while (vcn_iterator_has_more(iter)) {
		const bc_atom_t *bc = vcn_iterator_get_next(iter);
		bc_atom_t *bc_clone = bc_atom_clone(bc, N_dof);
		vcn_container_insert(cnt, bc_clone);
	}
	vcn_iterator_clear(iter);
}

void nb_bcond_clear(void *bcond_ptr)
{
	nb_bcond_t *bcond = bcond_ptr;
	vcn_container_clear(bcond->dirichlet_vtx);
	vcn_container_clear(bcond->neumann_vtx);
	vcn_container_clear(bcond->dirichlet_sgm);
	vcn_container_clear(bcond->neumann_sgm);
}

inline void* nb_bcond_create(uint8_t N_dof)
{
	void *bcond = malloc_bcond();
	nb_bcond_init(bcond_ptr, N_dof);
	return bcond;
}

static void* malloc_bcond(void)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = vcn_container_get_memsize(CONTAINER_ID);
	return malloc(size + 4 * size_cnt);
}

void* nb_bcond_clone(const void *const bcond)
{
	void *bcond_clone = malloc_bcond();
	nb_bcond_copy(bcond_clone, bcond);
	return bcond_clone;
}


void nb_bcond_destroy(nb_bcond_t *bcond)
{
	nb_bcond_clear(bcond);
	free(bcond);
}

inline uint8_t nb_bcond_get_N_dof(const nb_bcond_t *const bcond)
{
	return bcond->N_dof;
}

int nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader)
{
	int status = 1;
	if (0 != read_conditions(bcond->dirichlet_vtx, bcond->N_dof, cfreader))
		goto EXIT;
	if (0 != read_conditions(bcond->neumann_vtx, bcond->N_dof, cfreader))
		goto EXIT;
	if (0 != read_conditions(bcond->dirichlet_sgm, bcond->N_dof, cfreader))
		goto EXIT;
	if (0 != read_conditions(bcond->neumann_sgm, bcond->N_dof, cfreader))
		goto EXIT;
	status = 0;
EXIT:
	return status;
}

static int read_conditions(vcn_container_t *cnt, uint8_t N_dof,
			   vcn_cfreader_t *cfreader)
{
	int status = 1;
	unsigned int N;
	if (0 != vcn_cfreader_read_uint(cfreader, &N))
		goto EXIT;

	for (unsigned int i = 0; i < N; i++) {
		bc_atom_t *bc = bc_atom_create(N_dof);
		vcn_container_insert(cnt, bc);
		if (0 != vcn_cfreader_read_uint(cfreader, &(bc->id)))
			goto CLEANUP;
		for (uint8_t j = 0; j < N_dof; j++) {
			if (0 != vcn_cfreader_read_bool(cfreader,
							&(bc->mask[j])))
				goto CLEANUP;
		}
		for (uint8_t j = 0; j < N_dof; j++) {
			if (0 != vcn_cfreader_read_double(cfreader,
							  &(bc->val[j])))
				goto CLEANUP;
		}
	}
	status = 0;
	goto EXIT;
CLEANUP:
	vcn_container_clear(cnt);
EXIT:
	return status;
}

void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   const bool dof_mask[], const double value[])
{
	switch (type_id) {
	case NB_DIRICHLET:
		push_dirichlet(bcond, type_elem, elem_id, dof_mask, value);
		break;
	case NB_NEUMANN:
		push_neumann(bcond, type_elem, elem_id, dof_mask, value);
		break;
	}
}

static void push_dirichlet(nb_bcond_t *bcond,
			   nb_bcond_where type_elem, uint32_t elem_id,
			   const bool dof_mask[], const double value[])
{
	switch (type_elem) {
	case NB_BC_ON_POINT:
		bc_atom_t *bc = bc_atom_create(bcond->N_dof);
		bc_atom_set_data(bc, elem_id, dof_mask, value);
		vcn_container_insert(bcond->dirichlet_vtx, bc);
		break;
	case NB_BC_ON_SEGMENT:
		bc_atom_t *bc = bc_atom_create(bcond->N_dof);
		bc_atom_set_data(bc, elem_id, dof_mask, value);
		vcn_container_insert(bcond->dirichlet_sgm, bc);
		break;
	}
}

static void push_neumann(nb_bcond_t *bcond,
			 nb_bcond_where type_elem, uint32_t elem_id,
			 const bool dof_mask[], const double value[])
{
	switch (type_elem) {
	case NB_BC_ON_POINT:
		bc_atom_t *bc = bc_atom_create(bcond->N_dof);
		bc_atom_set_data(bc, elem_id, dof_mask, value);
		vcn_container_insert(bcond->neumann_vtx, bc);
		break;
	case NB_BC_ON_SEGMENT:
		bc_atom_t *bc = bc_atom_create(bcond->N_dof);
		bc_atom_set_data(bc, elem_id, dof_mask, value);
		vcn_container_insert(bcond->neumann_sgm, bc);
		break;
	}
}
