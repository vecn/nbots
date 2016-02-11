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
	vcn_container_t *dirichlet_on_vtx;
	vcn_container_t *neumann_on_vtx;
	vcn_container_t *dirichlet_on_sgm;
	vcn_container_t *neumann_on_sgm;
};

static void* allocate_bcond(void);
static void init_containers(nb_bcond_t *bcond);
static void copy_containers(nb_bcond_t *bcond,
			    const nb_bcond_t *const bcond_src);
static void clone_container_elements(vcn_container_t *cnt,
				     const vcn_container_t *const cnt_src,
				     uint8_t N_dof);

uint16_t nb_bcond_get_memsize(uint8_t N_dof)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = vcn_container_get_memsize(CONTAINER_ID);
	return size + 4 * size_cnt;
}

void nb_bcond_init(void *bcond_ptr);
void nb_bcond_copy(void *bcond_ptr, const void *const src_bcond_ptr);
void nb_bcond_clear(void *bcond_ptr);

inline nb_bcond_t* nb_bcond_create(uint8_t N_dof)
{
	nb_cond_t *bcond = allocate_bcond();
	bcond->N_dof = N_dof;
	init_containers(bcond);
return bcond;
}

static void* allocate_bcond(void)
{
	uint16_t size = sizeof(nb_bcond_t);
	uint16_t size_cnt = vcn_container_get_memsize(CONTAINER_ID);
	char *memblock = calloc(size + 4 * size_cnt, 1);

	nb_bcond_t *bcond = memblock;
	bcond->dirichlet_on_vtx = memblock + size;
	bcond->neumann_on_vtx = memblock + size + size_cnt;
	bcond->dirichlet_on_sgm = memblock + size + 2 * size_cnt;
	bcond->neumann_on_sgm = memblock + size + 3 * size_cnt;

	return bcond;
}

static void init_containers(nb_bcond_t *bcond)
{
	vcn_container_init(bcond->dirichlet_on_vtx, CONTAINER_ID);
	vcn_container_init(bcond->neumann_on_vtx, CONTAINER_ID);
	vcn_container_init(bcond->dirichlet_on_sgm, CONTAINER_ID);
	vcn_container_init(bcond->neumann_on_sgm, CONTAINER_ID);
}

nb_bcond_t* nb_bcond_clone(const nb_bcond_t *const bcond)
{
	nb_cond_t *bcond_clone = allocate_bcond();
	bcond_clone->N_dof = bcond->N_dof;
	copy_containers(bcond_clone, bcond);
	return bcond;
}

static void copy_containers(nb_bcond_t *bcond,
			    const nb_bcond_t *const bcond_src)
{
	clone_container_elements(bcond->dirichlet_on_vtx,
				 bcond_src->dirichlet_on_vtx, bcond->N_dof);
	clone_container_elements(bcond->neumann_on_vtx,
				 bcond_src->neumann_on_vtx, bcond->N_dof);
	clone_container_elements(bcond->dirichlet_on_sgm,
				 bcond_src->dirichlet_on_sgm, bcond->N_dof);
	clone_container_elements(bcond->neumann_on_sgm,
				 bcond_src->neumann_on_sgm, bcond->N_dof);
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

inline uint8_t nb_bcond_get_N_dof(const nb_bcond_t *const bcond)
{
	return bcond->N_dof;
}

void nb_bcond_read(nb_bcond_t *bcond, vcn_cfreader_t *cfreader);
void nb_bcond_copy(nb_bcond_t *dest, const nb_bcond_t *const src);

void nb_bcond_clear(nb_bcond_t *bcond)
{
	vcn_container_clear(bcond->dirichlet_on_vtx);
	vcn_container_clear(bcond->neumann_on_vtx);
	vcn_container_clear(bcond->dirichlet_on_sgm);
	vcn_container_clear(bcond->neumann_on_sgm);
}

void nb_bcond_destroy(nb_bcond_t *bcond)
{
	vcn_container_clear(bcond->dirichlet_on_vtx);
	vcn_container_clear(bcond->neumann_on_vtx);
	vcn_container_clear(bcond->dirichlet_on_sgm);
	vcn_container_clear(bcond->neumann_on_sgm);
	free(bcond);
}

void nb_bcond_push(nb_bcond_t *bcond, nb_bcond_id type_id,
		   nb_bcond_where type_elem, uint32_t elem_id,
		   const bool dof_mask[], const double value[]);
