#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "nb/memory_bot.h"
#include "nb/geometric_bot.h"

#include "integration_mesh.h"

#define INTEGRATOR_TYPE NB_TRIAN

uint32_t nb_cvfa_get_integration_mesh_memsize(void)
{
	return nb_partition_get_memsize(INTEGRATOR_TYPE);
}

void nb_cvfa_init_integration_mesh(nb_partition_t *intmsh)
{
	nb_partition_init(intmsh , INTEGRATOR_TYPE);
}

void nb_cvfa_load_integration_mesh(const nb_partition_t *part,
				   nb_partition_t *intmsh)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);

	uint32_t mesh_size = nb_mesh_get_memsize();
	uint32_t vtx_size = 2 * N_elems * sizeof(double);
	uint32_t perm_size = N_elems * sizeof(uint32_t);
	uint32_t memsize = mesh_size + vtx_size + perm_size;
	char *memblock = NB_SOFT_MALLOC(memsize);
	nb_mesh_t *mesh = (void*) memblock;
	double *vtx = (void*) (memblock + mesh_size);
	uint32_t *perm = (void*) (memblock + mesh_size + vtx_size);
	
	for (uint32_t i = 0; i < N_elems; i++) {
		vtx[i * 2] = nb_partition_elem_get_x(part, i);
		vtx[i*2+1] = nb_partition_elem_get_y(part, i);
	}	

	nb_mesh_init(mesh);
	nb_mesh_get_smallest_ns_alpha_complex(mesh, N_elems, vtx, 0.666);
	nb_partition_load_from_mesh(intmsh, mesh);
	nb_mesh_finish(mesh);

	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id = nb_partition_get_invtx(intmsh, i);
		perm[id] = i;
	}

	nb_partition_set_nodal_permutation(intmsh, perm);

	NB_SOFT_FREE(memsize, memblock);
}
