#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "nb/cfreader_cat.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/pde_bot/boundary_conditions.h"

typedef struct {
	uint32_t id;
	double val;
} bc1_t;

typedef struct {
	uint32_t id;
	bool dof_mask[2];
	double val[2];
} bc2_t;

typedef struct {
	uint32_t id;
	bool dof_mask[3];
	double val[3];
} bc3_t;

typedef struct {
	uint32_t id;
	bool *dof_mask;
	double *val;
} bcN_t;

struct nb_bcond_s {
	uint8_t N_dof; /* Degrees of freedom */
	vcn_container_t *dirichlet_on_vtx;
	vcn_container_t *neumann_on_vtx;
	vcn_container_t *dirichlet_on_sgm;
	vcn_container_t *neumann_on_sgm;
};

inline nb_bcond_t* nb_bcond_create(uint8_t N_dof)
{
	nb_cond_t *bc = calloc(1, sizeof(*bc));
	bc->N_dof = N_dof;
	return bc;
}
