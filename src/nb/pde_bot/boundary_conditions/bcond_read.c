#include "nb/cfreader_cat.h"
#include "nb/container_bot.h"
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_read.h"

#include "bc_atom.h"
#include "bcond_struct.h"

static int read_conditions(nb_container_t *cnt, uint8_t N_dof,
			   nb_cfreader_t *cfreader);

int nb_bcond_read(nb_bcond_t *bcond, nb_cfreader_t *cfreader)
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

static int read_conditions(nb_container_t *cnt, uint8_t N_dof,
			   nb_cfreader_t *cfreader)
{
	int status = 1;
	unsigned int N;
	if (0 != nb_cfreader_read_uint(cfreader, &N))
		goto EXIT;
    printf("N: %d\n", N); /* TEMPORAL */
	for (unsigned int i = 0; i < N; i++) {
		bc_atom_t *bc = bc_atom_create(N_dof);
		nb_container_insert(cnt, bc);
		if (0 != nb_cfreader_read_uint(cfreader, &(bc->id)))
			goto CLEANUP;
		for (uint8_t j = 0; j < N_dof; j++) {
			if (0 != nb_cfreader_read_bool(cfreader,
							&(bc->mask[j])))
				goto CLEANUP;
		}
		printf("N_dof: %d\n", N_dof); /* TEMPORAL */
		for (uint8_t j = 0; j < N_dof; j++) {
			if (bc->mask[j]) {
				double value = 0.0;
				if (0 != nb_cfreader_read_double(cfreader,
								  &value))
					goto CLEANUP;
				bc->val[j] = value;
				printf("%lf ", bc->val[j]); /* TEMPORAL */
			}
		}
	}
	status = 0;
	goto EXIT;
CLEANUP:
	nb_container_clear(cnt);
EXIT:
	return status;
}
