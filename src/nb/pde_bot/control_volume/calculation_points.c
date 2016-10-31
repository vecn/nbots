#include <stdint.h>

#include "nb/geometric_bot.h"

#include "calculation_points.h"

static void set_elem_cpoint_bnd(const nb_mesh2D_t *mesh, double *xc,
				uint32_t elem_id);
static void set_elem_cpoint_cen(const nb_mesh2D_t *mesh, double *xc,
				uint32_t elem_id);

void nb_cvfa_set_calculation_points(const nb_mesh2D_t *mesh, double *xc)
{
	uint32_t N = nb_mesh2D_get_N_elems(mesh);
	for (uint32_t i = 0; i < N; i++) {
		if (nb_mesh2D_elem_is_boundary(mesh, i))
			set_elem_cpoint_bnd(mesh, xc, i);
		else
			set_elem_cpoint_cen(mesh, xc, i);
	}
}

static void set_elem_cpoint_bnd(const nb_mesh2D_t *mesh, double *xc,
				uint32_t elem_id)
{
	uint16_t N = nb_mesh2D_elem_get_N_adj(mesh, elem_id);
	uint32_t id[2];
	uint16_t bnd_id = N;
	for (uint16_t i = 0; i < N; i++) {
		if (!nb_mesh2D_elem_has_ngb(mesh, elem_id, i)) {
			id[0] = nb_mesh2D_elem_get_adj(mesh, elem_id, i);
			id[1] = nb_mesh2D_elem_get_adj(mesh, elem_id,
						       (i + 1) % N);
			bnd_id = i;
			break;
		}
	}
	if (N != bnd_id) {
		uint16_t prev = (bnd_id == 0) ? (N - 1) : (bnd_id - 1);
		uint16_t next = (bnd_id + 1) % N;
		if (!nb_mesh2D_elem_has_ngb(mesh, elem_id, prev)) {
			xc[elem_id * 2] = nb_mesh2D_node_get_x(mesh, id[0]);
			xc[elem_id*2+1] = nb_mesh2D_node_get_y(mesh, id[0]);
		} else if (!nb_mesh2D_elem_has_ngb(mesh, elem_id, next)) {
			xc[elem_id * 2] = nb_mesh2D_node_get_x(mesh, id[1]);
			xc[elem_id*2+1] = nb_mesh2D_node_get_y(mesh, id[1]);
		} else {
			xc[elem_id * 2] = 0.5 *
				(nb_mesh2D_node_get_x(mesh, id[0]) +
				 nb_mesh2D_node_get_x(mesh, id[1]));
			xc[elem_id*2+1] = 0.5 *
				(nb_mesh2D_node_get_y(mesh, id[0]) +
				 nb_mesh2D_node_get_y(mesh, id[1]));
		}
	} else {
		set_elem_cpoint_cen(mesh, xc, elem_id);
	}
}

static void set_elem_cpoint_cen(const nb_mesh2D_t *mesh, double *xc,
				uint32_t elem_id)
{
	xc[elem_id * 2] = nb_mesh2D_elem_get_x(mesh, elem_id);
	xc[elem_id*2+1] = nb_mesh2D_elem_get_y(mesh, elem_id);
}
