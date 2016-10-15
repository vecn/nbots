#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/mesh2D.h"
#include "nb/geometric_bot/mesh/dewall.h"
#include "nb/geometric_bot/mesh/alpha_shape.h"

#include "mesh2D_structs.h"

#define POW2(a) ((a)*(a))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void get_alpha_complex(nb_tessellator2D_t *mesh, double alpha);
static double get_minmax_radius(const nb_tessellator2D_t *mesh);
static void get_max_length_x_trg(const nb_tessellator2D_t *mesh, double *max_t);
static double trg_get_max_length(msh_trg_t *trg);
static void get_min_trg_x_vtx(const nb_tessellator2D_t *mesh, const double *max_t,
			      msh_trg_t **min_v);
static void verify_min_on_vtx(const double *max_t, const msh_trg_t *trg,
			      msh_trg_t **min_v, uint32_t vid);
static msh_trg_t *get_max_trg(const double *max_t, uint32_t N_vtx,
			      msh_trg_t **min_v);

void nb_tessellator2D_get_alpha_complex(nb_tessellator2D_t *mesh, uint32_t N_vertices,
			       const double *const vertices,
			       double alpha)
{
	nb_tessellator2D_get_delaunay(mesh, N_vertices, vertices);
	get_alpha_complex(mesh, alpha);
}

static void get_alpha_complex(nb_tessellator2D_t *mesh, double alpha)
{
	uint32_t cnt_type = NB_QUEUE;

	uint32_t iter_size = nb_iterator_get_memsize();
	uint32_t cnt_size = nb_container_get_memsize(cnt_type);
	uint32_t memsize = iter_size + cnt_size;
	char *memblock = nb_soft_allocate_mem(memsize);
	nb_container_t *to_delete = (void*) memblock;
	nb_iterator_t *iter = (void*) (memblock + cnt_size);

	nb_container_init(to_delete, cnt_type);

	double r = 1.0 / alpha;
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		double cr = nb_utils2D_get_circumradius(trg->v1->x,
							 trg->v2->x,
							 trg->v3->x);
		if (cr > r)
			nb_container_insert(to_delete, trg);
	}
	nb_iterator_finish(iter);

	while(nb_container_is_not_empty(to_delete)) {
		msh_trg_t *trg = nb_container_delete_first(to_delete);
		mesh_substract_triangle(mesh, trg);
		mtrg_nb_free_mem(mesh, trg);
	}
	nb_container_finish(to_delete);

	nb_soft_free_mem(memsize, memblock);
}

double nb_tessellator2D_get_smallest_ns_alpha_complex(nb_tessellator2D_t *mesh,
					     uint32_t N_vertices,
					     const double *vertices,
					     double alpha_factor)
/* Smallest non-singular */
{
	nb_tessellator2D_get_delaunay(mesh, N_vertices, vertices);
	double rminmax = get_minmax_radius(mesh);
	double alpha = 1.0 / rminmax;
	get_alpha_complex(mesh, alpha * alpha_factor);
	return alpha;
}

static double get_minmax_radius(const nb_tessellator2D_t *mesh)
/* Max(Min(Max())) */
{
	uint32_t N_vtx = nb_tessellator2D_get_N_vtx(mesh);
	uint32_t N_trg = nb_tessellator2D_get_N_trg(mesh);
	uint32_t memsize = N_trg * sizeof(double) +
		N_vtx * sizeof(msh_trg_t*);
	char *memblock = nb_soft_allocate_mem(memsize);
	
	double *max_t = (void*) memblock;
	msh_trg_t **min_v = (void*) (memblock + N_trg * sizeof(double));
	
	mesh_enumerate_vtx((nb_tessellator2D_t*) mesh);
	mesh_enumerate_trg((nb_tessellator2D_t*) mesh);
	get_max_length_x_trg(mesh, max_t);

        get_min_trg_x_vtx(mesh, max_t, min_v);
	
	msh_trg_t *t_minmax = get_max_trg(max_t, N_vtx, min_v);

	double r_minmax = nb_utils2D_get_circumradius(t_minmax->v1->x,
						       t_minmax->v2->x,
						       t_minmax->v3->x);
	nb_soft_free_mem(memsize, memblock);
	return r_minmax;
}

static void get_max_length_x_trg(const nb_tessellator2D_t *mesh, double *max_t)
{
	nb_iterator_t *iter = nb_allocate_on_stack(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		uint32_t id = trg->id;
		max_t[id] = trg_get_max_length(trg);
	}
	nb_iterator_finish(iter);
}

static double trg_get_max_length(msh_trg_t *trg)
{
	double l1 = POW2(trg->v1->x[0] - trg->v2->x[0]) +
		POW2(trg->v1->x[1] - trg->v2->x[1]);
	double l2 = POW2(trg->v2->x[0] - trg->v3->x[0]) +
		POW2(trg->v2->x[1] - trg->v3->x[1]);
	double l3 = POW2(trg->v3->x[0] - trg->v1->x[0]) +
		POW2(trg->v3->x[1] - trg->v1->x[1]);
	return MAX(MAX(l1, l2), l3);
}

static void get_min_trg_x_vtx(const nb_tessellator2D_t *mesh, const double *max_t,
			      msh_trg_t **min_v)
{
	uint32_t N_vtx = nb_tessellator2D_get_N_vtx(mesh);
	memset(min_v, 0, N_vtx * sizeof(*min_v));

	nb_iterator_t *iter = nb_allocate_on_stack(nb_iterator_get_memsize());
	nb_iterator_init(iter);
	nb_iterator_set_container(iter, mesh->ht_trg);
	while (nb_iterator_has_more(iter)) {
		const msh_trg_t *trg = nb_iterator_get_next(iter);
		verify_min_on_vtx(max_t, trg, min_v, mvtx_get_id(trg->v1));
		verify_min_on_vtx(max_t, trg, min_v, mvtx_get_id(trg->v2));
		verify_min_on_vtx(max_t, trg, min_v, mvtx_get_id(trg->v3));
	}
	nb_iterator_finish(iter);	
}

static void verify_min_on_vtx(const double *max_t, const msh_trg_t *trg,
			      msh_trg_t **min_v, uint32_t vid)
{
	uint32_t tid = trg->id;
	if (NULL == min_v[vid]) {
		min_v[vid] = trg;
	} else {
		msh_trg_t *trg_check = min_v[vid];
		uint32_t id_check = trg_check->id;
		if (max_t[tid] < max_t[id_check])
			min_v[vid] = trg;
	}
}

static msh_trg_t *get_max_trg(const double *max_t, uint32_t N_vtx,
			      msh_trg_t **min_v)
{
	msh_trg_t *out = min_v[0];
	uint32_t id = min_v[0]->id;
	double max = max_t[id];
	for (uint32_t i = 1; i < N_vtx; i++) {
		msh_trg_t *trg = min_v[i];
		id = trg->id;
		if (max < max_t[id]) {
			max = max_t[id];
			out = trg;
		}
	}
	return out;
}
