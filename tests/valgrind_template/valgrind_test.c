/* valgrind --tool=memcheck 
 *          --track-origins=yes
 *          --leak-check=full
 *          --show-leak-kinds=all */

#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "vcn_bots.h"

static bool all_trg_are_delaunay(vcn_mesh_t *mesh);
static bool check_get_delaunay_spiral(int Ns, int Np);
static double* get_spiral(int Ns, int Np);
static void set_spiral(int Np, double init_r, 
		       double init_angle, double vertices[]);
static void set_vertex(int id, double vertices[], double x, double y);

int main()
{
	check_get_delaunay_spiral(16, 6);
	return 0;
}

static bool all_trg_are_delaunay(vcn_mesh_t *mesh)
{
	vcn_msh3trg_t *msh3trg = vcn_mesh_get_msh3trg(mesh, true, true,
						      false, false, false);
	bool (*inside)(const double v1[2], const double v2[2],
		       const double v3[2], const double p[2]) =
		vcn_utils2D_pnt_lies_strictly_in_circumcircle;
	bool all_delaunay = true;
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		uint32_t id1 = vcn_msh3trg_get_1st_vtx_id_of_trg(msh3trg, i);
		uint32_t id2 = vcn_msh3trg_get_2nd_vtx_id_of_trg(msh3trg, i);
		uint32_t id3 = vcn_msh3trg_get_3rd_vtx_id_of_trg(msh3trg, i);
		double *v1 = vcn_msh3trg_view_vtx(msh3trg, id1);
		double *v2 = vcn_msh3trg_view_vtx(msh3trg, id2);
		double *v3 = vcn_msh3trg_view_vtx(msh3trg, id3);
		for (uint32_t j = 0; j < msh3trg->N_vertices; j++) {
			if (id1 != j && id2 != j && id3 != j) {
				double *p = &(msh3trg->vertices[j*2]);
				if (inside(v1, v2, v3, p)) {
					all_delaunay = false;
					break;
				}
			}
		}
	}
	vcn_msh3trg_destroy(msh3trg);
	return all_delaunay;
}

static bool check_get_delaunay_spiral(int Ns, int Np)
{
	int N = Ns * Np + 2;
	double *vertices = get_spiral(Ns, Np);
	vcn_mesh_t *mesh = vcn_dewall_get_delaunay(N, vertices);
	free(vertices);
	int N_min_trg = (Ns * Np + 1);
	int N_min_edges = 2 * N_min_trg;
	bool N_trg_is_ok = (N_min_trg <= vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_min_edges <= vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static double* get_spiral(int Ns, int Np)
{
	int N = Ns * Np + 2;
	double *vertices = malloc(2 * N * sizeof(*vertices));
	double r = 1.0;
	double angle = 0.0;
	set_vertex(0, vertices, 0, 0);
	for (int i = 0; i < Ns; i++)  {
		int first_id = 1 + i * Np;
		set_spiral(Np, r, angle, &(vertices[first_id * 2]));
		r += 0.75 * r * VCN_PHI;
		angle += 1.5 * VCN_PI;
	}
	set_vertex(N - 1, vertices, r * cos(angle), r * sin(angle));
	return vertices;
}

static void set_spiral(int Np, double init_r,
		       double init_angle, double vertices[])
{
	double angle_step = (VCN_PI * 1.5) / Np;
	double r_step = (0.75 * init_r * VCN_PHI) / Np;
	for (uint32_t i = 0; i < Np; i++) {
		double angle = init_angle + i * angle_step;
		double r = init_r + i * r_step;
		vertices[i * 2] = r * cos(angle);
		vertices[i*2+1] = r * sin(angle);
	}
}

static inline void set_vertex(int id, double vertices[], double x, double y)
{
	vertices[id * 2] = x;
	vertices[id*2+1] = y;
}
