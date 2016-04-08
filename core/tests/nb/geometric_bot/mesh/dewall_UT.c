#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nb/geometric_bot/utils2D.h"
#include "nb/geometric_bot/mesh/elements2D/triangles.h"
#include "nb/geometric_bot/mesh/dewall.h"

#include "test_library.h"
#include "test_add.h"

#define POW2(a) ((a)*(a))

#define INPUTS_DIR "core/tests/nb/geometric_bot/mesh/dewall_UT_inputs"

static bool check_get_delaunay_polygon_0_center(void);
static bool check_get_delaunay_polygon_1_center(void);
static bool check_get_delaunay_polygon_2_center(void);
static bool check_get_delaunay_polygon_3_center(void);
static bool check_get_delaunay_polygon_polycenter(void);
static bool check_get_delaunay_2_polygonal_rings(void);
static bool check_get_delaunay_2_polygonal_rings_1_center(void);
static bool check_get_delaunay_2_polygonal_rings_2_center(void);
static bool check_get_delaunay_2_polygonal_rings_polycenter(void);
static bool check_get_delaunay_5_polygonal_rings_polycenter(void);
static bool check_get_delaunay_10_polygonal_rings_polycenter(void);
static bool check_get_delaunay_grid(void);
static bool check_get_delaunay_hexagonal_grid(void);
static bool check_get_delaunay_1_spiral_110p(void);
static bool check_get_delaunay_2_spiral_60p(void);
static bool check_get_delaunay_5_spiral_25p(void);
static bool check_get_delaunay_10_spiral_11p(void);
static bool check_get_delaunay_16_spiral_6p(void);
static bool check_get_delaunay_17_spiral_6p(void);
static bool check_get_delaunay_20_spiral_5p(void);
static bool check_get_delaunay_20_spiral_6p(void);
static bool check_get_delaunay_collinear(void);
static bool check_get_delaunay_quasi_collinear(void);
static bool check_get_delaunay_square(void);
static bool check_get_delaunay_1000_cloud(void);

static bool all_trg_are_delaunay(vcn_mesh_t *mesh);
static bool check_get_delaunay_polygon(int N, int N_centers);
static int get_expected_trg_of_polygon(int N, int N_centers);
static int get_expected_edg_of_polygon(int N, int N_centers);
static double* get_polygon(int N_sides, int N_centers, double r);
static void set_polygon(int N, double r, double vertices[]);
static void set_vertex(int id, double vertices[], double x, double y);
static bool check_get_delaunay_rings(int N_rings, int N_sides, int N_centers);
static int get_expected_trg_of_polygonal_rings(int N_rings, int N_sides,
					       int N_centers);
static int get_expected_edg_of_polygonal_rings(int N_rings, int N_sides,
					       int N_centers);
static double* get_rings(int N_rings, int N_sides, int N_centers, double r);
static double* get_grid(int N_width, int N_height, double size);
static double* get_hexagonal_grid(int N_width, int N_height, double size);
static bool check_get_delaunay_spiral(int Ns, int Np);
static double* get_spiral(int Ns, int Np);
static void set_spiral(int Np, double init_r, 
		       double init_angle, double vertices[]);
static double* get_collinear(int N);
static void set_line_without_end_point(double *vertices, int N,
				       double x0, double y0,
				       double x1, double y1);
static double* get_quasi_collinear(int N);
static double* get_square(int N_interior, double size);
static double* read_vertices(const char* filename, int* N_vertices);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_get_delaunay_polygon_0_center,
		     "Check get_delaunay() of polygon");
	vcn_test_add(tests_ptr, check_get_delaunay_polygon_1_center,
		     "Check get_delaunay() of polygon with central point");
	vcn_test_add(tests_ptr, check_get_delaunay_polygon_2_center,
		     "Check get_delaunay() of polygon with 2 central points");
	vcn_test_add(tests_ptr, check_get_delaunay_polygon_3_center,
		     "Check get_delaunay() of polygon with 3 central points");
	vcn_test_add(tests_ptr, check_get_delaunay_polygon_polycenter,
		     "Check get_delaunay() of polygon with poly. center");
	vcn_test_add(tests_ptr, check_get_delaunay_2_polygonal_rings,
		     "Check get_delaunay() of 2 polygonal rings");
	vcn_test_add(tests_ptr, check_get_delaunay_2_polygonal_rings_1_center,
		     "Check get_delaunay() of 2 polygonal rings with center");
	vcn_test_add(tests_ptr, check_get_delaunay_2_polygonal_rings_2_center,
		     "Check get_delaunay() of 2 polygonal rings with 2 centers");
	vcn_test_add(tests_ptr, check_get_delaunay_2_polygonal_rings_polycenter,
		     "Check get_delaunay() of 2 poly. rings with poly. center");
	vcn_test_add(tests_ptr, check_get_delaunay_5_polygonal_rings_polycenter,
		     "Check get_delaunay() of 5 poly. rings with poly. center");
	vcn_test_add(tests_ptr, check_get_delaunay_10_polygonal_rings_polycenter,
		     "Check get_delaunay() of 10 poly. rings with poly. center");
	vcn_test_add(tests_ptr, check_get_delaunay_grid,
		     "Check get_delaunay() of grid");
	vcn_test_add(tests_ptr, check_get_delaunay_hexagonal_grid,
		     "Check get_delaunay() of hexagonal grid");
	vcn_test_add(tests_ptr, check_get_delaunay_1_spiral_110p,
		     "Check get_delaunay() of 1 spiral with 110 points");
	vcn_test_add(tests_ptr, check_get_delaunay_2_spiral_60p,
		     "Check get_delaunay() of 2 spirals with 60 points");
	vcn_test_add(tests_ptr, check_get_delaunay_5_spiral_25p,
		     "Check get_delaunay() of 5 spirals with 25 points");
	vcn_test_add(tests_ptr, check_get_delaunay_10_spiral_11p,
		     "Check get_delaunay() of 10 spirals with 11 points");
	vcn_test_add(tests_ptr, check_get_delaunay_16_spiral_6p,
		     "Check get_delaunay() of 16 spirals with 6 points");
	vcn_test_add(tests_ptr, check_get_delaunay_17_spiral_6p,
		     "Check get_delaunay() of 17 spirals with 6 points");
	vcn_test_add(tests_ptr, check_get_delaunay_20_spiral_5p,
		     "Check get_delaunay() of 20 spirals with 5 points");
	vcn_test_add(tests_ptr, check_get_delaunay_20_spiral_6p,
		     "Check get_delaunay() of 20 spirals with 6 points");
	vcn_test_add(tests_ptr, check_get_delaunay_collinear,
		     "Check get_delaunay() of collinear points");
	vcn_test_add(tests_ptr, check_get_delaunay_quasi_collinear,
		     "Check get_delaunay() of quasi-collinear points");
	vcn_test_add(tests_ptr, check_get_delaunay_square,
		     "Check get_delaunay() of square");
	vcn_test_add(tests_ptr, check_get_delaunay_1000_cloud,
		     "Check get_delaunay() of cloud with 1000 vertices");
}

static inline bool check_get_delaunay_polygon_0_center(void)
{
	return check_get_delaunay_polygon(110, 0);
}

static inline bool check_get_delaunay_polygon_1_center(void)
{
	return check_get_delaunay_polygon(110, 1);
}

static inline bool check_get_delaunay_polygon_2_center(void)
{
	return check_get_delaunay_polygon(110, 2);
}

static inline bool check_get_delaunay_polygon_3_center(void)
{
	return check_get_delaunay_polygon(110, 3);
}

static inline bool check_get_delaunay_polygon_polycenter(void)
{
	return check_get_delaunay_polygon(110, 10);
}

static bool check_get_delaunay_2_polygonal_rings(void)
{
	return check_get_delaunay_rings(2, 60, 0);
}

static inline bool check_get_delaunay_2_polygonal_rings_1_center(void)
{
	return check_get_delaunay_rings(2, 60, 1);
}

static inline bool check_get_delaunay_2_polygonal_rings_2_center(void)
{
	return check_get_delaunay_rings(2, 60, 2);
}

static inline bool check_get_delaunay_2_polygonal_rings_polycenter(void)
{
	return check_get_delaunay_rings(2, 50, 10);
}

static inline bool check_get_delaunay_5_polygonal_rings_polycenter(void)
{
	return check_get_delaunay_rings(5, 20, 6);
}

static inline bool check_get_delaunay_10_polygonal_rings_polycenter(void)
{
	return check_get_delaunay_rings(10, 10, 5);
}

static bool check_get_delaunay_grid(void)
{
	int N = 12;
	double *vertices = get_grid(N, N, 10.0);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N * N, vertices);
	free(vertices);
	int N_squares = POW2(N - 1);
	bool N_trg_is_ok = (2 * N_squares == vcn_mesh_get_N_trg(mesh));
	int N_expected_edges = 2 * (N - 1) * N + N_squares;
	bool N_edg_is_ok = (N_expected_edges == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static bool check_get_delaunay_hexagonal_grid(void)
{
	int N = 12;
	double *vertices = get_hexagonal_grid(N, N, 10.0);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N * N, vertices);
	free(vertices);
	int N_boundary_trg = 0;
	if (N > 1)
		N_boundary_trg = (2*N - 3) / 2;
	int N_expected_trg = 2 * POW2(N - 1) + N_boundary_trg;
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	int N_expected_edges = (N - 1) * (3*N - 1) + N_boundary_trg;
	bool N_edg_is_ok = (N_expected_edges == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static inline bool check_get_delaunay_1_spiral_110p(void)
{
	return check_get_delaunay_spiral(1, 110);
}

static inline bool check_get_delaunay_2_spiral_60p(void)
{
	return check_get_delaunay_spiral(2, 60);
}

static inline bool check_get_delaunay_5_spiral_25p(void)
{
	return check_get_delaunay_spiral(5, 25);
}

static inline bool check_get_delaunay_10_spiral_11p(void)
{
	return check_get_delaunay_spiral(10, 11);
}

static inline bool check_get_delaunay_16_spiral_6p(void)
{
	return check_get_delaunay_spiral(16, 6);
}

static inline bool check_get_delaunay_17_spiral_6p(void)
{
	return check_get_delaunay_spiral(17, 6);
	/* TEMPORAL: Fails due to numerical error,
	 * the triangles in the middle are not Delaunay
	 */
}

static inline bool check_get_delaunay_20_spiral_6p(void)
{
	return check_get_delaunay_spiral(20, 6);
	/* TEMPORAL: Fails, create incomplete triangulation.
	 * I'm 90% sure that it is due to numerical error.
	 */
}

static inline bool check_get_delaunay_20_spiral_5p(void)
{
	/* TEMPORAL
	   return check_get_delaunay_spiral(20, 5);
	   Freeze the computer, memory leak? infinite loop?
	*/	
	return false;
}

static bool check_get_delaunay_collinear(void)
{
	int N = 100;
	double *vertices = get_collinear(N);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
	free(vertices);
	bool N_trg_is_ok = (0 == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (0 == vcn_mesh_get_N_edg(mesh));
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok;
}

static bool check_get_delaunay_quasi_collinear(void)
{
	int N = 100;
	double *vertices = get_quasi_collinear(N);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
	free(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, 0);
	int N_expected_edg = get_expected_edg_of_polygon(N, 0);
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static bool check_get_delaunay_square(void)
{
	int N_interior = 50;
	int N = 4 * N_interior;
	double *vertices = get_square(N_interior, 10.0);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
	free(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, 0);
	int N_expected_edg = get_expected_edg_of_polygon(N, 0);
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static bool check_get_delaunay_1000_cloud(void)
{
	char input_name[256];
	sprintf(input_name, "%s/cloud_1000.vtx", INPUTS_DIR);
	int N;
	double *vertices = read_vertices(input_name, &N);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
	free(vertices);
	int N_expected_trg = 1968;
	int N_expected_edg = 2967;
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
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

static bool check_get_delaunay_polygon(int N, int N_centers)
{
	double *vertices = get_polygon(N, N_centers, 10);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N + N_centers, vertices);
	free(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, N_centers);
	int N_expected_edg = get_expected_edg_of_polygon(N, N_centers);
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static int get_expected_trg_of_polygon(int N, int N_centers)
{
	int trg;
	if (3 > N) {
		trg = 0;
	} else {
		if (0 == N_centers)
			trg = N - 2;
		else if (1 == N_centers)
			trg = N;
		else if (2 == N_centers)
			trg = N + 2;
		else
			trg = N + N_centers + 
				get_expected_trg_of_polygon(N_centers, 0);
	}
	return trg;
}

static int get_expected_edg_of_polygon(int N, int N_centers)
{
	int edges;
	if (2 > N) {
		edges = 0;
	} else if (2 == N) {
		edges = 1;
	} else {
		if (0 == N_centers)
			edges = 2 * N - 3;
		else if (1 == N_centers)
			edges = 2 * N;
		else if (2 == N_centers)
			edges = 2 * N + 3;
		else
			edges = 2 * N + N_centers +
				get_expected_edg_of_polygon(N_centers, 0);
	}
	return edges;
}

static double* get_polygon(int N_sides, int N_centers, double r)
{
	int N = N_sides + N_centers;
	double *vertices = malloc(2 * N * sizeof(*vertices));
	set_polygon(N_sides, r, vertices);
	set_polygon(N_centers, r/4, &(vertices[N_sides*2]));
	return vertices;
}

static void set_polygon(int N, double r, double vertices[])
{
	if (1 == N) {
		set_vertex(0, vertices, 0, 0);
	} else if (1 < N) {
		double angle_step = (NB_PI * 2.0) / N;
		for (uint32_t i = 0; i < N; i++) {
			double a = i * angle_step;
			set_vertex(i, vertices, r * cos(a), r * sin(a));
		}
	}
}

static inline void set_vertex(int id, double vertices[], double x, double y)
{
	vertices[id * 2] = x;
	vertices[id*2+1] = y;
}

static bool check_get_delaunay_rings(int N_rings, int N_sides, int N_centers)
{
	int N = N_rings * N_sides + N_centers;
	double *vertices = get_rings(N_rings, N_sides, N_centers, 10);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
	free(vertices);
	int N_expected_trg =
		get_expected_trg_of_polygonal_rings(N_rings, N_sides,
						    N_centers);
	int N_expected_edg =
		get_expected_edg_of_polygonal_rings(N_rings, N_sides,
						    N_centers);
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == vcn_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static int get_expected_trg_of_polygonal_rings(int N_rings, int N_sides,
					       int N_centers)
{
	int trg;
	if (N_rings > 1) {
		trg = 2 * N_sides + 
			get_expected_trg_of_polygonal_rings(N_rings - 1,
							    N_sides,
							    N_centers);
	} else {
		trg = get_expected_trg_of_polygon(N_sides, N_centers);
	}
	return trg;
}

static int get_expected_edg_of_polygonal_rings(int N_rings, int N_sides,
					       int N_centers)
{
	int edges;
	if (N_rings > 1) {
		edges = 3 * N_sides + 
			get_expected_edg_of_polygonal_rings(N_rings - 1,
							    N_sides,
							    N_centers);
	} else {
		edges = get_expected_edg_of_polygon(N_sides, N_centers);
	}
	return edges;
}

static double* get_rings(int N_rings, int N_sides, int N_centers, double r)
{
	int N = N_rings * N_sides + N_centers;
	double *vertices = malloc(2 * N * sizeof(*vertices));
	set_polygon(N_centers, r/4, vertices);
	for (int i = 0; i < N_rings; i++) {
		int id = N_centers + i * N_sides;
		double radii = (i + 1) * r;
		set_polygon(N_sides, radii, &(vertices[id * 2]));
	}
	return vertices;
}


static double* get_grid(int N_width, int N_height, double size)
{
	double *vertices = malloc(2 * N_width * N_height * sizeof(*vertices));
	for (int i = 0; i < N_width; i++) {
		for (int j = 0; j < N_height; j++) {
			int id = i * N_width + j;
			vertices[id * 2] = i * size;
			vertices[id*2+1] = j * size;
		}
	}
	return vertices;
}

static double* get_hexagonal_grid(int N_width, int N_height, double size)
{
	double *vertices = malloc(2 * N_width * N_height * sizeof(*vertices));
	double h_factor = sqrt(0.75);
	for (int i = 0; i < N_width; i++) {
		for (int j = 0; j < N_height; j++) {
			int id = i * N_width + j;
			vertices[id * 2] = j * size + (0.5 * size * (i%2));
			vertices[id*2+1] = i * size * h_factor;
		}
	}
	return vertices;
}

static bool check_get_delaunay_spiral(int Ns, int Np)
{
	int N = Ns * Np + 2;
	double *vertices = get_spiral(Ns, Np);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_delaunay(mesh, N, vertices);
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
		r += 0.75 * r * NB_PHI;
		angle += 1.5 * NB_PI;
	}
	set_vertex(N - 1, vertices, r * cos(angle), r * sin(angle));
	return vertices;
}

static void set_spiral(int Np, double init_r,
		       double init_angle, double vertices[])
{
	double angle_step = (NB_PI * 1.5) / Np;
	double r_step = (0.75 * init_r * NB_PHI) / Np;
	for (uint32_t i = 0; i < Np; i++) {
		double angle = init_angle + i * angle_step;
		double r = init_r + i * r_step;
		vertices[i * 2] = r * cos(angle);
		vertices[i*2+1] = r * sin(angle);
	}
}

static double* get_collinear(int N)
{
	double *vertices = malloc(2 * N * sizeof(*vertices));
	set_line_without_end_point(vertices, N, -N, -N, N, N);
	return vertices;
}

static void set_line_without_end_point(double *vertices, int N,
				       double x0, double y0,
				       double x1, double y1)
{
	double x_step = (x1 - x0) / N;
	double y_step = (y1 - y0) / N;
	for (int i = 0; i < N; i++) {
		vertices[i * 2] = x0 + i * x_step;
		vertices[i*2+1] = y0 + i * y_step;
	}
}

static double* get_quasi_collinear(int N)
{
	double *vertices = malloc(2 * N * sizeof(*vertices));
	for (int i = 0; i < N; i++) {
		vertices[i * 2] = i;
		vertices[i*2+1] = POW2(i);
	}
	return vertices;
}

static double* get_square(int N_interior, double size)
{
	int N = 4 * N_interior;
	double half_size = 0.5 * size;
	double *vertices = malloc(2 * N * sizeof(*vertices));
	set_line_without_end_point(vertices, N_interior,
				   -half_size, -half_size,
				   half_size, -half_size);
	set_line_without_end_point(&(vertices[N_interior * 2]), N_interior, 
				   half_size, -half_size,
				   half_size, half_size);
	set_line_without_end_point(&(vertices[2*N_interior * 2]), N_interior, 
				   half_size, half_size,
				   -half_size, half_size);
	set_line_without_end_point(&(vertices[3*N_interior * 2]), N_interior, 
				   -half_size, half_size,
				   -half_size, -half_size);
	return vertices;
}

static double* read_vertices(const char* filename, int* N_vertices)
{
	double *vertices = NULL;
	*N_vertices = 0;
	FILE* fp = fopen(filename, "r");
	if (NULL != fp) {
		if (1 == fscanf(fp, "%i", N_vertices)) {
			if (0 < *N_vertices) {
				vertices = malloc(*N_vertices * 2 * sizeof(*vertices));
				for (uint32_t i = 0; i < *N_vertices; i++) {
					if (2 != fscanf(fp, "%lf %lf",
							&(vertices[i * 2]),
							&(vertices[i*2+1]))) {
						free(vertices);
						vertices = NULL;
						*N_vertices = 0;
						break;
					}
				}
			}
		}
		fclose(fp);
	}
	return vertices;
}
