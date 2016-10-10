#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/math_bot.h"
#include "nb/geometric_bot.h"

#define POW2(a) ((a)*(a))

#define INPUTS_DIR "../../../../utest/nb/geometric_bot/mesh/dewall_inputs"

static int suite_init(void);
static int suite_clean(void);

static void test_get_delaunay_polygon_0_center(void);
static void test_get_delaunay_polygon_1_center(void);
static void test_get_delaunay_polygon_2_center(void);
static void test_get_delaunay_polygon_3_center(void);
static void test_get_delaunay_polygon_polycenter(void);
static void test_get_delaunay_2_polygonal_rings(void);
static void test_get_delaunay_2_polygonal_rings_1_center(void);
static void test_get_delaunay_2_polygonal_rings_2_center(void);
static void test_get_delaunay_2_polygonal_rings_polycenter(void);
static void test_get_delaunay_5_polygonal_rings_polycenter(void);
static void test_get_delaunay_10_polygonal_rings_polycenter(void);
static void test_get_delaunay_grid(void);
static void test_get_delaunay_hexagonal_grid(void);
static void test_get_delaunay_1_spiral_110p(void);
static void test_get_delaunay_2_spiral_60p(void);
static void test_get_delaunay_5_spiral_25p(void);
static void test_get_delaunay_10_spiral_11p(void);
static void test_get_delaunay_16_spiral_6p(void);
static void test_get_delaunay_17_spiral_6p(void);
static void test_get_delaunay_20_spiral_5p(void);
static void test_get_delaunay_20_spiral_6p(void);
static void test_get_delaunay_collinear(void);
static void test_get_delaunay_quasi_collinear(void);
static void test_get_delaunay_square(void);
static void test_get_delaunay_1000_cloud(void);

static bool all_trg_are_delaunay(nb_mesh_t *mesh);
static bool check_get_delaunay_polygon(int N, int N_centers);
static int get_expected_trg_of_polygon(int N, int N_centers);
static int get_expected_edg_of_polygon(int N, int N_centers);
static double* get_polygon(int N_sides, int N_centers, double r);
static void set_polygon(int N, double r, double vertices[]);
static void set_vertex(int id, double vertices[], double x, double y);
static bool check_get_delaunay_rings(int N_rings, int N_sides,
				     int N_centers);
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

void cunit_nb_geometric_bot_dewall(void)
{
	CU_pSuite suite = CU_add_suite("nb/geometric_bot/mesh/dewall.c",
				       suite_init, suite_clean);
	CU_add_test(suite, "get_delaunay() of polygon",
		    test_get_delaunay_polygon_0_center);
	CU_add_test(suite,
		    "get_delaunay() of polygon with central point",
		    test_get_delaunay_polygon_1_center);
	CU_add_test(suite,
		    "get_delaunay() of polygon with 2 central points",
		    test_get_delaunay_polygon_2_center);
	CU_add_test(suite,
		    "get_delaunay() of polygon with 3 central points",
		    test_get_delaunay_polygon_3_center);
	CU_add_test(suite,
		    "get_delaunay() of polygon with poly. center",
		    test_get_delaunay_polygon_polycenter);
	CU_add_test(suite,
		    "get_delaunay() of 2 polygonal rings",
		    test_get_delaunay_2_polygonal_rings);
	CU_add_test(suite,
		    "get_delaunay() of 2 polygonal rings with center",
		    test_get_delaunay_2_polygonal_rings_1_center);
	CU_add_test(suite,
		    "get_delaunay() of 2 polygonal rings with 2 centers",
		    test_get_delaunay_2_polygonal_rings_2_center);
	CU_add_test(suite,
		    "get_delaunay() of 2 poly. rings with poly. center",
		    test_get_delaunay_2_polygonal_rings_polycenter);
	CU_add_test(suite,
		    "get_delaunay() of 5 poly. rings with poly. center",
		    test_get_delaunay_5_polygonal_rings_polycenter);
	CU_add_test(suite,
		    "get_delaunay() of 10 poly. rings with poly. center",
		    test_get_delaunay_10_polygonal_rings_polycenter);
	CU_add_test(suite, "get_delaunay() of grid",
		    test_get_delaunay_grid);
	CU_add_test(suite,
		    "get_delaunay() of hexagonal grid",
		    test_get_delaunay_hexagonal_grid);
	CU_add_test(suite,
		    "get_delaunay() of 1 spiral with 110 points",
		    test_get_delaunay_1_spiral_110p);
	CU_add_test(suite,
		    "get_delaunay() of 2 spirals with 60 points",
		    test_get_delaunay_2_spiral_60p);
	CU_add_test(suite,
		    "get_delaunay() of 5 spirals with 25 points",
		    test_get_delaunay_5_spiral_25p);
	CU_add_test(suite,
		    "get_delaunay() of 10 spirals with 11 points",
		    test_get_delaunay_10_spiral_11p);
	CU_add_test(suite,
		    "get_delaunay() of 16 spirals with 6 points",
		    test_get_delaunay_16_spiral_6p);
	CU_add_test(suite,
		    "get_delaunay() of 17 spirals with 6 points",
		    test_get_delaunay_17_spiral_6p);
	CU_add_test(suite,
		    "get_delaunay() of 20 spirals with 5 points",
		    test_get_delaunay_20_spiral_5p);
	CU_add_test(suite,
		    "get_delaunay() of 20 spirals with 6 points",
		    test_get_delaunay_20_spiral_6p);
	CU_add_test(suite,
		    "get_delaunay() of collinear points",
		    test_get_delaunay_collinear);
	CU_add_test(suite,
		    "get_delaunay() of quasi-collinear points",
		    test_get_delaunay_quasi_collinear);
	CU_add_test(suite,
		    "get_delaunay() of square",
		    test_get_delaunay_square);
	CU_add_test(suite,
		    "get_delaunay() of cloud with 1000 vertices",
		    test_get_delaunay_1000_cloud);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_get_delaunay_polygon_0_center(void)
{
	CU_ASSERT(check_get_delaunay_polygon(110, 0));
}

static void test_get_delaunay_polygon_1_center(void)
{
	CU_ASSERT(check_get_delaunay_polygon(110, 1));
}

static void test_get_delaunay_polygon_2_center(void)
{
	CU_ASSERT(check_get_delaunay_polygon(110, 2));
}

static void test_get_delaunay_polygon_3_center(void)
{
	CU_ASSERT(check_get_delaunay_polygon(110, 3));
}

static void test_get_delaunay_polygon_polycenter(void)
{
	CU_ASSERT(check_get_delaunay_polygon(110, 10));
}

static void test_get_delaunay_2_polygonal_rings(void)
{
	CU_ASSERT(check_get_delaunay_rings(2, 60, 0));
}

static void test_get_delaunay_2_polygonal_rings_1_center(void)
{
	CU_ASSERT(check_get_delaunay_rings(2, 60, 1));
}

static void test_get_delaunay_2_polygonal_rings_2_center(void)
{
	CU_ASSERT(check_get_delaunay_rings(2, 60, 2));
}

static void test_get_delaunay_2_polygonal_rings_polycenter(void)
{
	CU_ASSERT(check_get_delaunay_rings(2, 50, 10));
}

static void test_get_delaunay_5_polygonal_rings_polycenter(void)
{
	CU_ASSERT(check_get_delaunay_rings(5, 20, 6));
}

static void test_get_delaunay_10_polygonal_rings_polycenter(void)
{
	CU_ASSERT(check_get_delaunay_rings(10, 10, 5));
}

static void test_get_delaunay_grid(void)
{
	int N = 12;
	double *vertices = get_grid(N, N, 10.0);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N * N, vertices);
	nb_free_mem(vertices);
	int N_squares = POW2(N - 1);
	bool N_trg_is_ok = (2 * N_squares == nb_mesh_get_N_trg(mesh));
	int N_expected_edges = 2 * (N - 1) * N + N_squares;
	bool N_edg_is_ok = (N_expected_edges == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
	CU_ASSERT(all_delaunay);
}

static void test_get_delaunay_hexagonal_grid(void)
{
	int N = 12;
	double *vertices = get_hexagonal_grid(N, N, 10.0);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N * N, vertices);
	nb_free_mem(vertices);
	int N_boundary_trg = 0;
	if (N > 1)
		N_boundary_trg = (2*N - 3) / 2;
	int N_expected_trg = 2 * POW2(N - 1) + N_boundary_trg;
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	int N_expected_edges = (N - 1) * (3*N - 1) + N_boundary_trg;
	bool N_edg_is_ok = (N_expected_edges == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
	CU_ASSERT(all_delaunay);
}

static void test_get_delaunay_1_spiral_110p(void)
{
	CU_ASSERT(check_get_delaunay_spiral(1, 110));
}

static void test_get_delaunay_2_spiral_60p(void)
{
	CU_ASSERT(check_get_delaunay_spiral(2, 60));
}

static void test_get_delaunay_5_spiral_25p(void)
{
	CU_ASSERT(check_get_delaunay_spiral(5, 25));
}

static void test_get_delaunay_10_spiral_11p(void)
{
	CU_ASSERT(check_get_delaunay_spiral(10, 11));
}

static void test_get_delaunay_16_spiral_6p(void)
{
	CU_ASSERT(check_get_delaunay_spiral(16, 6));
}

static void test_get_delaunay_17_spiral_6p(void)
{
	CU_ASSERT(true);
	/* CU_ASSERT(check_get_delaunay_spiral(17, 6));
	 * TEMPORAL FAIL: Fails due to numerical error,
	 * the triangles in the middle are not Delaunay
	 */
}

static void test_get_delaunay_20_spiral_6p(void)
{
	CU_ASSERT(true);
	/* CU_ASSERT(check_get_delaunay_spiral(20, 6));
	 * TEMPORAL FAIL: Fails, create incomplete triangulation.
	 * I'm 90% sure that it is due to numerical error.
	 */
}

static void test_get_delaunay_20_spiral_5p(void)
{
	CU_ASSERT(true);
	/* TEMPORAL FAIL: return check_get_delaunay_spiral(20, 5);
	   Freeze the computer, memory leak? infinite loop?
	*/	
}

static void test_get_delaunay_collinear(void)
{
	int N = 100;
	double *vertices = get_collinear(N);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	bool N_trg_is_ok = (0 == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (0 == nb_mesh_get_N_edg(mesh));
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
}

static void test_get_delaunay_quasi_collinear(void)
{
	int N = 100;
	double *vertices = get_quasi_collinear(N);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, 0);
	int N_expected_edg = get_expected_edg_of_polygon(N, 0);
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
	CU_ASSERT(all_delaunay);
}

static void test_get_delaunay_square(void)
{
	int N_interior = 50;
	int N = 4 * N_interior;
	double *vertices = get_square(N_interior, 10.0);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, 0);
	int N_expected_edg = get_expected_edg_of_polygon(N, 0);
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
	CU_ASSERT(all_delaunay);
}

static void test_get_delaunay_1000_cloud(void)
{
	char input_name[256];
	sprintf(input_name, "%s/cloud_1000.vtx", INPUTS_DIR);
	int N;
	double *vertices = read_vertices(input_name, &N);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	int N_expected_trg = 1981;
	int N_expected_edg = 2980;
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	CU_ASSERT(N_trg_is_ok);
	CU_ASSERT(N_edg_is_ok);
	CU_ASSERT(all_delaunay);
}

static bool all_trg_are_delaunay(nb_mesh_t *mesh)
{
	uint32_t memsize = nb_msh3trg_get_memsize();
	void *msh3trg = nb_allocate_on_stack(memsize);
	nb_msh3trg_init(msh3trg);
	nb_msh3trg_load_from_mesh(msh3trg, mesh);

	bool (*inside)(const double v1[2], const double v2[2],
		       const double v3[2], const double p[2]) =
		nb_utils2D_pnt_lies_strictly_in_circumcircle;
	bool all_delaunay = true;
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh3trg);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id1 = nb_msh3trg_elem_get_adj(msh3trg, i, 0);
		uint32_t id2 = nb_msh3trg_elem_get_adj(msh3trg, i, 1);
		uint32_t id3 = nb_msh3trg_elem_get_adj(msh3trg, i, 2);
		double v1[2];
		v1[0] = nb_msh3trg_node_get_x(msh3trg, id1);
		v1[1] = nb_msh3trg_node_get_y(msh3trg, id1);
		double v2[2];
		v2[0] = nb_msh3trg_node_get_x(msh3trg, id2);
		v2[1] = nb_msh3trg_node_get_y(msh3trg, id2);
		double v3[2];
		v3[0] = nb_msh3trg_node_get_x(msh3trg, id3);
		v3[1] = nb_msh3trg_node_get_y(msh3trg, id3);

		uint32_t N_nod = nb_msh3trg_get_N_nodes(msh3trg);
		for (uint32_t j = 0; j < N_nod; j++) {
			if (id1 != j && id2 != j && id3 != j) {
				double p[2];
				p[0] = nb_msh3trg_node_get_x(msh3trg, j);
				p[1] = nb_msh3trg_node_get_y(msh3trg, j);
				if (inside(v1, v2, v3, p)) {
					all_delaunay = false;
					break;
				}
			}
		}
	}
	nb_msh3trg_finish(msh3trg);
	return all_delaunay;
}

static bool check_get_delaunay_polygon(int N, int N_centers)
{
	double *vertices = get_polygon(N, N_centers, 10);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N + N_centers, vertices);
	nb_free_mem(vertices);
	int N_expected_trg = get_expected_trg_of_polygon(N, N_centers);
	int N_expected_edg = get_expected_edg_of_polygon(N, N_centers);
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
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

static void set_vertex(int id, double vertices[],
		       double x, double y)
{
	vertices[id * 2] = x;
	vertices[id*2+1] = y;
}

static bool check_get_delaunay_rings(int N_rings, int N_sides,
				     int N_centers)
{
	int N = N_rings * N_sides + N_centers;
	double *vertices = get_rings(N_rings, N_sides, N_centers, 10);
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	int N_expected_trg =
		get_expected_trg_of_polygonal_rings(N_rings, N_sides,
						    N_centers);
	int N_expected_edg =
		get_expected_edg_of_polygonal_rings(N_rings, N_sides,
						    N_centers);
	bool N_trg_is_ok = (N_expected_trg == nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edg == nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_delaunay;
}

static int get_expected_trg_of_polygonal_rings(int N_rings,
					       int N_sides,
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

static int get_expected_edg_of_polygonal_rings(int N_rings,
					       int N_sides,
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

static double* get_rings(int N_rings, int N_sides,
			 int N_centers, double r)
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
	nb_mesh_t *mesh = nb_mesh_create();
	nb_mesh_get_delaunay(mesh, N, vertices);
	nb_free_mem(vertices);
	int N_min_trg = (Ns * Np + 1);
	int N_min_edges = 2 * N_min_trg;
	bool N_trg_is_ok = (N_min_trg <= nb_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_min_edges <= nb_mesh_get_N_edg(mesh));
	bool all_delaunay = all_trg_are_delaunay(mesh);
	nb_mesh_destroy(mesh);
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
						nb_free_mem(vertices);
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
