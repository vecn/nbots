#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/math_bot.h"
#include "nb/geometric_bot.h"

#define POW2(a) ((a)*(a))

typedef struct {
	uint32_t N_vtx;
	uint32_t N_sgm;
	double *vertices;
	uint32_t *segments;
} input_t;

static int suite_init(void);
static int suite_clean(void);

static void test_get_cdelaunay_twist_square(void);
static void test_get_cdelaunay_twist_polygon_5(void);
static void test_get_cdelaunay_twist_polygon_50(void);
static void test_get_cdelaunay_twist_polygon_100(void);
static void test_get_cdelaunay_60_strips(void);
static void test_get_cdelaunay_centipede_5(void);
static void test_get_cdelaunay_centipede_50(void);

static bool check_get_cdelaunay_twist_polygon(uint32_t N_sides);
static void set_twist_polygon(input_t *input, uint32_t N_sides, double size);
static int get_expected_trg_of_polygon(int N, int N_centers);
static int get_expected_edg_of_polygon(int N, int N_centers);
static void set_polygon(int N, double r, double vertices[]);
static void set_polygon_bissections(int N_sides,
				    const double vtx_polygon[],
				    int N_biss, double vertices[]);
static void set_vertex(int id, double vertices[], double x, double y);
static bool check_get_cdelaunay_strips(int N);
static void set_strips(input_t *input, uint32_t N_strips, double size);
static void input_set_sgm(input_t *input, uint32_t id,
			  uint32_t v1, uint32_t v2);
static void input_alloc(input_t *input);
static void input_clear(input_t *input);
static bool all_trg_are_cdelaunay(vcn_mesh_t *mesh, input_t *input);
static bool is_not_constrained(const double v1[2], const double v2[2],
			       const double v3[2], const double p[2],
			       const input_t *const input);
static bool input_constrains(const input_t *const input,
			     const double v1[2], const double v2[2],
			     const double v3[2], const double p[2]);
static bool intersects(const double a1[2], const double a2[2],
		       const double b1[2], const double b2[2]);
static bool check_get_cdelaunay_centipede(int N_pairs);
static void set_centipede(input_t *input, uint32_t N_pairs);

void cunit_nb_geometric_bot_constrained_delaunay(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/geometric_bot/mesh/constrained_delaunay.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "get_constrained_delaunay() of twisted square",
		    test_get_cdelaunay_twist_square);
	CU_add_test(suite, "get_constrained_delaunay() of twisted pentagon",
		    test_get_cdelaunay_twist_polygon_5);
	CU_add_test(suite, "get_constrained_delaunay() tw. poly. of 50 vtx",
		    test_get_cdelaunay_twist_polygon_50);
	CU_add_test(suite, "get_constrained_delaunay() tw. poly. of 100 vtx",
		    test_get_cdelaunay_twist_polygon_100);
	CU_add_test(suite,
		    "get_constrained_delaunay() of 60 parallel strips",
		    test_get_cdelaunay_60_strips);
	CU_add_test(suite,
		    "get_constrained_delaunay() of 5 pair centipede",
		    test_get_cdelaunay_centipede_5);
	CU_add_test(suite,
		    "get_constrained_delaunay() of 50 pair centipede",
		    test_get_cdelaunay_centipede_50);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_get_cdelaunay_twist_square(void)
{
	CU_ASSERT(check_get_cdelaunay_twist_polygon(4));
}

static void test_get_cdelaunay_twist_polygon_5(void)
{
	CU_ASSERT(check_get_cdelaunay_twist_polygon(5));
}

static void test_get_cdelaunay_twist_polygon_50(void)
{
	CU_ASSERT(check_get_cdelaunay_twist_polygon(50));
}

static void test_get_cdelaunay_twist_polygon_100(void)
{
	CU_ASSERT(check_get_cdelaunay_twist_polygon(100));
}

static void test_get_cdelaunay_60_strips(void)
{
	CU_ASSERT(check_get_cdelaunay_strips(60));
}

static void test_get_cdelaunay_centipede_5(void)
{
	CU_ASSERT(check_get_cdelaunay_centipede(5));
}

static void test_get_cdelaunay_centipede_50(void)
{
	CU_ASSERT(check_get_cdelaunay_centipede(50));
}

static bool check_get_cdelaunay_twist_polygon(uint32_t N_sides)
{
	input_t input;
	set_twist_polygon(&input, N_sides, 20);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_constrained_delaunay(mesh, input.N_vtx, input.vertices,
					  input.N_sgm, input.segments);
	int N_expected_trg = 3 * N_sides +
		get_expected_trg_of_polygon(N_sides, 0);
	int N_expected_edges = 5 * N_sides +
		get_expected_edg_of_polygon(N_sides, 0);
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edges == vcn_mesh_get_N_edg(mesh));
	bool all_cdelaunay = all_trg_are_cdelaunay(mesh, &input);
	input_clear(&input);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_cdelaunay;
}

static void set_twist_polygon(input_t *input, uint32_t N_sides, double size)
{
	input->N_vtx = 3 * N_sides;
	input->N_sgm = N_sides;
	input_alloc(input);	
	set_polygon(N_sides, size/3.0, input->vertices);
	set_polygon(N_sides, size, &(input->vertices[N_sides * 2]));
	set_polygon_bissections(N_sides, &(input->vertices[N_sides * 2]),
				1, &(input->vertices[(2 * N_sides) * 2]));
	for (int i = 0; i < N_sides; i++) {
		int v1 = (i + 1) % N_sides;
		int v2 = N_sides + i;
		input_set_sgm(input, i, v1, v2);
	}
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

static void set_polygon_bissections(int N_sides, const double vtx_polygon[],
				    int N_biss, double vertices[])
{
	if (1 < N_sides) {
		double step = 1.0 / (1.0 + N_biss);
		for (int i = 0; i < N_sides; i++) {
			int p1 = i;
			int p2 = (i + 1) % N_sides;
			for (int j = 0; j < N_biss; j++) {
				int id = i * N_biss + j;
				double w = (j + 1) * step;
				vertices[id * 2] = w * vtx_polygon[p2 * 2] +
					(1.0 - w) * vtx_polygon[p1 * 2];
				vertices[id*2+1] = w * vtx_polygon[p2*2+1] +
					(1.0 - w) * vtx_polygon[p1*2+1];
			}
		}
	}
}

static inline void set_vertex(int id, double vertices[], double x, double y)
{
	vertices[id * 2] = x;
	vertices[id*2+1] = y;
}

static bool check_get_cdelaunay_strips(int N)
{
	input_t input;
	set_strips(&input, N, 10);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_constrained_delaunay(mesh, input.N_vtx, input.vertices,
					  input.N_sgm, input.segments);
	int N_boundary_trg = 0;
	if (N > 1)
		N_boundary_trg = (2*N - 3) / 2;
	int N_expected_trg = 2 * (N - 1) + N_boundary_trg;
	int N_expected_edges = 3 * (N - 1) + N + N_boundary_trg;
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edges == vcn_mesh_get_N_edg(mesh));
	bool all_cdelaunay = all_trg_are_cdelaunay(mesh, &input);
	input_clear(&input);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_cdelaunay;
}

static void set_strips(input_t *input, uint32_t N_strips, double size)
{
	input->N_vtx = 2 * N_strips;
	input->N_sgm = N_strips;
	input_alloc(input);
	for (int i = 0; i < N_strips; i++) {
		double x0 = (i % 2) * (size / 2);
		double y = i * (size / 10);
		set_vertex(i * 2, input->vertices, x0, y);
		set_vertex(i*2+1, input->vertices, x0 + size, y);
		input_set_sgm(input, i, i * 2, i*2+1);
	}
}

static inline void input_set_sgm(input_t *input, uint32_t id,
				 uint32_t v1, uint32_t v2)
{
	input->segments[id * 2] = v1;
	input->segments[id*2+1] = v2;
}

static void input_alloc(input_t *input)
{
	if (0 < input->N_vtx) {
		input->vertices = malloc(2 * input->N_vtx * 
					 sizeof(*(input->vertices)));
	}
	if (0 < input->N_sgm) {
		input->segments = malloc(2 * input->N_sgm * 
					 sizeof(*(input->segments)));
	}
}

static void input_clear(input_t *input)
{
	if (0 < input->N_vtx)
		free (input->vertices);
	if (0 < input->N_sgm)
		free(input->segments);
}

static bool all_trg_are_cdelaunay(vcn_mesh_t *mesh, input_t *input)
{
	uint32_t memsize = nb_msh3trg_get_memsize();
	void *msh3trg = alloca(memsize);
	nb_msh3trg_init(msh3trg);
	nb_msh3trg_load_from_mesh(msh3trg, mesh);

	bool all_delaunay = true;
	uint32_t N_elems = nb_msh3trg_get_N_elems(msh3trg);
	for (uint32_t i = 0; i < N_elems; i++) {
		uint32_t id1 = nb_msh3trg_elem_get_adj(msh3trg, i, 0);
		uint32_t id2 = nb_msh3trg_elem_get_adj(msh3trg, i, 1);
		uint32_t id3 = nb_msh3trg_elem_get_adj(msh3trg, i, 2);
		
		double v1[2];
		v1[0] = nb_msh3trg_get_x_node(msh3trg, id1);
		v1[1] = nb_msh3trg_get_y_node(msh3trg, id1);
		double v2[2];
		v2[0] = nb_msh3trg_get_x_node(msh3trg, id2);
		v2[1] = nb_msh3trg_get_y_node(msh3trg, id2);
		double v3[2];
		v3[0] = nb_msh3trg_get_x_node(msh3trg, id3);
		v3[1] = nb_msh3trg_get_y_node(msh3trg, id3);

		uint32_t N_nodes = nb_msh3trg_get_N_nodes(msh3trg);
		for (uint32_t j = 0; j < N_nodes; j++) {
			if (id1 != j && id2 != j && id3 != j) {
				double p[2];
				p[0] = nb_msh3trg_get_x_node(msh3trg, id3);
				p[1] = nb_msh3trg_get_y_node(msh3trg, id3);
				if (is_not_constrained(v1, v2, v3, p, input)) {
					all_delaunay = false;
					break;
				}
			}
		}
	}
	nb_msh3trg_finish(msh3trg);
	return all_delaunay;
}

static bool is_not_constrained(const double v1[2], const double v2[2],
			       const double v3[2], const double p[2],
			       const input_t *const input)
{
	bool (*inside)(const double v1[2], const double v2[2],
		       const double v3[2], const double p[2]) =
		vcn_utils2D_pnt_lies_strictly_in_circumcircle;
	bool is_constrained = true;
	if (inside(v1, v2, v3, p))
		is_constrained = input_constrains(input, v1, v2, v3, p);
	return ! is_constrained;
}

static bool input_constrains(const input_t *const input,
			     const double v1[2], const double v2[2],
			     const double v3[2], const double p[2])
			     
{
	bool constrains = false;
	for (int i = 0; i < input->N_sgm; i++) {
		int id1 = input->segments[i * 2];
		int id2 = input->segments[i*2+1];
		double *s1 = &(input->vertices[id1 * 2]);
		double *s2 = &(input->vertices[id2 * 2]);
		if (intersects(v1, p, s1, s2)) {
			constrains = true;
			break;
		} else if (intersects(v2, p, s1, s2)) {
			constrains = true;
			break;
		} else if (intersects(v3, p, s1, s2)) {
			constrains = true;
			break;
		}
	}
	return constrains;
}

static inline bool intersects(const double a1[2], const double a2[2],
			      const double b1[2], const double b2[2])
{
	return (NB_INTERSECTED == 
		vcn_utils2D_are_sgm_intersected(a1, a2, b1, b2, NULL));
}

static bool check_get_cdelaunay_centipede(int N_pairs)
{
	input_t input;
	set_centipede(&input, N_pairs);
	vcn_mesh_t *mesh = vcn_mesh_create();
	vcn_mesh_get_constrained_delaunay(mesh, input.N_vtx, input.vertices,
					  input.N_sgm, input.segments);
	int N_expected_trg = 2 * N_pairs;
	int N_expected_edges = 4 * N_pairs + 1;
	bool N_trg_is_ok = (N_expected_trg == vcn_mesh_get_N_trg(mesh));
	bool N_edg_is_ok = (N_expected_edges == vcn_mesh_get_N_edg(mesh));
	bool all_cdelaunay = all_trg_are_cdelaunay(mesh, &input);
	input_clear(&input);
	vcn_mesh_destroy(mesh);
	return N_trg_is_ok && N_edg_is_ok && all_cdelaunay;
}

static void set_centipede(input_t *input, uint32_t N_pairs)
{
	input->N_vtx = 2 * N_pairs + 2;
	input->N_sgm = 1;
	input_alloc(input);
	
	double half_size = 20 + (N_pairs - 1)/2.0;
	set_vertex(0, input->vertices, -half_size, 0);
	set_vertex(1, input->vertices, half_size, 0);
	for (uint32_t i = 0; i < N_pairs; i++) {
		uint32_t id1 = 2 + (i * 2);
		uint32_t id2 = 2 + (i*2+1);
		double x = (double) i - (N_pairs - 1)/2.0;
		set_vertex(id1, input->vertices, x, 1.0);
		set_vertex(id2, input->vertices, x, -1.0);		
	}
	input_set_sgm(input, 0, 0, 1);
}
