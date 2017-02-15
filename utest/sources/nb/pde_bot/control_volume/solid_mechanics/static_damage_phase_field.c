#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "cunit/Basic.h"

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/io_bot.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

#define INPUTS_DIR "../utest/sources/nb/pde_bot/damage_inputs"
#define OUTPUT_DIR "dmg_tmp"

#define POW2(a) ((a)*(a))
#define CHECK_ZERO(a) ((fabs(a)<1e-25)?1:(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct {
	uint32_t N_faces;
	uint32_t N_elems;
	double *disp;
	double *strain;
  	double *stress;
	double *damage;
	char *boundary_mask;
	bool allocated;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_mode_I(void);
static void test_mode_I_perfored_strip(void);
static void test_mode_II(void);
static void test_mode_II_asym_notched_3point_bending(void);
static void test_count_steps_in_results(void);
static void test_draw_results(void);

static void check_mode_I(const void *mesh,
			 const results_t *results);
static void check_mode_I_perfored_strip(const void *mesh,
					const results_t *results);
static void check_mode_II(const void *mesh,
			  const results_t *results);
static void check_mode_II_asym_notched_3point_bending(const void *mesh,
						      const results_t *results);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_mesh2D_type mesh_type,
		     void (*check_results)(const void*,
					   const results_t*));
static int simulate(const char *problem_data, nb_mesh2D_t *mesh,
		    results_t *results, uint32_t N_vtx);
static void get_mesh(const nb_model_t *model, void *mesh,
		     uint32_t N_vtx);
static void results_init(results_t *results, uint32_t N_faces,
			 uint32_t N_elems);
static void results_finish(results_t *results);
static int read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(nb_cfreader_t *cfr, nb_model_t *model);
static int read_material(nb_cfreader_t *cfr, nb_material_t *mat);
static int read_elasticity2D_params(nb_cfreader_t *cfr,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static void show_drawing_progress(float prog);
static void stdout_show_progress(const char *task, float prog);
static void show_cvfa_error_msg(int cvfa_status);

void cunit_nb_pde_bot_cvfa_sm_static_damage_phase_field(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/" \
			     "static_elasticity.c",
			     suite_init, suite_clean);
	//CU_add_test(suite, "Mode I Phase field", test_mode_I);
	CU_add_test(suite, "Mode I Perfored Strip under tension",
		    test_mode_I_perfored_strip);
	//CU_add_test(suite, "Mode II Phase field", test_mode_II);
	//CU_add_test(suite, "Mode II Asym notched 3 point bending",
	//	    test_mode_II_asym_notched_3point_bending);
	CU_add_test(suite, "Count steps in results",
		    test_count_steps_in_results);
	CU_add_test(suite, "Drawing results", test_draw_results);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_mode_I(void)
{
	run_test("%s/Mode_I_3point_bending.txt", 11000, NB_POLY,
		 check_mode_I);
}

static void check_mode_I(const void *mesh,
			 const results_t *results)
{
	CU_ASSERT(true);/* TEMPORAL */
}

static void test_mode_I_perfored_strip(void)
{
	run_test("%s/Mode_I_perfored_strip_under_tension.txt", 11000, NB_POLY,
		 check_mode_I_perfored_strip);
}

static void check_mode_I_perfored_strip(const void *mesh,
					const results_t *results)
{
	CU_ASSERT(true);/* TEMPORAL */
}

static void test_mode_II(void)
{
	run_test("%s/Mode_II_4point_bending.txt", 6000, NB_POLY,
		 check_mode_II);
}

static void check_mode_II(const void *mesh,
			 const results_t *results)
{
	CU_ASSERT(true);/* TEMPORAL */
}

static void test_mode_II_asym_notched_3point_bending(void)
{
	run_test("%s/Mode_II_Asym_notched_3point_bending.txt", 10000, NB_POLY,
		 check_mode_II_asym_notched_3point_bending);
}

static void check_mode_II_asym_notched_3point_bending(const void *mesh,
						      const results_t *results)
{
	CU_ASSERT(true);/* TEMPORAL */
}


static void TEMPORAL1(nb_mesh2D_t *mesh, results_t *results)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	double *disp = malloc(N_elems * sizeof(*disp));

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	double *disp_nodes = malloc(N_nodes * sizeof(*disp_nodes));

	double box[4];
	nb_mesh2D_get_enveloping_box(mesh, box);
	double max_dist = 0.05 * MAX(box[2]-box[0], box[3]-box[1]);
	nb_mesh2D_distort_with_field(mesh, NB_ELEMENT, results->disp, max_dist);

	for (uint32_t i = 0; i < N_elems; i++)
		disp[i] = results->disp[i*2];

	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, disp, disp_nodes);
	nb_mesh2D_export_draw(mesh, "./CVFA_dx.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 disp_nodes, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_elems; i++)
		disp[i] = results->disp[i*2+1];

	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, disp, disp_nodes);
	nb_mesh2D_export_draw(mesh, "./CVFA_dy.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 disp_nodes, true);/* TEMPORAL */

	nb_free_mem(disp);
	nb_free_mem(disp_nodes);
}

static void TEMPORAL2(nb_mesh2D_t *mesh, results_t *results)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	uint32_t memsize = N_nodes * (7 * sizeof(double) + sizeof(uint16_t));
	char *memblock = nb_soft_allocate_mem(memsize);
	double *strain = (void*) memblock;
	double *stress = (void*) (memblock + N_nodes * 3 * sizeof(double));
	double *vm_stress = (void*) (memblock + N_nodes * 6 * sizeof(double));
	uint16_t *counter = (void*) (memblock + N_nodes * 7 * sizeof(double));

	memset(stress, 0, 3 * N_nodes * sizeof(*stress));
	memset(strain, 0, 3 * N_nodes * sizeof(*strain));
	memset(counter, 0, N_nodes * sizeof(*counter));
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		uint32_t v1 = nb_mesh2D_edge_get_1n(mesh, i);
		uint32_t v2 = nb_mesh2D_edge_get_2n(mesh, i);
		if (!(results->boundary_mask[i])) {
			stress[v1 * 3] += results->stress[i * 3];
			stress[v1*3+1] += results->stress[i*3+1];
			stress[v1*3+2] += results->stress[i*3+2];
			strain[v1 * 3] += results->strain[i * 3];
			strain[v1*3+1] += results->strain[i*3+1];
			strain[v1*3+2] += results->strain[i*3+2];
			counter[v1] += 1;

			stress[v2 * 3] += results->stress[i * 3];
			stress[v2*3+1] += results->stress[i*3+1];
			stress[v2*3+2] += results->stress[i*3+2];
			strain[v2 * 3] += results->strain[i * 3];
			strain[v2*3+1] += results->strain[i*3+1];
			strain[v2*3+2] += results->strain[i*3+2];
			counter[v2] += 1;
		}
	}
	for (uint32_t i = 0; i < N_nodes; i++) {
		if (counter[i] > 0) {
			stress[i * 3] /= counter[i];
			stress[i*3+1] /= counter[i];
			stress[i*3+2] /= counter[i];
			strain[i * 3] /= counter[i];
			strain[i*3+1] /= counter[i];
			strain[i*3+2] /= counter[i];
		}
	}

	for (uint32_t i = 0; i < N_nodes; i++) {
		vm_stress[i] = nb_pde_get_vm_stress(stress[i * 3],
						    stress[i*3+1],
						    stress[i*3+2]);
	}

	nb_mesh2D_export_draw(mesh, "./CVFA_vm.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3];

	nb_mesh2D_export_draw(mesh, "./CVFA_Sxx.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3+1];

	nb_mesh2D_export_draw(mesh, "./CVFA_Syy.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3+2];

	nb_mesh2D_export_draw(mesh, "./CVFA_Sxy.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++) {
		double tr = strain[i * 3] + strain[i*3+1];
		double norm2 = POW2(strain[i * 3]) +
		  2 * POW2(0.5 * strain[i*3+2]) + POW2(strain[i*3+1]);
		double mu = 8.333333e+09;
		double lambda = 4.166667e+09;
		vm_stress[i] = mu * norm2 + 0.5 * lambda * POW2(tr);
	}

	nb_mesh2D_export_draw(mesh, "./CVFA_energy.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */

	nb_soft_free_mem(memsize, memblock);
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_mesh2D_type mesh_type,
		     void (*check_results)(const void*,
					   const results_t*))
{
	results_t results;
	results.allocated = false;

	uint32_t mesh_memsize = nb_mesh2D_get_memsize(mesh_type);
	nb_mesh2D_t *mesh = nb_soft_allocate_mem(mesh_memsize);
	nb_mesh2D_init(mesh, mesh_type);

	int status = simulate(problem_data, mesh, &results, N_vtx);
	
	CU_ASSERT(0 == status);

	check_results(mesh, &results);

	TEMPORAL2(mesh, &results);
	TEMPORAL1(mesh, &results);

	nb_mesh2D_finish(mesh);
	nb_soft_free_mem(mesh_memsize, mesh);
	results_finish(&results);
}

static int simulate(const char *problem_data, nb_mesh2D_t *mesh,
		    results_t *results, uint32_t N_vtx)
{
	int status = 1;
	nb_model_t* model = nb_allocate_on_stack(nb_model_get_memsize());
	nb_model_init(model);
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = nb_allocate_on_stack(bcond_size);
	nb_bcond_init(bcond, 2);
	nb_material_t* material = nb_material_create();
	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;

	char input[255];
	sprintf(input, problem_data, INPUTS_DIR);
	int read_status =
		read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);

	if (0 != read_status)
		goto CLEANUP;

	get_mesh(model, mesh, N_vtx);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	results_init(results, N_faces, N_elems);

	sprintf(input, "%s/mesh.nbt", OUTPUT_DIR);
	nb_mesh2D_save_nbt(mesh, input);

	int status_cvfa =
		nb_cvfa_compute_2D_damage_phase_field(mesh, material,
						      bcond,
						      false, NULL,
						      analysis2D,
						      &params2D,
						      OUTPUT_DIR,
						      results->disp,
						      results->strain,
						      results->damage,
						      results->boundary_mask);

	if (0 != status_cvfa) {
		show_cvfa_error_msg(status_cvfa);
		goto CLEANUP;
	}

	nb_cvfa_compute_stress_from_damage_and_strain(mesh, material,
						      analysis2D,
						      results->strain,
						      results->damage,
						      results->stress);

	status = 0;
CLEANUP:
	nb_model_finish(model);
	nb_bcond_finish(bcond);
	nb_material_destroy(material);

	return status;
}

static void get_mesh(const nb_model_t *model, void *mesh,
		     uint32_t N_vtx)
{
	uint32_t t2d_memsize = nb_tessellator2D_get_memsize();
	nb_tessellator2D_t* t2d = nb_allocate_on_stack(t2d_memsize);
	nb_tessellator2D_init(t2d);
	nb_tessellator2D_set_size_constraint(t2d,
					     NB_MESH_SIZE_CONSTRAINT_MAX_TRG,
					     N_vtx);
	nb_tessellator2D_set_geometric_constraint(t2d,
						  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
						  NB_GEOMETRIC_TOL);
	nb_tessellator2D_generate_from_model(t2d, model);

	nb_mesh2D_load_from_tessellator2D(mesh, t2d);
	nb_tessellator2D_finish(t2d);

	nb_mesh2D_centroid_iteration(mesh, 500, NULL, NULL);
	nb_mesh2D_export_draw(mesh, "./mesh.png", 1000, 800,
			      NB_NULL, NB_NULL, NULL, true);/* TEMPORAL */
	nb_cvfa_draw_integration_mesh(mesh, "./CVFA_alpha_x.eps",/*T*/
				      1000, 800);              /* TEMPORAL */
}

static void results_init(results_t *results, uint32_t N_faces,
			 uint32_t N_elems)
{
	uint32_t size_disp = N_elems * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_faces * 3 * sizeof(*(results->strain));
	uint32_t size_damage = N_faces * sizeof(*(results->damage));
	uint32_t size_mask = N_faces * sizeof(*(results->boundary_mask));
	uint32_t total_size = size_disp + 2 * size_strain +
		size_damage + size_mask;
	char *memblock = nb_allocate_mem(total_size);

	results->N_faces = N_faces;
	results->N_elems = N_elems;
	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
	results->damage = (void*)(memblock + size_disp + 2 * size_strain);
	results->boundary_mask = (void*)(memblock + size_disp +
					 2 * size_strain + size_damage);
	results->allocated = true;
}

static inline void results_finish(results_t *results)
{
	if (results->allocated)
		nb_free_mem(results->disp);
}

static int read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D)
{
	nb_cfreader_t* cfr = nb_cfreader_create();
	nb_cfreader_add_line_comment_token(cfr, "#");
	int status = nb_cfreader_open_file(cfr, filename);
	if (0 != status) {
		printf("\nERROR: Can not open file %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_geometry(cfr, model)) {
		printf("\nERROR: Geometry contains errors in %s.\n",
		       filename);
		goto CLOSE;
	}
	if (0 != nb_bcond_read(bcond, cfr)) {
		printf("\nERROR: Boundary C. contain errors in %s.\n",
		       filename);
		goto CLOSE;
	}
	if (0 != read_material(cfr, mat)) {
		printf("\nERROR: Material contains errors in %s.\n",
		       filename);
		goto CLOSE;
	}
	if (0 != read_elasticity2D_params(cfr, analysis2D, params2D)) {
		printf("\nERROR: Reading numerical params in %s.\n",
		       filename);
		goto CLOSE;
	}
	status = 0;
CLOSE:
	nb_cfreader_close_file(cfr);
EXIT:
	nb_cfreader_destroy(cfr);
	return status;
}

static int read_geometry(nb_cfreader_t *cfr, nb_model_t *model)
{
	int status = 1;
	/* Read modele vertices */
	uint32_t N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto EXIT;

	if (1 > N)
		goto EXIT;
	model->N = N;

	model->vertex = nb_allocate_mem(2 * model->N * sizeof(*(model->vertex)));
	for (uint32_t i = 0; i < 2 * model->N; i++) {
		if (0 != nb_cfreader_read_double(cfr, &(model->vertex[i])))
			goto CLEANUP_VERTICES;
	}
	/* Read model segments */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto CLEANUP_VERTICES;

	if (1 > N)
		goto CLEANUP_VERTICES;
	model->M = N;

	model->edge = nb_allocate_mem(2 * model->M * sizeof(*model->edge));
	for (uint32_t i = 0; i < 2 * model->M; i++) {
		if (0 != nb_cfreader_read_uint(cfr, &(model->edge[i])))
			goto CLEANUP_SEGMENTS;
	}
	/* Read model holes */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto CLEANUP_SEGMENTS;
	model->H = N;

	model->holes = NULL;
	if (0 < model->H) {
		model->holes = nb_allocate_mem(2 * model->H *
					       sizeof(*(model->holes)));
		for (uint32_t i = 0; i < 2 * model->H; i++) {
			if (0 != nb_cfreader_read_double(cfr,
							  &(model->holes[i])))
				goto CLEANUP_HOLES;
		}
	}

	status = 0;
	goto EXIT;

CLEANUP_HOLES:
	if (0 < model->H) {
		nb_free_mem(model->holes);
		model->H = 0;
		model->holes = NULL;
	}
CLEANUP_SEGMENTS:
	model->M = 0;
	nb_free_mem(model->edge);
	model->edge = 0;
CLEANUP_VERTICES:
	model->N = 0;
	nb_free_mem(model->vertex);
	model->vertex = NULL;
EXIT:
	return status;
}


static int read_material(nb_cfreader_t *cfr, nb_material_t *mat)
{
	int status = 1;
	double poisson_module;
	if (0 != nb_cfreader_read_double(cfr, &poisson_module))
		goto EXIT;
	nb_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (0 != nb_cfreader_read_double(cfr, &elasticity_module))
		goto EXIT;
	nb_material_set_elasticity_module(mat, elasticity_module);

	double energy_release_rate;
	if (0 != nb_cfreader_read_double(cfr, &energy_release_rate))
		goto EXIT;
	nb_material_set_energy_release_rate(mat, energy_release_rate);

	double damage_length_scale;
	if (0 != nb_cfreader_read_double(cfr, &damage_length_scale))
		goto EXIT;
	nb_material_set_damage_length_scale(mat, damage_length_scale);
	status = 0;
EXIT:
	return status;
}


static int read_elasticity2D_params(nb_cfreader_t *cfr,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D)
{
	int status = 1;
	int iaux;
	if (0 != nb_cfreader_read_int(cfr, &iaux))
		goto EXIT;
	
	switch (iaux) {
	case 0:
		*analysis2D = NB_PLANE_STRESS;
		break;
	case 1:
		*analysis2D = NB_PLANE_STRAIN;
		break;
	case 2:
		*analysis2D = NB_SOLID_OF_REVOLUTION;
		break;
	default:
		*analysis2D = NB_PLANE_STRESS;
	}

	/* FIX: Usable only for plane stress */
	if (0 != nb_cfreader_read_double(cfr, &(params2D->thickness)))
		goto EXIT;

	status = 0;
EXIT:
	return status;
}

static void test_count_steps_in_results(void)
{
	char name[100];
	sprintf(name, "%s/Mode_I_results.nbt", INPUTS_DIR);
	uint32_t N_steps;
	int status = nb_mesh2D_field_get_N_steps(name, &N_steps);
	CU_ASSERT(0 == status);
	CU_ASSERT(14 == N_steps);
}

static void test_draw_results(void)
{
	int status = nb_mesh2D_field_read_and_draw(OUTPUT_DIR, OUTPUT_DIR,
						   show_drawing_progress);
	CU_ASSERT(0 == status);/* TEMPORAL */
}

static void show_drawing_progress(float prog)
{
	stdout_show_progress("Drawing CFVA", prog);
}

static void stdout_show_progress(const char *task, float prog)
{
	int prog_x100 = (int)(prog * 100 + 0.5);
	printf("%s [", task);
	for (int i = 0; i < 40; i++) {
		if (prog_x100 > 2.49 * (i+1))
			printf("#");
		else
			printf(" ");
	}
	printf("] %i %%\r", prog_x100);
	fflush(stdout);
}
static void show_cvfa_error_msg(int cvfa_status)
{
	switch (cvfa_status) {
	case 1:
		printf("-- NB_ERROR: CVFA assembly fails.\n");
		break;
	case 2:
		printf("-- NB_ERROR: CVFA solver fails.\n");
		break;
	default:
		printf("-- NB_ERROR: CVFA fails.\n");
		break;
	}
}
