#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "cunit/Basic.h"

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

#define INPUTS_DIR "../utest/sources/nb/pde_bot/damage_inputs"

#define POW2(a) ((a)*(a))
#define CHECK_ZERO(a) ((fabs(a)<1e-25)?1:(a))
#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	uint32_t N_faces;
	uint32_t N_elems;
	double *disp;
	double *strain;
  	double *stress;
	double *damage;
	char *boundary_mask;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_mode_I(void);
static void check_mode_I(const void *mesh,
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
static int read_geometry(nb_cfreader_t *cfreader, nb_model_t *model);
static int read_material(nb_cfreader_t *cfreader, nb_material_t *mat);
static int read_elasticity2D_params(nb_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static void show_cvfa_error_msg(int cvfa_status);

void cunit_nb_pde_bot_cvfa_sm_static_damage_phase_field(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/" \
			     "static_elasticity.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Mode I Phase field", test_mode_I);
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
	run_test("%s/Mode_I_3point_bending.txt", 5000, NB_QUAD,
		 check_mode_I);
}

static void check_mode_I(const void *mesh,
			 const results_t *results)
{
	CU_ASSERT(true);/* TEMPORAL */
}


static void TEMPORAL1(nb_mesh2D_t *mesh, results_t *results)
{
	nb_mesh2D_export_draw(mesh, "./mesh.eps", 1000, 800,
			      NB_NULL, NB_NULL, NULL, true);/* TEMPORAL */

	nb_cvfa_draw_integration_mesh(mesh, "./CVFA_alpha_x.eps",/*T*/
				      1000, 800);              /* TEMPORAL */

	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	double *disp = malloc(N_elems * sizeof(*disp));

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	double *disp_nodes = malloc(N_nodes * sizeof(*disp_nodes));

	nb_mesh2D_distort_with_field(mesh, NB_ELEMENT, results->disp, 20.0);

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
	uint32_t memsize = N_nodes * (4 * sizeof(double) + sizeof(uint16_t));
	char *memblock = nb_soft_allocate_mem(memsize);
	double *stress = (void*) memblock;
	double *vm_stress = (void*) (memblock + N_nodes * 3 * sizeof(double));
	uint16_t *counter = (void*) (memblock + N_nodes * 4 * sizeof(double));

	memset(stress, 0, 3 * N_nodes * sizeof(*stress));
	memset(counter, 0, N_nodes * sizeof(*counter));
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		uint32_t v1 = nb_mesh2D_edge_get_1n(mesh, i);
		uint32_t v2 = nb_mesh2D_edge_get_2n(mesh, i);
		if (!(results->boundary_mask[i])) {
			stress[v1 * 3] += results->stress[i * 3];
			stress[v1*3+1] += results->stress[i*3+1];
			stress[v1*3+2] += results->stress[i*3+2];
			counter[v1] += 1;

			stress[v2 * 3] += results->stress[i * 3];
			stress[v2*3+1] += results->stress[i*3+1];
			stress[v2*3+2] += results->stress[i*3+2];
			counter[v2] += 1;
		}
	}
	for (uint32_t i = 0; i < N_nodes; i++) {
		stress[i * 3] /= counter[i];
		stress[i*3+1] /= counter[i];
		stress[i*3+2] /= counter[i];
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

	double max_dmg = 0;
	for (uint32_t i = 0; i < N_faces; i++) {
		if (max_dmg < results->damage[i])
			max_dmg = results->damage[i];
	}
	printf("====> MAX DAMAGE: %e\n", max_dmg);/* TEMPORAL */
	for (uint32_t i = 0; i < N_faces && 0; i++)
		results->damage[i] /= max_dmg;
	nb_mesh2D_export_draw(mesh, "./CVFA_dmg.png", 1000, 800,
			      NB_FACE, NB_FIELD,
			      results->damage, true);/* TEMPORAL */

	nb_soft_free_mem(memsize, memblock);
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_mesh2D_type mesh_type,
		     void (*check_results)(const void*,
					   const results_t*))
{
	results_t results;
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

	int status_cvfa = nb_cvfa_compute_2D_damage_phase_field(mesh, material,
								bcond,
								false, NULL,
								analysis2D,
								&params2D,
								results->disp,
								results->strain,
								results->damage,
								results->boundary_mask);

	if (0 != status_cvfa) {
		show_cvfa_error_msg(status_cvfa);
		goto CLEANUP;
	}

	nb_cvfa_compute_stress_from_strain(mesh,  material, analysis2D,
					   results->strain,
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
}

static inline void results_finish(results_t *results)
{
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
	int status = 1;
	/* Initialize custom format to read file */
	nb_cfreader_t* cfreader = nb_cfreader_create(filename, "#");
	if (NULL == cfreader) {
		printf("\nERROR: Can not open file %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_geometry(cfreader, model)) {
		printf("\nERROR: Geometry contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != nb_bcond_read(bcond, cfreader)) {
		printf("\nERROR: Boundary C. contain errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_material(cfreader, mat)) {
		printf("\nERROR: Material contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_elasticity2D_params(cfreader, analysis2D, params2D)) {
		printf("\nERROR: Reading numerical params in %s.\n",
		       filename);
		goto EXIT;
	}
	status = 0;
EXIT:
	nb_cfreader_destroy(cfreader);
	return status;
}

static int read_geometry(nb_cfreader_t *cfreader, nb_model_t *model)
{
	int status = 1;
	/* Read modele vertices */
	uint32_t N = 0;
	if (0 != nb_cfreader_read_uint(cfreader, &N))
		goto EXIT;

	if (1 > N)
		goto EXIT;
	model->N = N;

	model->vertex = nb_allocate_mem(2 * model->N * sizeof(*(model->vertex)));
	for (uint32_t i = 0; i < 2 * model->N; i++) {
		if (0 != nb_cfreader_read_double(cfreader, &(model->vertex[i])))
			goto CLEANUP_VERTICES;
	}
	/* Read model segments */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfreader, &N))
		goto CLEANUP_VERTICES;

	if (1 > N)
		goto CLEANUP_VERTICES;
	model->M = N;

	model->edge = nb_allocate_mem(2 * model->M * sizeof(*model->edge));
	for (uint32_t i = 0; i < 2 * model->M; i++) {
		if (0 != nb_cfreader_read_uint(cfreader, &(model->edge[i])))
			goto CLEANUP_SEGMENTS;
	}
	/* Read model holes */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfreader, &N))
		goto CLEANUP_SEGMENTS;
	model->H = N;

	model->holes = NULL;
	if (0 < model->H) {
		model->holes = nb_allocate_mem(2 * model->H * sizeof(*(model->holes)));
		for (uint32_t i = 0; i < 2 * model->H; i++) {
			if (0 != nb_cfreader_read_double(cfreader,
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


static int read_material(nb_cfreader_t *cfreader, nb_material_t *mat)
{
	int status = 1;
	double poisson_module;
	if (0 != nb_cfreader_read_double(cfreader, &poisson_module))
		goto EXIT;
	nb_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (0 != nb_cfreader_read_double(cfreader, &elasticity_module))
		goto EXIT;
	nb_material_set_elasticity_module(mat, elasticity_module);

	double fracture_energy;
	if (0 != nb_cfreader_read_double(cfreader, &fracture_energy))
		goto EXIT;
	nb_material_set_fracture_energy(mat, fracture_energy);

	double compression_limit_stress;
	if (0 != nb_cfreader_read_double(cfreader, &compression_limit_stress))
		goto EXIT;
	nb_material_set_compression_limit_stress(mat,
						      compression_limit_stress);

	double traction_limit_stress;
	if (0 != nb_cfreader_read_double(cfreader, &traction_limit_stress))
		goto EXIT;
	nb_material_set_traction_limit_stress(mat, traction_limit_stress);
	status = 0;
EXIT:
	return status;
}


static int read_elasticity2D_params(nb_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D)
{
	int status = 1;
	int iaux;
	if (0 != nb_cfreader_read_int(cfreader, &iaux))
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
	if (0 != nb_cfreader_read_double(cfreader, &(params2D->thickness)))
		goto EXIT;
	status = 0;
EXIT:
	return status;
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
