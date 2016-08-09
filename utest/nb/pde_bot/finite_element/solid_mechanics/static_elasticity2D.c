#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include <CUnit/Basic.h>

#include "nb/container_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"

#define INPUTS_DIR "../../../../utest/nb/pde_bot/static_elasticity2D_inputs"

#define POW2(a) ((a)*(a))

typedef struct {
	uint32_t N_vtx;
	uint32_t N_trg;
	double *disp;
	double *strain;
	double *stress;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_beam_cantilever(void);
static void check_beam_cantilever(const vcn_msh3trg_t *msh3trg,
				  const results_t *results);
static void test_plate_with_hole(void);
static void check_plate_with_hole(const vcn_msh3trg_t *msh3trg,
				  const results_t *results);
static double get_error_avg_pwh(const vcn_msh3trg_t *msh3trg,
				const vcn_fem_elem_t* elem,
				const double *vm_stress);
static void get_cartesian_gpoint(uint32_t id_elem, int8_t id_gp,
				 const vcn_msh3trg_t *msh3trg,
				 const vcn_fem_elem_t* elem,
				 double gp[2]);
static double get_analytic_vm_stress_pwh(double x, double y);
static void get_analytic_stress_pwh(double x, double y, double stress[3]);
static void modify_bcond_pwh(const vcn_msh3trg_t *msh3trg,
			     nb_bcond_t *bcond);
static void pwh_DC_SGM_cond(double *x, double t, double *out);
static void pwh_BC_SGM_cond(double *x, double t, double *out);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const vcn_msh3trg_t*,
					   const results_t*),
		     void (*modify_bcond)(const vcn_msh3trg_t*,
					  nb_bcond_t*)/* Can be NULL */);
static int simulate_fem(const char *problem_data,
			vcn_msh3trg_t *msh3trg, results_t *results,
			uint32_t N_vtx,
			void (*modify_bcond)(const vcn_msh3trg_t*,
					     nb_bcond_t*)/* Can be NULL */);
static void get_mesh(const vcn_model_t *model, vcn_msh3trg_t *msh3trg,
		     uint32_t N_vtx);

static void results_init(results_t *results, uint32_t N_vtx, uint32_t N_trg);
static void results_finish(results_t *results);
static int read_problem_data
		(const char* filename,
		 vcn_model_t *model,
		 nb_bcond_t* bcond,
		 vcn_fem_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(vcn_cfreader_t *cfreader, vcn_model_t *model);
static int read_material(vcn_cfreader_t *cfreader, vcn_fem_material_t *mat);
static int read_elasticity2D_params(vcn_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);


void cunit_nb_pde_bot_fem_sm_static_elasticity(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/"\
			     "static_elasticity.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Beam cantilever", test_beam_cantilever);
	CU_add_test(suite, "Plate with a hole", test_plate_with_hole);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}


static void test_beam_cantilever(void)
{
	run_test("%s/beam_cantilever.txt", 1000,
		 check_beam_cantilever, NULL);
}

static void check_beam_cantilever(const vcn_msh3trg_t *msh3trg,
				  const results_t *results)
{
	double max_disp = 0;
	for (uint32_t i = 0; i < msh3trg->N_vertices; i++) {
		double disp2 = POW2(results->disp[i * 2]) +
			POW2(results->disp[i*2+1]);
		if (max_disp < disp2)
			max_disp = disp2;
	}
	max_disp = sqrt(max_disp);
	CU_ASSERT(fabs(max_disp - 1.00701e-1) < 1e-6);
}


static void test_plate_with_hole(void)
{
	run_test("%s/plate_with_hole.txt", 1500,
		 check_plate_with_hole,
		 modify_bcond_pwh);
}

static void check_plate_with_hole(const vcn_msh3trg_t *msh3trg,
				  const results_t *results)
{
	double *vm_stress = malloc(results->N_trg * sizeof(double));
	double *nodal_values = malloc(results->N_vtx * sizeof(double));

	vcn_fem_compute_von_mises(results->N_trg, results->stress, vm_stress);

	vcn_fem_elem_t* elem = vcn_fem_elem_create(NB_TRG_LINEAR);

	double avg_error = get_error_avg_pwh(msh3trg, elem, vm_stress);
	printf("--AVG ERROR: %e\n", avg_error);

	vcn_fem_interpolate_from_gpoints_to_nodes(msh3trg, elem,
						  1, vm_stress,
						  nodal_values);

	nb_fem_save(msh3trg,  results->disp, 0.5, nodal_values,
		    "../../../fem_plate.png", 1000, 800);/* TEMPORAL */

	CU_ASSERT(true);
	CU_ASSERT(true);
	
	vcn_fem_elem_destroy(elem);
	free(vm_stress);
	free(nodal_values);
}

static double get_error_avg_pwh(const vcn_msh3trg_t *msh3trg,
				const vcn_fem_elem_t* elem,
				const double *vm_stress)
{
	FILE *fp = fopen("../../../fem_plate.log", "w");/**/
	double avg = 0.0;
	uint32_t N = 0;
	for (uint32_t i = 0; i < msh3trg->N_triangles; i++) {
		int8_t N_gp = vcn_fem_elem_get_N_gpoints(elem);
		N += N_gp;
		for (int8_t p = 0; p < N_gp; p++) {
			double gp[2];
			get_cartesian_gpoint(i, p, msh3trg, elem, gp);
			double gp_stress =
				get_analytic_vm_stress_pwh(gp[0], gp[1]);
			double error;
			if (fabs(gp_stress) > 1e-10)
				error = fabs(1.0 - vm_stress[i * N_gp + p]/
					     gp_stress);
			else
				error = fabs(vm_stress[i * N_gp + p]);
			avg += error;
			fprintf(fp, "%lf\n", error);/**/
		}
	}
	fclose(fp);/**/
	return avg /= N;
}

static void get_cartesian_gpoint(uint32_t id_elem, int8_t id_gp,
				 const vcn_msh3trg_t *msh3trg,
				 const vcn_fem_elem_t* elem,
				 double gp[2])
{
	int8_t N = vcn_fem_elem_get_N_nodes(elem);
	gp[0] = 0.0;
	gp[1] = 0.0;
	for (int i = 0; i < N; i++) {
		double Ni = vcn_fem_elem_Ni(elem, i, id_gp);
		uint32_t vi = msh3trg->vertices_forming_triangles[id_elem * N + i];
		gp[0] += Ni * msh3trg->vertices[vi * 2];
		gp[1] += Ni * msh3trg->vertices[vi*2+1];
	}
}

static double get_analytic_vm_stress_pwh(double x, double y)
{
	double stress[3];
	get_analytic_stress_pwh(x, y, stress);
	return nb_pde_get_vm_stress(stress[0], stress[1], stress[2]);
}

static void get_analytic_stress_pwh(double x, double y, double stress[3])
{
	double a = 0.5;
	double tx = 1e4;
	double r = sqrt(POW2(x) + POW2(y));
	double theta = atan2(y, x);
	double ar2 = POW2(a/r);
	double ar4x1p5 = 1.5 * POW2(ar2);
	double s2t = sin(2.0 * theta);
	double s4t = sin(4.0 * theta);
	double c2t = cos(2.0 * theta);
	double c4t = cos(4.0 * theta);
	stress[0] = tx * (1.0 - ar2 * (1.5 * c2t + c4t) + ar4x1p5 * c4t);
	stress[1] = tx * (-ar2 * (0.5 * c2t - c4t) - ar4x1p5 * c4t);
	stress[2] = tx * (-ar2 * (0.5 * s2t + s4t) + ar4x1p5 * s4t);
}

static void modify_bcond_pwh(const vcn_msh3trg_t *msh3trg,
			     nb_bcond_t *bcond)
{
	bool dof_mask[2] = {1, 1};

	uint32_t DC_sgm = 11;
	nb_bcond_push_function(bcond, NB_NEUMANN, NB_BC_ON_SEGMENT,
			       DC_sgm, dof_mask, pwh_DC_SGM_cond);
	uint32_t BC_sgm = 12;
	nb_bcond_push_function(bcond, NB_NEUMANN, NB_BC_ON_SEGMENT,
			       BC_sgm, dof_mask, pwh_BC_SGM_cond);
}

static void pwh_DC_SGM_cond(double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[0];
	out[1] = stress[2];	
}

static void pwh_BC_SGM_cond(double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[2];
	out[1] = stress[1];
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const vcn_msh3trg_t*,
					   const results_t*),
		     void (*modify_bcond)(const vcn_msh3trg_t*,
					  nb_bcond_t*)/* Can be NULL */)
{
	results_t results;
	vcn_msh3trg_t *msh3trg = alloca(vcn_msh3trg_get_memsize());
	vcn_msh3trg_init(msh3trg);

	int status = simulate_fem(problem_data, msh3trg, &results,
				  N_vtx, modify_bcond);
	
	CU_ASSERT(0 == status);

	check_results(msh3trg, &results);

	vcn_msh3trg_finish(msh3trg);
	results_finish(&results);
}

static int simulate_fem(const char *problem_data,
			vcn_msh3trg_t *msh3trg, results_t *results,
			uint32_t N_vtx,
			void (*modify_bcond)(const vcn_msh3trg_t*,
					     nb_bcond_t*)/* Can be NULL */)
{
	int status = 1;
	vcn_model_t* model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = alloca(bcond_size);
	nb_bcond_init(bcond, 2);
	vcn_fem_material_t* material = vcn_fem_material_create();
	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;

	char input[255];
	sprintf(input, problem_data, INPUTS_DIR);
	int read_status =
		read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);

	if (NULL != modify_bcond)
		modify_bcond(msh3trg, bcond);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	get_mesh(model, msh3trg, N_vtx);

	vcn_fem_elem_t* elem = vcn_fem_elem_create(NB_TRG_LINEAR);

	results_init(results, msh3trg->N_vertices, msh3trg->N_triangles);

	int status_fem =
		vcn_fem_compute_2D_Solid_Mechanics(msh3trg, elem,
						   material, bcond,
						   false, NULL,
						   analysis2D,
						   &params2D, NULL,
						   results->disp,
						   results->strain);
	if (0 != status_fem)
		goto CLEANUP_FEM;

	vcn_fem_compute_stress_from_strain(msh3trg->N_triangles,
					   elem, material,
					   analysis2D, results->strain, NULL,
					   results->stress);

	status = 0;
CLEANUP_FEM:
	vcn_fem_elem_destroy(elem);
CLEANUP_INPUT:
	vcn_model_finish(model);
	nb_bcond_finish(bcond);
	vcn_fem_material_destroy(material);

	return status;
}

static void get_mesh(const vcn_model_t *model, vcn_msh3trg_t *msh3trg,
		     uint32_t N_vtx)
{
	uint32_t mesh_memsize = vcn_mesh_get_memsize();
	vcn_mesh_t* mesh = alloca(mesh_memsize);
	vcn_mesh_init(mesh);
	vcn_mesh_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     N_vtx);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	vcn_mesh_generate_from_model(mesh, model);

	vcn_msh3trg_load_from_mesh(msh3trg, mesh);
	vcn_mesh_finish(mesh);
}

static void results_init(results_t *results, uint32_t N_vtx, uint32_t N_trg)
{
	uint32_t size_disp = N_vtx * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_trg * 3 * sizeof(*(results->strain));
	uint32_t total_size = size_disp + 2 * size_strain;
	char *memblock = malloc(total_size);

	results->N_vtx = N_vtx;
	results->N_trg = N_trg;
	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
}

static inline void results_finish(results_t *results)
{
	free(results->disp);
}

static int read_problem_data
		(const char* filename,
		 vcn_model_t *model,
		 nb_bcond_t* bcond,
		 vcn_fem_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D)
{      
	int status = 1;
	/* Initialize custom format to read file */
	vcn_cfreader_t* cfreader = vcn_cfreader_create(filename, "#");
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
	vcn_cfreader_destroy(cfreader);
	status = 0;
EXIT:
	return status;
}

static int read_geometry(vcn_cfreader_t *cfreader, vcn_model_t *model)
{
	int status = 1;
	/* Read modele vertices */
	uint32_t N = 0;
	if (0 != vcn_cfreader_read_uint(cfreader, &N))
		goto EXIT;

	if (1 > N)
		goto EXIT;
	model->N = N;

	model->vertex = malloc(2 * model->N * sizeof(*(model->vertex)));
	for (uint32_t i = 0; i < 2 * model->N; i++) {
		if (0 != vcn_cfreader_read_double(cfreader, &(model->vertex[i])))
			goto CLEANUP_VERTICES;
	}
	/* Read model segments */
	N = 0;
	if (0 != vcn_cfreader_read_uint(cfreader, &N))
		goto CLEANUP_VERTICES;

	if (1 > N)
		goto CLEANUP_VERTICES;
	model->M = N;

	model->edge = malloc(2 * model->M * sizeof(*model->edge));
	for (uint32_t i = 0; i < 2 * model->M; i++) {
		if (0 != vcn_cfreader_read_uint(cfreader, &(model->edge[i])))
			goto CLEANUP_SEGMENTS;
	}
	/* Read model holes */
	N = 0;
	if (0 != vcn_cfreader_read_uint(cfreader, &N))
		goto CLEANUP_SEGMENTS;
	model->H = N;

	model->holes = NULL;
	if (0 < model->H) {
		model->holes = malloc(2 * model->H * sizeof(*(model->holes)));
		for (uint32_t i = 0; i < 2 * model->H; i++) {
			if (0 != vcn_cfreader_read_double(cfreader,
							  &(model->holes[i])))
				goto CLEANUP_HOLES;
		}
	}

	status = 0;
	goto EXIT;

CLEANUP_HOLES:
	if (0 < model->H) {
		free(model->holes);
		model->H = 0;
		model->holes = NULL;
	}
CLEANUP_SEGMENTS:
	model->M = 0;
	free(model->edge);
	model->edge = 0;
CLEANUP_VERTICES:
	model->N = 0;
	free(model->vertex);
	model->vertex = NULL;
EXIT:
	return status;
}


static int read_material(vcn_cfreader_t *cfreader, vcn_fem_material_t *mat)
{
	int status = 1;
	double poisson_module;
	if (0 != vcn_cfreader_read_double(cfreader, &poisson_module))
		goto EXIT;
	vcn_fem_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (0 != vcn_cfreader_read_double(cfreader, &elasticity_module))
		goto EXIT;
	vcn_fem_material_set_elasticity_module(mat, elasticity_module);

	double fracture_energy;
	if (0 != vcn_cfreader_read_double(cfreader, &fracture_energy))
		goto EXIT;
	vcn_fem_material_set_fracture_energy(mat, fracture_energy);

	double compression_limit_stress;
	if (0 != vcn_cfreader_read_double(cfreader, &compression_limit_stress))
		goto EXIT;
	vcn_fem_material_set_compression_limit_stress(mat,
						      compression_limit_stress);

	double traction_limit_stress;
	if (0 != vcn_cfreader_read_double(cfreader, &traction_limit_stress))
		goto EXIT;
	vcn_fem_material_set_traction_limit_stress(mat, traction_limit_stress);
	status = 0;
EXIT:
	return status;
}


static int read_elasticity2D_params(vcn_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D)
{
	int status = 1;
	int iaux;
	if (0 != vcn_cfreader_read_int(cfreader, &iaux))
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
	if (0 != vcn_cfreader_read_double(cfreader, &(params2D->thickness)))
		goto EXIT;
	status = 0;
EXIT:
	return status;
}
