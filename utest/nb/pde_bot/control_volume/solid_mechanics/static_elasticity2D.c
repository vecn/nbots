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
	uint32_t N_nodes;
	uint32_t N_elems;
	double *disp;
	double *strain;
	double *stress;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_beam_cantilever(void);
static void check_beam_cantilever(const void *part,
				  const results_t *results);
static void test_plate_with_hole(void);
static void check_plate_with_hole(const void *part,
				  const results_t *results);
static double get_error_avg_pwh(const void *part,
				const double *vm_stress);
static double get_analytic_vm_stress_pwh(double x, double y);
static void get_analytic_stress_pwh(double x, double y, double stress[3]);
static void modify_bcond_pwh(const void *part,
			     nb_bcond_t *bcond);
static void pwh_DC_SGM_cond(const double *x, double t, double *out);
static void pwh_BC_SGM_cond(const double *x, double t, double *out);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_partition_type part_type,
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */);
static int simulate(const char *problem_data,
		    nb_partition_t *part, results_t *results,
		    uint32_t N_vtx,
		    void (*modify_bcond)(const void*,
					     nb_bcond_t*)/* Can be NULL */);
static void get_mesh(const vcn_model_t *model, void *part,
		     uint32_t N_vtx);

static void results_init(results_t *results, uint32_t N_nodes, uint32_t N_elems);
static void results_finish(results_t *results);
static int read_problem_data
		(const char* filename,
		 vcn_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(vcn_cfreader_t *cfreader, vcn_model_t *model);
static int read_material(vcn_cfreader_t *cfreader, nb_material_t *mat);
static int read_elasticity2D_params(vcn_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static void show_cvfa_error_msg(int cvfa_status);

void cunit_nb_pde_bot_cvfa_sm_static_elasticity(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/"\
			     "static_elasticity.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Beam cantilever", test_beam_cantilever);
	//CU_add_test(suite, "Plate with a hole", test_plate_with_hole);
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
	uint8_t N = 6;/********++ TEMPORAL ************/
	double xi[12] = {0,0,
			 1.5, 0,
			 0.464, 1.4265,
			 -1.2135, 0.8817,
			 -1.2135, -0.8817,
			 0.4635, -1.4265};
	double eval[6];
	double grad[12];
	FILE *fp =  fopen("../../../phi_1.txt", "w");
	fprintf(fp, "#g(x) = x\n");
	uint16_t N_grid = 50;
	double step = 6.0/(N_grid-1);
	for (uint16_t i = 0; i < N_grid; i++) {
		double x[2];
		x[0] = i * step - 3;
		for (uint16_t j = 0; j < N_grid; j++) {
			x[1] = j * step - 3;
			nb_nonpolynomial_simple_eval(N, 2, xi,
						     x, eval);
			nb_nonpolynomial_simple_eval_grad(N, 2, xi,
							  x, grad);
			fprintf(fp, "%e %e", x[0], x[1]);
			for (uint8_t k = 0; k < N; k++)
				fprintf(fp, " %e %e %e", eval[k],
					grad[k*2], grad[k*2+1]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);/********++ TEMPORAL ************/

	run_test("%s/beam_cantilever.txt", 1000, NB_TRIAN,
		 check_beam_cantilever, NULL);
}

static void check_beam_cantilever(const void *part,
				  const results_t *results)
{
	double max_disp = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double disp2 = POW2(results->disp[i * 2]) +
			POW2(results->disp[i*2+1]);
		if (max_disp < disp2)
			max_disp = disp2;
	}
	max_disp = sqrt(max_disp);
	printf("--- MAX DISP: %e\n", max_disp); /* TEMPORAL */
	CU_ASSERT(fabs(max_disp - 1.000109e-1) < 1e-6);
}


static void test_plate_with_hole(void)
{
	run_test("%s/plate_with_hole.txt", 1000, NB_QUAD,
		 check_plate_with_hole,
		 modify_bcond_pwh);
}

static void check_plate_with_hole(const void *part,
				  const results_t *results)
{
	double *vm_stress = malloc(results->N_elems * sizeof(double));
	double *nodal_values = malloc(results->N_nodes * sizeof(double));

	nb_pde_compute_von_mises(results->N_elems, results->stress, vm_stress);

	double avg_error = get_error_avg_pwh(part, vm_stress);

	printf("--- AVG ERROR: %e\n", avg_error); /* TEMPORAL */
	CU_ASSERT(avg_error < 9.7e-3);

	free(vm_stress);
	free(nodal_values);
}

static double get_error_avg_pwh(const void *part,
				const double *vm_stress)
{
	double avg = 0.0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		double x = nb_partition_elem_get_x(part, i);
		double y = nb_partition_elem_get_y(part, i);
		double stress =
			get_analytic_vm_stress_pwh(x, y);
		double error;
		if (fabs(stress) > 1e-10)
			error = fabs(1.0 - vm_stress[i] / stress);
		else
			error = fabs(vm_stress[i]);
		avg += error;
	}
	return avg /= N_elems;
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

static void modify_bcond_pwh(const void *part,
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

static void pwh_DC_SGM_cond(const double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[0];
	out[1] = stress[2];	
}

static void pwh_BC_SGM_cond(const double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[2];
	out[1] = stress[1];
}

static void TEMPORAL1(nb_partition_t *part, results_t *results)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	double *total_disp = malloc(N_elems * sizeof(*total_disp));

	for (uint32_t i = 0; i < N_elems; i++) {
		total_disp[i] = sqrt(POW2(results->disp[i*2]) +
				     POW2(results->disp[i*2+1]));
	}

	uint32_t N_nodes = nb_partition_get_N_nodes(part);
	double *disp_nodes = malloc(N_nodes * sizeof(*disp_nodes));

	nb_partition_extrapolate_elems_to_nodes(part, 1, total_disp,
						disp_nodes);
	nb_partition_distort_with_field(part, NB_ELEMENT, results->disp, 0.5);

	nb_partition_export_draw(part, "../../../CVFA.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 disp_nodes, true);/* TEMPORAL */
	free(total_disp);
	free(disp_nodes);
}

static void TEMPORAL2(nb_partition_t *part, results_t *results)
{
	uint32_t N_elems = nb_partition_get_N_elems(part);
	double *vm_stress = malloc(N_elems * sizeof(*vm_stress));

	for (uint32_t i = 0; i < N_elems; i++) {
		vm_stress[i] = nb_pde_get_vm_stress(results->stress[i * 3],
						    results->stress[i*3+1],
						    results->stress[i*3+2]);
	}

	uint32_t N_nodes = nb_partition_get_N_nodes(part);
	double *vm_nodes = malloc(N_nodes * sizeof(*vm_nodes));

	nb_partition_extrapolate_elems_to_nodes(part, 1, vm_stress,
						vm_nodes);

	nb_partition_export_draw(part, "../../../CVFA_stress_elems.png",
				 1000, 800,
				 NB_ELEMENT, NB_FIELD,
				 vm_stress, true);/* TEMPORAL */
	nb_partition_export_draw(part, "../../../CVFA_stress.png", 1000, 800,
				 NB_NODE, NB_FIELD,
				 vm_nodes, true);/* TEMPORAL */
	free(vm_stress);
	free(vm_nodes);
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_partition_type part_type,
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */)
{
	results_t results;
	nb_partition_t *part = alloca(nb_partition_get_memsize(part_type));
	nb_partition_init(part, part_type);

	int status = simulate(problem_data, part, &results,
			      N_vtx, modify_bcond);
	
	CU_ASSERT(0 == status);

	check_results(part, &results);

	TEMPORAL1(part, &results);
	TEMPORAL2(part, &results);

	nb_partition_finish(part);
	results_finish(&results);
}

static int simulate(const char *problem_data,
		    nb_partition_t *part, results_t *results,
		    uint32_t N_vtx,
		    void (*modify_bcond)(const void*,
					 nb_bcond_t*)/* Can be NULL */)
{
	int status = 1;
	vcn_model_t* model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = alloca(bcond_size);
	nb_bcond_init(bcond, 2);
	nb_material_t* material = nb_material_create();
	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;

	char input[255];
	sprintf(input, problem_data, INPUTS_DIR);
	int read_status =
		read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);

	if (NULL != modify_bcond)
		modify_bcond(part, bcond);

	if (0 != read_status)
		goto CLEANUP;

	get_mesh(model, part, N_vtx);

	uint32_t N_nodes = nb_partition_get_N_nodes(part);
	uint32_t N_elems = nb_partition_get_N_elems(part);
	results_init(results, N_nodes, N_elems);

	int status_cvfa =
		nb_cvfa_compute_2D_Solid_Mechanics(part, material, bcond,
						   false, NULL,
						   analysis2D,
						   &params2D,
						   results->disp,
						   results->strain);

	if (0 != status_cvfa) {
		show_cvfa_error_msg(status_cvfa);
		goto CLEANUP;
	}

	nb_cvfa_compute_stress_from_strain(part,  material, analysis2D,
					   results->strain,
					   results->stress);

	status = 0;
CLEANUP:
	vcn_model_finish(model);
	nb_bcond_finish(bcond);
	nb_material_destroy(material);

	return status;
}

static void get_mesh(const vcn_model_t *model, void *part,
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

	nb_partition_load_from_mesh(part, mesh);
	vcn_mesh_finish(mesh);
}

static void results_init(results_t *results, uint32_t N_nodes,
			 uint32_t N_elems)
{
	uint32_t size_disp = N_elems * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_elems * 3 * sizeof(*(results->strain));
	uint32_t total_size = size_disp + 2 * size_strain;
	char *memblock = malloc(total_size);

	results->N_nodes = N_nodes;
	results->N_elems = N_elems;
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
		 nb_material_t* mat,
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


static int read_material(vcn_cfreader_t *cfreader, nb_material_t *mat)
{
	int status = 1;
	double poisson_module;
	if (0 != vcn_cfreader_read_double(cfreader, &poisson_module))
		goto EXIT;
	nb_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (0 != vcn_cfreader_read_double(cfreader, &elasticity_module))
		goto EXIT;
	nb_material_set_elasticity_module(mat, elasticity_module);

	double fracture_energy;
	if (0 != vcn_cfreader_read_double(cfreader, &fracture_energy))
		goto EXIT;
	nb_material_set_fracture_energy(mat, fracture_energy);

	double compression_limit_stress;
	if (0 != vcn_cfreader_read_double(cfreader, &compression_limit_stress))
		goto EXIT;
	nb_material_set_compression_limit_stress(mat,
						      compression_limit_stress);

	double traction_limit_stress;
	if (0 != vcn_cfreader_read_double(cfreader, &traction_limit_stress))
		goto EXIT;
	nb_material_set_traction_limit_stress(mat, traction_limit_stress);
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
