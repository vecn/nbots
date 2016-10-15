#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
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
static void check_beam_cantilever(const void *part,
				  const results_t *results);
static void test_plate_with_hole(void);
static void check_plate_with_hole(const void *part,
				  const results_t *results);
static double get_error_avg_pwh(const void *part,
				const nb_fem_elem_t* elem,
				const double *stress);
static void get_cartesian_gpoint(uint32_t id_elem, int8_t id_gp,
				 const void *part,
				 const nb_fem_elem_t* elem,
				 double gp[2]);
static void get_analytic_stress_pwh(double x, double y, double stress[3]);
static void modify_bcond_pwh(const void *part,
			     nb_bcond_t *bcond);
static void pwh_DC_SGM_cond(double *x, double t, double *out);
static void pwh_BC_SGM_cond(double *x, double t, double *out);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */);
static int simulate(const char *problem_data,
		    nb_partition_t *part, results_t *results,
		    uint32_t N_vtx,
		    void (*modify_bcond)(const void*,
					     nb_bcond_t*)/* Can be NULL */);
static void get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx);

static void results_init(results_t *results, uint32_t N_vtx, uint32_t N_trg);
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

static void check_beam_cantilever(const void *part,
				  const results_t *results)
{
	double max_disp = 0;
	uint32_t N_nodes = nb_partition_get_N_nodes(part);
	for (uint32_t i = 0; i < N_nodes; i++) {
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
	run_test("%s/plate_with_hole.txt", 1000,
		 check_plate_with_hole,
		 modify_bcond_pwh);
}

static void check_plate_with_hole(const void *part,
				  const results_t *results)
{
	nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);

	double avg_error = get_error_avg_pwh(part, elem, results->stress);

	CU_ASSERT(avg_error < 9.7e-3);
	
	nb_fem_elem_destroy(elem);
}

static double get_error_avg_pwh(const void *part,
				const nb_fem_elem_t* elem,
				const double *stress)
{
	double avg = 0.0;
	uint32_t N = 0;
	uint32_t N_elems = nb_partition_get_N_elems(part);
	for (uint32_t i = 0; i < N_elems; i++) {
		int8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
		N += N_gp;
		for (int8_t p = 0; p < N_gp; p++) {
			double gp[2];
			get_cartesian_gpoint(i, p, part, elem, gp);
			double analytic_stress[3];
			get_analytic_stress_pwh(gp[0], gp[1], analytic_stress);
			uint32_t id = i * N_gp + p;
			double vm_stress =
			  nb_pde_get_vm_stress(stress[id * 3],
					       stress[id*3+1],
					       stress[id*3+2]);
			double analytic_vm_stress =
			  nb_pde_get_vm_stress(analytic_stress[0],
					       analytic_stress[1],
					       analytic_stress[2]);

			double error = fabs(1.0 - vm_stress /
					    analytic_vm_stress);
			avg += error;
		}
	}
	return avg /= N;
}

static void get_cartesian_gpoint(uint32_t id_elem, int8_t id_gp,
				 const void *part,
				 const nb_fem_elem_t* elem,
				 double gp[2])
{
	int8_t N = nb_fem_elem_get_N_nodes(elem);
	gp[0] = 0.0;
	gp[1] = 0.0;
	for (int i = 0; i < N; i++) {
		double Ni = nb_fem_elem_Ni(elem, i, id_gp);
		uint32_t vi = nb_partition_elem_get_adj(part, id_elem, i);
		gp[0] += Ni * nb_partition_node_get_x(part, vi);
		gp[1] += Ni * nb_partition_node_get_y(part, vi);
	}
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
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */)
{
	results_t results;
	nb_partition_t *part = nb_allocate_on_stack(nb_partition_get_memsize(NB_TRIAN));
	nb_partition_init(part, NB_TRIAN);

	int status = simulate(problem_data, part, &results,
			      N_vtx, modify_bcond);
	
	CU_ASSERT(0 == status);

	check_results(part, &results);

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

	if (NULL != modify_bcond)
		modify_bcond(part, bcond);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	get_mesh(model, part, N_vtx);

	nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);

	uint32_t N_nodes = nb_partition_get_N_nodes(part);
	uint32_t N_elems = nb_partition_get_N_elems(part);
	results_init(results, N_nodes, N_elems);

	int status_fem =
		nb_fem_compute_2D_Solid_Mechanics(part, elem,
						   material, bcond,
						   false, NULL,
						   analysis2D,
						   &params2D, NULL,
						   results->disp,
						   results->strain);
	if (0 != status_fem)
		goto CLEANUP_FEM;

	nb_fem_compute_stress_from_strain(N_elems, elem, material,
					   analysis2D, results->strain, NULL,
					   results->stress);

	status = 0;
CLEANUP_FEM:
	nb_fem_elem_destroy(elem);
CLEANUP_INPUT:
	nb_model_finish(model);
	nb_bcond_finish(bcond);
	nb_material_destroy(material);

	return status;
}

static void get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx)
{
	uint32_t mesh_memsize = nb_tessellator2D__get_memsize();
	nb_tessellator2D__t* mesh = nb_allocate_on_stack(mesh_memsize);
	nb_tessellator2D__init(mesh);
	nb_tessellator2D__set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     N_vtx);
	nb_tessellator2D__set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	nb_tessellator2D__generate_from_model(mesh, model);

	nb_partition_load_from_mesh(part, mesh);
	nb_tessellator2D__finish(mesh);
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
	nb_cfreader_destroy(cfreader);
	status = 0;
EXIT:
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

	model->vertex = malloc(2 * model->N * sizeof(*(model->vertex)));
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

	model->edge = malloc(2 * model->M * sizeof(*model->edge));
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
		model->holes = malloc(2 * model->H * sizeof(*(model->holes)));
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
