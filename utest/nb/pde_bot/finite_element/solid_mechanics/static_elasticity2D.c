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
#include "nb/pde_bot/boundary_conditions/bcond.h"
#include "nb/pde_bot/boundary_conditions/bcond_read.h"
#include "nb/pde_bot/material.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#define INPUTS_DIR "../../../../utest/nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D_inputs"
#define OUTPUT "../../../"

#include "nb/geometric_bot/mesh/modules2D/exporter_cairo.h" /* TEMPORAL */
#include "nb/pde_bot/finite_element/modules/exporter_cairo.h" /* TEMPORAL */

#define POW2(a) ((a)*(a))

typedef struct {
	double *disp;
	double *strain;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_static_elasticity2D(void);

static void results_init(results_t *results,
			 const vcn_msh3trg_t *msh3trg);
static void results_finish(results_t *results);
static double* get_total_disp(uint32_t N, const double *displacement);
static int read_initial_conditions
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
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/static_elasticity.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "static_elasticity2D()",
		    test_static_elasticity2D);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_static_elasticity2D(void)
{
	vcn_model_t* model = alloca(vcn_model_get_memsize());
	vcn_model_init(model);
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = alloca(bcond_size);
	nb_bcond_init(bcond, 2);
	vcn_fem_material_t* material = vcn_fem_material_create();
	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;

	char input[255];
	sprintf(input, "%s/beam_cantilever.txt", INPUTS_DIR);
	int read_status =
		read_initial_conditions(input, model, bcond, material,
					&analysis2D, &params2D);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     100);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	vcn_mesh_generate_from_model(mesh, model);

	vcn_msh3trg_t* delaunay = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);
	vcn_mesh_destroy(mesh);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(NB_TRG_LINEAR);

	results_t results;
	results_init(&results, delaunay);

	int status_fem =
		vcn_fem_compute_2D_Solid_Mechanics(delaunay, elemtype,
						   material, bcond,
						   false, NULL,
						   analysis2D,
						   &params2D, NULL,
						   results.disp, results.strain);
	if (0 != status_fem)
		goto CLEANUP_FEM;
	
	double *total_disp = get_total_disp(delaunay->N_vertices,
					    results.disp);
	
	char filename[100];
	sprintf(filename, "%s/TEMPORAL_FEM.png", OUTPUT);
	nb_fem_save_png(delaunay, total_disp, filename, 1000, 800);/* TEMPORAL */

	free(total_disp);
CLEANUP_FEM:
	vcn_msh3trg_destroy(delaunay);
	vcn_fem_elem_destroy(elemtype);
	results_finish(&results);
CLEANUP_INPUT:
	vcn_model_finish(model);
	nb_bcond_finish(bcond);
	vcn_fem_material_destroy(material);
	CU_ASSERT(false);
}


static void results_init(results_t *results,
			 const vcn_msh3trg_t *msh3trg)
{
	uint32_t size_disp = msh3trg->N_vertices * 2 *
		sizeof(*(results->disp));
	uint32_t size_strain = msh3trg->N_triangles * 3 *
		sizeof(*(results->strain));
	uint32_t total_size = size_disp + size_strain;
	char *memblock = malloc(total_size);

	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
}

static inline void results_finish(results_t *results)
{
	free(results->disp);
}

static double* get_total_disp(uint32_t N, const double *displacement)
{
	double *total_disp = malloc(N * sizeof(*total_disp));
	for (uint32_t i = 0; i < N; i++) {
		total_disp[i] = sqrt(POW2(displacement[i * 2]) +
				     POW2(displacement[i*2+1]));
	}
	
	return total_disp;	
}

static int read_initial_conditions
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

	uint32_t *model_sgm = malloc(2 * model->M * sizeof(*model->edge));
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
