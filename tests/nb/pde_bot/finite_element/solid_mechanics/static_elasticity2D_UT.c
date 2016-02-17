#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <alloca.h>

#include "nb/container_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "test_library.h"
#include "test_add.h"

#define INPUTS_DIR "../tests/nb/pde_bot/finite_element/solid_mechanics/static_elasticity2D_UT_inputs"

#include "nb/geometric_bot/mesh/modules2D/exporter_cairo.h" /* TEMPORAL */
#include "nb/pde_bot/finite_element/modules/exporter_cairo.h" /* TEMPORAL */

#define POW2(a) ((a)*(a))

static bool check_static_elasticity2D(void);

static double* get_total_disp(uint32_t N, const double *displacement);

static int read_initial_conditions
		(const char* filename,
		 vcn_model_t *model,
		 nb_bcond_t* bcond,
		 vcn_fem_material_t* mat,
		 char* enable_plane_stress_analysis,
		 double *thickness);
static int read_geometry(vcn_cfreader_t *cfreader, vcn_model_t *model);
static int read_material(vcn_cfreader_t *cfreader, vcn_fem_material_t *mat);
static int read_elasticity2D_params(vcn_cfreader_t *cfreader,
				     char* enable_plane_stress_analysis,
				     double *thickness);

inline int vcn_test_get_driver_id(void)
{
	return NB_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_static_elasticity2D,
		     "Check static_elasticity2D()");
}

static bool check_static_elasticity2D(void)
{
	vcn_model_t* model = vcn_model_create();
	uint16_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = alloca(bcond_size);
	nb_bcond_init(bcond, 2);
	vcn_fem_material_t* material = vcn_fem_material_create();
	char enable_plane_stress_analysis = 1;
	double thickness;

	char input[255];
	sprintf(input, "%s/beam_cantilever.txt", INPUTS_DIR);
	int read_status =
		read_initial_conditions(input, model, bcond, material,
					&enable_plane_stress_analysis,
					&thickness);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     200);
	vcn_mesh_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_mesh_save_png(mesh, "TEMPORAL.png", 1000, 800);/* TEMPORAL */
	vcn_msh3trg_t* delaunay = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);
	vcn_mesh_destroy(mesh);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(NB_TRG_LINEAR);

	double* displacement = 
		malloc(delaunay->N_vertices * 2 * sizeof(*displacement));
	double* strain = 
		malloc(delaunay->N_triangles * 3 * sizeof(*strain));

	int status_fem =
		vcn_fem_compute_2D_Solid_Mechanics(delaunay, elemtype,
						   material, bcond,
						   false, NULL,
						   enable_plane_stress_analysis,
						   thickness, NULL,
						   displacement, strain);
	if (0 != status_fem)
		goto CLEANUP_FEM;
	
	double *total_disp = get_total_disp(delaunay->N_vertices,
					    displacement);
	
	nb_fem_save_png(delaunay, total_disp, "TEMPORAL_FEM.png", 1000, 800);/* TEMPORAL */

	free(total_disp);
CLEANUP_FEM:
	vcn_msh3trg_destroy(delaunay);
	vcn_fem_elem_destroy(elemtype);
	free(displacement);
	free(strain);
CLEANUP_INPUT:
	vcn_model_destroy(model);
	nb_bcond_finish(bcond);
	vcn_fem_material_destroy(material);
	return false;
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
		 char* enable_plane_stress_analysis,
		 double *thickness)
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
	if (0 != read_elasticity2D_params(cfreader, 
					  enable_plane_stress_analysis,
					  thickness)) {
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
	uint32_t N_model_vtx = 0;
	if (0 != vcn_cfreader_read_uint(cfreader, &N_model_vtx))
		goto EXIT;

	if (1 > N_model_vtx)
		goto EXIT;

	double *model_vtx = malloc(2 * N_model_vtx * sizeof(double));
	for (uint32_t i = 0; i < 2 * N_model_vtx; i++) {
		if (0 != vcn_cfreader_read_double(cfreader, &(model_vtx[i])))
			goto CLEANUP_VERTICES;
	}
	/* Read model segments */
	uint32_t N_model_sgm = 0;
	if (0 != vcn_cfreader_read_uint(cfreader, &N_model_sgm))
		goto CLEANUP_VERTICES;

	if (1 > N_model_sgm)
		goto CLEANUP_VERTICES;

	uint32_t *model_sgm = malloc(2 * N_model_sgm * sizeof(*model_sgm));
	for (uint32_t i = 0; i < 2 * N_model_sgm; i++) {
		if (0 != vcn_cfreader_read_uint(cfreader, &(model_sgm[i])))
			goto CLEANUP_SEGMENTS;
	}
	/* Read model holes */
	uint32_t N_model_holes;
	if (0 != vcn_cfreader_read_uint(cfreader, &N_model_holes))
		goto CLEANUP_SEGMENTS;

	double *model_holes = NULL;
	if (0 < N_model_holes) {
		model_holes = malloc(2 * N_model_holes * sizeof(double));
		for (uint32_t i = 0; i < 2 * N_model_holes; i++) {
			if (0 != vcn_cfreader_read_double(cfreader,
							  &(model_holes[i])))
				goto CLEANUP_HOLES;
		}
	}
	vcn_model_load_from_arrays(model,
				   model_vtx, N_model_vtx,
				   model_sgm, N_model_sgm,
				   model_holes, N_model_holes);
	status = 0;
CLEANUP_HOLES:
	if (0 < N_model_holes)
		free(model_holes);
CLEANUP_SEGMENTS:
	free(model_sgm);
CLEANUP_VERTICES:
	free(model_vtx);
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
				    char* enable_plane_stress_analysis,
				    double *thickness)
{
	int status = 1;
	int iaux;
	if (0 != vcn_cfreader_read_int(cfreader, &iaux))
		goto EXIT;
	*enable_plane_stress_analysis = (char)iaux;

	if (0 != vcn_cfreader_read_double(cfreader, thickness))
		goto EXIT;
	status = 0;
EXIT:
	return status;
}
