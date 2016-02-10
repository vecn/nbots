#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "vcn/container_bot.h"
#include "vcn/cfreader_cat.h"
#include "vcn/eigen_bot.h"
#include "vcn/geometric_bot.h"
#include "vcn/pde_bot/finite_element/solid_mechanics/static_elasticity2D.h"

#include "test_library.h"
#include "test_add.h"

#define INPUTS_DIR "../tests/vcn/pde_bot/finite_element/solid_mechanics/static_elasticity2D_UT_inputs"

#include "vcn/geometric_bot/mesh/modules2D/exporter_cairo.h" /* TEMPORAL */
#include "vcn/pde_bot/finite_element/modules/exporter_cairo.h" /* TEMPORAL */

#define POW2(a) ((a)*(a))

static bool check_static_elasticity2D(void);

static double* get_total_disp(uint32_t N, const double *displacement);

static int read_initial_conditions
		(const char* filename,
		 vcn_model_t *model,
		 vcn_bcond_t* bcond,
		 vcn_fem_material_t* mat,
		 char* enable_plane_stress_analysis,
		 double *thickness);
static int read_geometry(vcn_cfreader_t *cfreader, vcn_model_t *model);
static int read_boundary_conditions(vcn_cfreader_t *cfreader,
				     vcn_bcond_t *bcond);
static void read_dirichlet_on_vtx(vcn_cfreader_t *cfreader,
				  vcn_bcond_t *bcond);
static void read_bcond_row(vcn_cfreader_t *cfreader,
			   vcn_bcond_t *bcond, uint32_t id);
static int read_material(vcn_cfreader_t *cfreader, vcn_fem_material_t *mat);
static int read_elasticity2D_params(vcn_cfreader_t *cfreader,
				     char* enable_plane_stress_analysis,
				     double *thickness);

inline int vcn_test_get_driver_id(void)
{
	return VCN_DRIVER_UNIT_TEST;
}

void vcn_test_load_tests(void *tests_ptr)
{
	vcn_test_add(tests_ptr, check_static_elasticity2D,
		     "Check static_elasticity2D()");
}

static bool check_static_elasticity2D(void)
{
	vcn_model_t* model = vcn_model_create();
	vcn_bcond_t* bcond = vcn_fem_bcond_create();
	vcn_fem_material_t* material = vcn_fem_material_create();
	char enable_plane_stress_analysis = 1;
	double thickness;

	char input[255];
	//sprintf(input, "%s/beam_fixed_on_sides.txt", INPUTS_DIR);
	sprintf(input, "%s/beam_cantilever.txt", INPUTS_DIR);
	int read_status =
		read_initial_conditions(input, model, bcond, material,
					&enable_plane_stress_analysis,
					&thickness);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_set_geometric_constraint(mesh,
					  VCN_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  0.4);
	vcn_mesh_generate_from_model(mesh, model);
	vcn_mesh_save_png(mesh, "TEMPORAL.png", 1000, 800);/* TEMPORAL */
	vcn_msh3trg_t* delaunay = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, true, true);
	vcn_mesh_destroy(mesh);

	vcn_bcond_t* bmeshcond =
		vcn_fem_bcond_create_from_model_to_mesh(delaunay, bcond);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(VCN_TRG_LINEAR);

	double* displacement = 
		malloc(delaunay->N_vertices * 2 * sizeof(*displacement));
	double* strain = 
		malloc(delaunay->N_triangles * 3 * sizeof(*strain));

	int status_fem =
		vcn_fem_compute_2D_Solid_Mechanics(delaunay, elemtype, material,
						   bmeshcond,
						   false, NULL,
						   vcn_sparse_solve_Cholesky,
						   enable_plane_stress_analysis,
						   thickness,
						   1, NULL,
						   displacement,
						   strain,
						   "static_elasticity2D_UT.log");
	if (0 != status_fem)
		goto CLEANUP_FEM;
	
	double *total_disp = get_total_disp(delaunay->N_vertices,
					    displacement);
	
	nb_fem_save_png(delaunay, total_disp, "TEMPORAL_FEM.png", 1000, 800);/* TEMPORAL */

	free(total_disp);
CLEANUP_FEM:
	vcn_fem_bcond_destroy(bmeshcond);
	vcn_msh3trg_destroy(delaunay);
	vcn_fem_elem_destroy(elemtype);
	free(displacement);
	free(strain);
CLEANUP_INPUT:
	vcn_model_destroy(model);
	vcn_fem_bcond_destroy(bcond);
	vcn_fem_material_destroy(material);
	return false;
}


static double* get_total_disp(uint32_t N, const double *displacement)
{
	double *total_disp = malloc(N * sizeof(*total_disp));
	for (uint32_t i = 0; i < N; i++)
		total_disp[i] = sqrt(POW2(displacement[i * 2]) +
				     POW2(displacement[i*2+1]));
	return total_disp;	
}

static int read_initial_conditions
		(const char* filename,
		 vcn_model_t *model,
		 vcn_bcond_t* bcond,
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
	if (0 != read_boundary_conditions(cfreader, bcond)) {
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

static int read_boundary_conditions(vcn_cfreader_t *cfreader,
				    vcn_bcond_t *bcond)
{
	bcond->N_dof = 2;
	read_dirichlet_on_vtx(cfreader, bcond);
	
	/* Read Neuman conditions upon vertices */
	if (vcn_cfreader_read_uint(cfreader, &(bcond->N_Neuman_on_vtx)) != 0) {
		vcn_cfreader_destroy(cfreader);
		printf("Error: The 'number of Neuman conditions on \n");
		printf("       vertices' can not be readed.\n");
		return 1;
	}
	if (bcond->N_Neuman_on_vtx > 0) {
		bcond->Neuman_on_vtx_idx = 
			malloc(bcond->N_Neuman_on_vtx * sizeof(uint32_t));
		bcond->Neuman_on_vtx_dof_mask =
			calloc(bcond->N_dof * 
			       bcond->N_Neuman_on_vtx, sizeof(bool));
		bcond->Neuman_on_vtx_val =
			calloc(bcond->N_dof * 
			       bcond->N_Neuman_on_vtx, 
			       sizeof(double));
	}
	for (uint32_t i = 0; i < bcond->N_Neuman_on_vtx; i++) {
		/* Read vertex id of Neuman condition */
		if(vcn_cfreader_read_uint(cfreader,
					      &(bcond->Neuman_on_vtx_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			printf("Error: Can't read 'vertex index' of Neuman conditions. \n");
			return 1;
		}
		/* Read mask of Neuman conditions */
		for (uint32_t j = 0; j < bcond->N_dof; j++){
			int mask;
			if(vcn_cfreader_read_int(cfreader, &mask) != 0){
				vcn_cfreader_destroy(cfreader);
				printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
				return 1;
			}
			bcond->Neuman_on_vtx_dof_mask
				[i * bcond->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Neuman condition components */
		for (uint32_t j = 0; j < bcond->N_dof; j++){
			if (bcond->Neuman_on_vtx_dof_mask
			   [i * bcond->N_dof + j]) {
				if (vcn_cfreader_read_double(cfreader,
							     &(bcond->Neuman_on_vtx_val[i * bcond->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					printf("Error: Can't read 'values' of Neuman conditions. \n");
					return 1;
				}
			}
		}
	}

	/* Read Dirichlet conditions upon segments */
	if (vcn_cfreader_read_uint(cfreader, &(bcond->N_Dirichlet_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		printf("Error: The 'number of Dirichlet conditions on \n");
		printf("       segments' can not be readed.\n");
		return 1;
	}
	if (bcond->N_Dirichlet_on_sgm > 0) {
		bcond->Dirichlet_on_sgm_idx = 
			malloc(bcond->N_Dirichlet_on_sgm * sizeof(uint32_t));
		bcond->Dirichlet_on_sgm_dof_mask =
			calloc(bcond->N_dof * 
			       bcond->N_Dirichlet_on_sgm, sizeof(bool));
		bcond->Dirichlet_on_sgm_val =
			calloc(bcond->N_dof * 
			       bcond->N_Dirichlet_on_sgm, 
			       sizeof(double));
	}
	for (uint32_t i = 0; i < bcond->N_Dirichlet_on_sgm; i++) {
		/* Read vertex id of Dirichlet condition */
		if(vcn_cfreader_read_uint(cfreader,
					      &(bcond->Dirichlet_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			printf("Error: Can't read 'segment index' of Dirichlet conditions. \n");
			return 1;
		}
		/* Read mask of Dirichlet conditions */
		for (uint32_t j = 0; j < bcond->N_dof; j++){
			int mask;
			if (vcn_cfreader_read_int(cfreader, &mask) != 0) {
				vcn_cfreader_destroy(cfreader);
				printf("Error: Can't read 'DoF mask' of Dirichlet conditions. \n");
				return 1;
			}
			bcond->Dirichlet_on_sgm_dof_mask
				[i * bcond->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Dirichlet condition components */
		for (uint32_t j = 0; j < bcond->N_dof; j++) {
			if (bcond->Dirichlet_on_sgm_dof_mask
			    [i * bcond->N_dof + j]) {
				if(vcn_cfreader_read_double(cfreader,
							&(bcond->Dirichlet_on_sgm_val[i * bcond->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					printf("Error: Can't read 'values' of Dirichlet conditions. \n");
					return 1;
				}
			}
		}
	}

	/* Read Neuman conditions upon segments */
	if (vcn_cfreader_read_uint(cfreader, &(bcond->N_Neuman_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		printf("Error: The 'number of Neuman conditions on \n");
		printf("       segments' can not be readed.\n");
		return 1;
	}
	if (bcond->N_Neuman_on_sgm > 0) {
		bcond->Neuman_on_sgm_idx = 
			malloc(bcond->N_Neuman_on_sgm * sizeof(uint32_t));
		bcond->Neuman_on_sgm_dof_mask =
			calloc(bcond->N_dof * 
			       bcond->N_Neuman_on_sgm, sizeof(bool));
		bcond->Neuman_on_sgm_val =
			calloc(bcond->N_dof * 
			       bcond->N_Neuman_on_sgm, 
			       sizeof(double));
	}
	for (uint32_t i = 0; i < bcond->N_Neuman_on_sgm; i++) {
		/* Read vertex id of Neuman condition */
		if (vcn_cfreader_read_uint(cfreader,
				       &(bcond->Neuman_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			printf("Error: Can't read 'segment index' of Neuman conditions. \n");
			return 1;
		}
		/* Read mask of Neuman conditions */
		for (uint32_t j = 0; j < bcond->N_dof; j++) {
			int mask;
			if (vcn_cfreader_read_int(cfreader, &mask) != 0) {
				vcn_cfreader_destroy(cfreader);
				printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
				return 1;
			}
			bcond->Neuman_on_sgm_dof_mask
				[i * bcond->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Neuman condition components */
		for (uint32_t j = 0; j < bcond->N_dof; j++) {
			if (bcond->Neuman_on_sgm_dof_mask
			    [i * bcond->N_dof + j]) {
				if(vcn_cfreader_read_double(cfreader,
							&(bcond->Neuman_on_sgm_val[i * bcond->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					printf("Error: Can't read 'values' of Neuman conditions. \n");
					return 1;
				}
			}
		}
	}
	return 0;
}

static void read_dirichlet_on_vtx(vcn_cfreader_t *cfreader,
				  vcn_bcond_t *bcond)
{
	if (0 != vcn_cfreader_read_uint(cfreader,
					&(bcond->N_Dirichlet_on_vtx)))
		goto EXIT;

	if (0 < bcond->N_Dirichlet_on_vtx) {
		bcond->Dirichlet_on_vtx_idx = 
			malloc(bcond->N_Dirichlet_on_vtx *
			       sizeof(*(bcond->Dirichlet_on_vtx_idx)));
		bcond->Dirichlet_on_vtx_dof_mask =
			calloc(bcond->N_dof * bcond->N_Dirichlet_on_vtx,
			       sizeof(*(bcond->Dirichlet_on_vtx_dof_mask)));
		bcond->Dirichlet_on_vtx_val =
			calloc(bcond->N_dof * bcond->N_Dirichlet_on_vtx, 
			       sizeof(*(bcond->Dirichlet_on_vtx_val)));
		for (uint32_t i = 0; i < bcond->N_Dirichlet_on_vtx; i++)
			read_bcond_row(cfreader, bcond, i);
	}
EXIT:
	return;
}

static void read_bcond_row(vcn_cfreader_t *cfreader,
			   vcn_bcond_t *bcond, uint32_t id)
{
	/* Read vertex id  */
	uint32_t elem_id;
	if(0 != vcn_cfreader_read_uint(cfreader, &elem_id))
		goto EXIT;
	bcond->Dirichlet_on_vtx_idx[id] = elem_id;
	/* Read mask */
	for (uint32_t j = 0; j < bcond->N_dof; j++) {
		int mask;
		if (0 != vcn_cfreader_read_int(cfreader, &mask))
			goto EXIT;
				
		bcond->Dirichlet_on_vtx_dof_mask[id * bcond->N_dof + j] =
			(mask == 1);
	}
	/* Read components */
	for (uint32_t j = 0; j < bcond->N_dof; j++) {
		if (bcond->Dirichlet_on_vtx_dof_mask[id * bcond->N_dof + j]) {
			double val;
			if (0 != vcn_cfreader_read_double(cfreader, &val))
				goto EXIT;
			bcond->Dirichlet_on_vtx_val[id * bcond->N_dof + j] = val;
		}
	}
EXIT:
	return;
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
