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

static bool check_static_elasticity2D(void);

vcn_model_t* read_initial_conditions(
	const char* filename,
	vcn_bcond_t* bconditions,  /* Output */
	vcn_fem_material_t* mat,                     /* Output */
	char* enable_plane_stress_analysis,  /* Output */
	double *thickness);                  /* Output */

void output_save_dma(const char* filename,
		     const char* author,
		     const char* project_name,
		     bool enable_double_precision,
		     uint32_t N_vertices, 
		     double* vertices,
		     uint32_t N_elements,
		     uint32_t* elements_connectivity_matrix,
		     double* displacement,
		     double* strain);

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
	/* Read optimization parameters */
	vcn_bcond_t* bconditions = vcn_fem_bcond_create();
	vcn_fem_material_t* material = vcn_fem_material_create();
	char enable_plane_stress_analysis;
	double thickness;

	char input[255];
	sprintf(input, "%s/beam_fixed_on_sides.in", INPUTS_DIR);
	vcn_model_t* model = 
		read_initial_conditions(input, bconditions, 
					material,
					&enable_plane_stress_analysis,
					&thickness);

	if (NULL == model) {
		vcn_fem_bcond_destroy(bconditions);
		vcn_fem_material_destroy(material);
		printf("Error: Reading Input file.\n");
		return 1;
	}

	/* Mesh domain */
	vcn_mesh_t* mesh = vcn_mesh_create();
	vcn_mesh_generate_from_model(mesh, model);
	vcn_msh3trg_t* delaunay = 
		vcn_mesh_get_msh3trg(mesh, true, true, true, false, NULL);
	vcn_mesh_destroy(mesh);

	/* Write logfile */
	printf("Mesh nodes: %i\n", delaunay->N_vertices);
	printf("Mesh elements: %i\n", delaunay->N_triangles);

	vcn_bcond_t* bmeshcond =
		vcn_fem_bcond_create_from_model_to_mesh(delaunay, bconditions);

	/* FEM Analysis */
	vcn_fem_elem_t* elemtype = vcn_fem_elem_create(VCN_TRG_LINEAR);

	double* displacement = 
		malloc(delaunay->N_vertices * 2 * sizeof(*displacement));
	double* strain = 
		malloc(delaunay->N_triangles * 3 * sizeof(*strain));

	vcn_fem_compute_2D_Solid_Mechanics(delaunay, elemtype, material,
					   bmeshcond,
					   false, NULL,
					   vcn_sparse_solve_Cholesky,
					   enable_plane_stress_analysis,
					   thickness,
					   2, NULL,
					   displacement,
					   strain,
					   "static_elasticity2D_UT.log");

	output_save_dma("static_elasticity2D_UT.out",
			"Victor Eduardo Cardoso Nungaray",
			"Benchmark of SMFEM",
			false,
			delaunay->N_vertices,
			delaunay->vertices,
			delaunay->N_triangles,
			delaunay->vertices_forming_triangles,
			displacement,
			strain);

	/* Free memory */
	vcn_fem_bcond_destroy(bconditions);
	vcn_fem_bcond_destroy(bmeshcond);
	vcn_model_destroy(model);
	vcn_msh3trg_destroy(delaunay);
	vcn_fem_material_destroy(material);
	vcn_fem_elem_destroy(elemtype);
	free(displacement);
	free(strain);

	/* Successful exit */
	return 0;
}

vcn_model_t* read_initial_conditions(
                   const char* filename,
		   vcn_bcond_t* bconditions, /* Output */
		   vcn_fem_material_t* mat,                    /* Output */
		   char* enable_plane_stress_analysis, /* Output */
		   double *thickness)                  /* Output */
{      
	uint32_t i, j; /* Iterative variables */

	/* Initialize custom format to read file */
	vcn_cfreader_t* cfreader = vcn_cfreader_create(filename, "#");
	if (NULL == cfreader)
		return NULL;
	/* Read modele vertices */
	uint32_t N_model_vtx = 0;
	if(vcn_cfreader_read_uint32_t(cfreader, &N_model_vtx) != 0){
		vcn_cfreader_destroy(cfreader);
		printf("Error: The 'number of vertices' in the modeletry can't be readed.\n");
		return NULL;
	}
	if (N_model_vtx < 1) {
		vcn_cfreader_destroy(cfreader);
		printf("Error: There are not vertices to read in the input file.\n");
		return NULL;
	}
	double *model_vtx = malloc(2 * N_model_vtx * sizeof(double));
	for (i=0; i < 2 * N_model_vtx; i++) {
		if(vcn_cfreader_read_double(cfreader, &(model_vtx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			printf("Error: The 'vertices' of the modeletry can't be readed.\n");
			return NULL;
		}
	}
	/* Read model segments */
	uint32_t N_model_sgm = 0;
	if (vcn_cfreader_read_uint32_t(cfreader, &N_model_sgm) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		printf("Error: The 'number of segments' in the modeletry can't be readed.\n");
		return NULL;
	}
	if (N_model_sgm < 1) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		printf("Error: There are not segments to read in the input file\n");
		return NULL;    
	}
	uint32_t *model_sgm = malloc(2 * N_model_sgm * sizeof(*model_sgm));
	for (i = 0; i < 2 * N_model_sgm; i++) {
		if (vcn_cfreader_read_uint32_t(cfreader, &(model_sgm[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			free(model_sgm);
			printf("Error: The 'segments' of the modeletry can't be readed.\n");
			return NULL;
		}
	}
	/* Read model holes */
	uint32_t N_model_holes;
	if(vcn_cfreader_read_uint32_t(cfreader, &N_model_holes) != 0){
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		printf("Error: The 'number of holes' in the modeletry can't be readed.\n");
		return NULL;
	}
	double *model_holes = NULL;
	if (N_model_holes > 0) {
		model_holes = (double*)malloc(2 * N_model_holes * sizeof(double));
		for (i = 0; i < 2 * N_model_holes; i++) {
			if (vcn_cfreader_read_double(cfreader, &(model_holes[i])) != 0) {
				vcn_cfreader_destroy(cfreader);
				free(model_vtx);
				free(model_sgm);
				free(model_holes);
				printf("Error: The 'holes' of the modeletry can't be readed.\n");
				return NULL;
			}
		}
	}
	/* Read boundary conditions */
	bconditions->N_dof = 2; /* Solid Mechanics just solve displacements */
	/* Read Dirichlet conditions upon vertices */
	if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Dirichlet_on_vtx)) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'number of Dirichlet conditions on \n");
		printf("       vertices' can not be readed.\n");
		return NULL;
	}
	if (bconditions->N_Dirichlet_on_vtx > 0) {
		bconditions->Dirichlet_on_vtx_idx = 
			malloc(bconditions->N_Dirichlet_on_vtx * sizeof(uint32_t));
		bconditions->Dirichlet_on_vtx_dof_mask =
			calloc(bconditions->N_dof * 
			       bconditions->N_Dirichlet_on_vtx, sizeof(bool));
		bconditions->Dirichlet_on_vtx_val =
			calloc(bconditions->N_dof * 
			       bconditions->N_Dirichlet_on_vtx, 
			       sizeof(double));
	}
	for (i = 0; i < bconditions->N_Dirichlet_on_vtx; i++) {
		/* Read vertex id of Dirichlet condition */
		if(vcn_cfreader_read_uint32_t(cfreader,
					      &(bconditions->Dirichlet_on_vtx_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			free(model_sgm);
			if (N_model_holes > 0)
				free(model_holes);
			printf("Error: Can't read 'vertex index' of Dirichlet conditions. \n");
			return NULL;
    }
		/* Read mask of Dirichlet conditions */
		for (j=0; j < bconditions->N_dof; j++) {
			int mask;
			if (vcn_cfreader_read_int(cfreader, &mask) != 0){
				vcn_cfreader_destroy(cfreader);
				free(model_vtx);
				free(model_sgm);
				if(N_model_holes > 0)
					free(model_holes);
				printf("Error: Can't read 'DoF mask' of Dirichlet conditions. \n");
				return NULL;
			}
			bconditions->Dirichlet_on_vtx_dof_mask
				[i * bconditions->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Dirichlet condition components */
		for (j = 0; j < bconditions->N_dof; j++) {
			if (bconditions->Dirichlet_on_vtx_dof_mask
			    [i * bconditions->N_dof + j]) {
				if (vcn_cfreader_read_double(cfreader,
							     &(bconditions->Dirichlet_on_vtx_val[i * bconditions->N_dof + j])) != 0){
					vcn_cfreader_destroy(cfreader);
					free(model_vtx);
					free(model_sgm);
					if (N_model_holes > 0)
						free(model_holes);
					printf("Error: Can't read 'values' of Dirichlet conditions. \n");
					return NULL;
				}
			}
		}
	}
	
	/* Read Neuman conditions upon vertices */
	if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Neuman_on_vtx)) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'number of Neuman conditions on \n");
		printf("       vertices' can not be readed.\n");
		return NULL;
	}
	if (bconditions->N_Neuman_on_vtx > 0) {
		bconditions->Neuman_on_vtx_idx = 
			malloc(bconditions->N_Neuman_on_vtx * sizeof(uint32_t));
		bconditions->Neuman_on_vtx_dof_mask =
			calloc(bconditions->N_dof * 
			       bconditions->N_Neuman_on_vtx, sizeof(bool));
		bconditions->Neuman_on_vtx_val =
			calloc(bconditions->N_dof * 
			       bconditions->N_Neuman_on_vtx, 
			       sizeof(double));
	}
	for (i = 0; i < bconditions->N_Neuman_on_vtx; i++) {
		/* Read vertex id of Neuman condition */
		if(vcn_cfreader_read_uint32_t(cfreader,
					      &(bconditions->Neuman_on_vtx_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			free(model_sgm);
			if (N_model_holes > 0)
				free(model_holes);
			printf("Error: Can't read 'vertex index' of Neuman conditions. \n");
			return NULL;
		}
		/* Read mask of Neuman conditions */
		for (j = 0; j < bconditions->N_dof; j++){
			int mask;
			if(vcn_cfreader_read_int(cfreader, &mask) != 0){
				vcn_cfreader_destroy(cfreader);
				free(model_vtx);
				free(model_sgm);
				if (N_model_holes > 0)
					free(model_holes);
				printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
				return NULL;
			}
			bconditions->Neuman_on_vtx_dof_mask
				[i * bconditions->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Neuman condition components */
		for (j = 0; j < bconditions->N_dof; j++){
			if (bconditions->Neuman_on_vtx_dof_mask
			   [i * bconditions->N_dof + j]) {
				if (vcn_cfreader_read_double(cfreader,
							     &(bconditions->Neuman_on_vtx_val[i * bconditions->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					free(model_vtx);
					free(model_sgm);
					if (N_model_holes > 0)
						free(model_holes);
					printf("Error: Can't read 'values' of Neuman conditions. \n");
					return NULL;
				}
			}
		}
	}

	/* Read Dirichlet conditions upon segments */
	if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Dirichlet_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'number of Dirichlet conditions on \n");
		printf("       segments' can not be readed.\n");
		return NULL;
	}
	if (bconditions->N_Dirichlet_on_sgm > 0) {
		bconditions->Dirichlet_on_sgm_idx = 
			malloc(bconditions->N_Dirichlet_on_sgm * sizeof(uint32_t));
		bconditions->Dirichlet_on_sgm_dof_mask =
			calloc(bconditions->N_dof * 
			       bconditions->N_Dirichlet_on_sgm, sizeof(bool));
		bconditions->Dirichlet_on_sgm_val =
			calloc(bconditions->N_dof * 
			       bconditions->N_Dirichlet_on_sgm, 
			       sizeof(double));
	}
	for (i = 0; i < bconditions->N_Dirichlet_on_sgm; i++) {
		/* Read vertex id of Dirichlet condition */
		if(vcn_cfreader_read_uint32_t(cfreader,
					      &(bconditions->Dirichlet_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			free(model_sgm);
			if (N_model_holes > 0)
				free(model_holes);
			printf("Error: Can't read 'segment index' of Dirichlet conditions. \n");
			return NULL;
		}
		/* Read mask of Dirichlet conditions */
		for (j=0; j < bconditions->N_dof; j++){
			int mask;
			if (vcn_cfreader_read_int(cfreader, &mask) != 0) {
				vcn_cfreader_destroy(cfreader);
				free(model_vtx);
				free(model_sgm);
				if (N_model_holes > 0)
					free(model_holes);
				printf("Error: Can't read 'DoF mask' of Dirichlet conditions. \n");
				return NULL;
			}
			bconditions->Dirichlet_on_sgm_dof_mask
				[i * bconditions->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Dirichlet condition components */
		for (j = 0; j < bconditions->N_dof; j++) {
			if (bconditions->Dirichlet_on_sgm_dof_mask
			    [i * bconditions->N_dof + j]) {
				if(vcn_cfreader_read_double(cfreader,
							&(bconditions->Dirichlet_on_sgm_val[i * bconditions->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					free(model_vtx);
					free(model_sgm);
					if (N_model_holes > 0)
						free(model_holes);
					printf("Error: Can't read 'values' of Dirichlet conditions. \n");
					return NULL;
				}
			}
		}
	}

	/* Read Neuman conditions upon segments */
	if (vcn_cfreader_read_uint32_t(cfreader, &(bconditions->N_Neuman_on_sgm)) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'number of Neuman conditions on \n");
		printf("       segments' can not be readed.\n");
		return NULL;
	}
	if (bconditions->N_Neuman_on_sgm > 0) {
		bconditions->Neuman_on_sgm_idx = 
			malloc(bconditions->N_Neuman_on_sgm * sizeof(uint32_t));
		bconditions->Neuman_on_sgm_dof_mask =
			calloc(bconditions->N_dof * 
			       bconditions->N_Neuman_on_sgm, sizeof(bool));
		bconditions->Neuman_on_sgm_val =
			calloc(bconditions->N_dof * 
			       bconditions->N_Neuman_on_sgm, 
			       sizeof(double));
	}
	for (i = 0; i < bconditions->N_Neuman_on_sgm; i++) {
		/* Read vertex id of Neuman condition */
		if (vcn_cfreader_read_uint32_t(cfreader,
				       &(bconditions->Neuman_on_sgm_idx[i])) != 0) {
			vcn_cfreader_destroy(cfreader);
			free(model_vtx);
			free(model_sgm);
			if (N_model_holes > 0)
				free(model_holes);
			printf("Error: Can't read 'segment index' of Neuman conditions. \n");
			return NULL;
		}
		/* Read mask of Neuman conditions */
		for (j=0; j < bconditions->N_dof; j++) {
			int mask;
			if (vcn_cfreader_read_int(cfreader, &mask) != 0) {
				vcn_cfreader_destroy(cfreader);
				free(model_vtx);
				free(model_sgm);
				if (N_model_holes > 0)
					free(model_holes);
				printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
				return NULL;
			}
			bconditions->Neuman_on_sgm_dof_mask
				[i * bconditions->N_dof + j] = (mask==1)?true:false;
		}
		/* Read Neuman condition components */
		for (j = 0; j < bconditions->N_dof; j++) {
			if (bconditions->Neuman_on_sgm_dof_mask
			    [i * bconditions->N_dof + j]) {
				if(vcn_cfreader_read_double(cfreader,
							&(bconditions->Neuman_on_sgm_val[i * bconditions->N_dof + j])) != 0) {
					vcn_cfreader_destroy(cfreader);
					free(model_vtx);
					free(model_sgm);
					if (N_model_holes > 0)
						free(model_holes);
					printf("Error: Can't read 'values' of Neuman conditions. \n");
					return NULL;
				}
			}
		}
	}
  
	/* Read materials properties */
	double poisson_module;
	if (vcn_cfreader_read_double(cfreader, &poisson_module) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'poisson module' of the material can't be readed.\n");
		return NULL;
	}
	vcn_fem_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (vcn_cfreader_read_double(cfreader, &elasticity_module) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'elasticity module' of the material can't be readed.\n");
		return NULL;
	}
	vcn_fem_material_set_elasticity_module(mat, elasticity_module);

	double fracture_energy;
	if (vcn_cfreader_read_double(cfreader, &fracture_energy) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'fracture energy' of the material can't be readed.\n");
		return NULL;
	}
	vcn_fem_material_set_fracture_energy(mat, fracture_energy);

	double compression_limit_stress;
	if (vcn_cfreader_read_double(cfreader, &compression_limit_stress) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'max strength' of the material can't be readed.\n");
		return NULL;
	}
	vcn_fem_material_set_compression_limit_stress(mat, compression_limit_stress);

	double traction_limit_stress;
	if (vcn_cfreader_read_double(cfreader, &traction_limit_stress) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'max strength' of the material can't be readed.\n");
		return NULL;
	}
	vcn_fem_material_set_traction_limit_stress(mat, traction_limit_stress);

	/* Read analysis params and output flags */
	int iaux;
	if (vcn_cfreader_read_int(cfreader, &iaux) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'plane stress enabling flag' \
can't be readed.\n");
		return NULL;
	}
	*enable_plane_stress_analysis = iaux;

	if (vcn_cfreader_read_double(cfreader, thickness) != 0) {
		vcn_cfreader_destroy(cfreader);
		free(model_vtx);
		free(model_sgm);
		if (N_model_holes > 0)
			free(model_holes);
		printf("Error: The 'thickness' of the material can't be readed.\n");
		return NULL;
	}

	/* Create domain's modeletry structure */
	vcn_model_t* model = vcn_model_load_from_arrays(
		model_vtx, N_model_vtx,
		model_sgm, N_model_sgm,
		model_holes, N_model_holes);
	/* Free memory */
	vcn_cfreader_destroy(cfreader);
	free(model_vtx);
	free(model_sgm);
	if (N_model_holes > 0)
		free(model_holes);

	/* Return domain's modeletry structure */
	return model;
}

void output_save_dma(const char* filename,
		     const char* author,
		     const char* project_name,
		     bool enable_double_precision,
		     uint32_t N_vertices, 
		     double* vertices,
		     uint32_t N_elements,
		     uint32_t* elements_connectivity_matrix,
		     double* displacement,
		     double* strain)
{
	int32_t aux; /* Variable to write using 4-byte integers */
	/* Open  file */
	FILE* fp = fopen(filename, "wb");
	/****************************************************************/
	/********************** Write DMA Header ************************/
	/****************************************************************/
	/* Write signature */
	const char signature[21] = "DynamicMeshAttributes";
	fwrite(signature, 1, 21, fp);
	/* Write version */
	uint8_t version = 1;
	fwrite(&version, 1, 1, fp);
	/* Write author's name */
	uint8_t author_name_length = strlen(author);
	fwrite(&author_name_length, 1, 1, fp);
	fwrite(author, 1, author_name_length, fp);
	/* Write mesh name */
	uint8_t project_name_length = strlen(project_name);
	fwrite(&project_name_length, 1, 1, fp);
	fwrite(project_name, 1, project_name_length, fp);
	/* Write configuration number */
	uint8_t configuration = 1;
	if(enable_double_precision)
		configuration = 5;
	/* Configuration: 
	 *   > 2D
	 *   > [Precision by parameter]
	 *   > Only dynamic attributes (Fixed vertices)
	 *   > Variable time step (We don't use the time step)
	 */
	fwrite(&configuration, 1, 1, fp);
	/* Write number of vertices */
	aux = N_vertices;
	fwrite(&aux, 4, 1, fp);

	/****************************************************************/
	/************ Write General information of attributes ***********/
	/****************************************************************/
	uint8_t N_attributes = 2;
	fwrite(&N_attributes, 1, 1, fp);

	/* Write displacement data */
	{
		uint8_t attribute_name_length = 12;
		fwrite(&attribute_name_length, 1, 1, fp);
		char attribute_name[12] = "Displacement";
		fwrite(attribute_name, 1, 12, fp);
		char attribute_type = 'p'; /* Dynamic vector of vertices */
		fwrite(&attribute_type, 1, 1, fp);
		uint8_t vec_length = 2;
		fwrite(&vec_length, 1, 1, fp);
	}

	/* Write strain data */
	{
		uint8_t attribute_name_length = 6;
		fwrite(&attribute_name_length, 1, 1, fp);
		char attribute_name[6] = "Strain";
		fwrite(attribute_name, 1, 6, fp);
		char attribute_type = 'P'; /* Dynamic vector on vertices */
		fwrite(&attribute_type, 1, 1, fp);
		uint8_t vec_length = 3;
		fwrite(&vec_length, 1, 1, fp);
	}

	/****************************************************************/
	/******** Write vertices coordinates and displacements **********/
	/****************************************************************/
	if (enable_double_precision) {
		fwrite(vertices, 8, N_vertices * 2, fp);
		fwrite(displacement, 8, N_vertices * 2, fp);
	} else {
		for (uint32_t i=0; i < N_vertices * 2; i++) {
			float v = vertices[i]; /* Cast to simple precision */
			fwrite(&v, 4, 1, fp);
		}
		for (uint32_t i=0; i < N_vertices * 2; i++) {
			float v = displacement[i]; /* Cast to simple precision */
			fwrite(&v, 4, 1, fp);
		}
	}

  
	/****************************************************************/
	/****************** Write mesh connectivity *********************/
	/****************************************************************/
	aux = N_elements; /* Number of triangular elements */
	fwrite(&aux, 4, 1, fp);
	for (uint32_t i=0; i < N_elements; i++) {
		aux = i;
		fwrite(&aux, 4, 1, fp); /* Write element ID */
		/* Write vertices ids conforming the element */
		for (int j = 0; j < 3; j++) {
			aux = elements_connectivity_matrix[i*3+j];
			fwrite(&aux, 4, 1, fp);
		}
	}
	aux = 0;          /* Number of cuadrilateral elements */
	fwrite(&aux, 4, 1, fp);

	/* Write strains */
	if (enable_double_precision) {
		fwrite(strain, 8, N_elements * 3, fp);
	} else {
		for (uint32_t i=0; i < N_elements * 3; i++) {
			float v = strain[i]; /* Cast to simple precision */
			fwrite(&v, 4, 1, fp);
		}
	}
	/******************* End Writing DMA format *********************/
	fclose(fp);
}
