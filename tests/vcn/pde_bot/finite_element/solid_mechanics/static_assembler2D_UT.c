#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "libre_utils/libre_utils.h"
#include "libre_cfreader/libre_cfreader.h"
#include "libre_dstructs/libre_dstructs.h"
#include "libre_solvers/libre_solvers.h"
#include "libre_grid/libre_grid.h"
#include "libre_mesh/libre_mesh.h"
#include "libre_smfem/libre_smfem.h"

model_t* read_initial_conditions(
		   const char* filename,
		   boundary_conditions_t* bconditions,  /* Output */
		   material_t* mat,                     /* Output */
		   char* enable_plane_stress_analysis,  /* Output */
		   double *thickness);                  /* Output */

void output_save_dma(const char* filename,
		     const char* author,
		     const char* project_name,
		     bool enable_double_precision,
		     uint N_vertices, 
		     double* vertices,
		     uint N_elements,
		     uint* elements_connectivity_matrix,
		     double* displacement,
		     double* strain);

int main(int argc, char* argv[]){
  /* Validate input */
  if(argc < 3){
    printf("Usage options:\n > program.bin input_parameters.in output.dma\n");
    return 1;
  }

  /* Read optimization parameters */
  boundary_conditions_t* bconditions = bconditions_create();
  material_t* material = material_create();
  char enable_plane_stress_analysis;
  double thickness;
  model_t* model = 
    read_initial_conditions(argv[1], bconditions, 
			    material,
			    &enable_plane_stress_analysis,
			    &thickness);

  if(model == NULL){
    bconditions_destroy(bconditions);
    material_destroy(material);
    printf("Error: Reading Input file.\n");
    return 1;
  }

  /* Mesh domain */
  size_control_params_t sc_params;
  sc_params.min_size = 0.1;
  sc_params.max_size = 0.1;
  sc_params.max_sgm_size = 0.0;
  mesh_t* mesh = 
    mesh_create_from_model(model, 2000, NULL, &sc_params);
  mesh_delaunay_t* delaunay = 
    mesh_get_delaunay(mesh, true, true, true, false, NULL);
  mesh_destroy(mesh);

  /* Write logfile */
  printf("Mesh nodes: %i\n", delaunay->N_vertices);
  printf("Mesh elements: %i\n", delaunay->N_triangles);
    
  boundary_conditions_t* bmeshcond =
     bconditions_create_from_model_to_mesh(delaunay, bconditions);

  /* FEM Analysis */
  element_type_t* elemtype = elemtype_create_triangle();

  double* displacement = 
    (double*) malloc(delaunay->N_vertices * 2 * sizeof(double));
  double* strain = 
    (double*) malloc(delaunay->N_triangles * 3 * sizeof(double));

  FEM_compute_2D_Solid_Mechanics
    (delaunay->N_vertices,
     delaunay->vertices,
     delaunay->N_triangles,
     delaunay->vertices_forming_triangles,
     elemtype,
     material,
     bmeshcond,
     false, NULL, /* Disable considering self weight */
     sparse_solve_Cholesky,
     enable_plane_stress_analysis,
     thickness,
     2, NULL,
     displacement,
     strain,
     "benchmark.smfem.logfile");

  /* Save in dma format */
  output_save_dma(argv[2],
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
  bconditions_destroy(bconditions);
  bconditions_destroy(bmeshcond);
  model_destroy(model);
  mesh_delaunay_destroy(delaunay);
  material_destroy(material);
  elemtype_destroy(elemtype);
  free(displacement);
  free(strain);

  /* Successful exit */
  return 0;
}

model_t* read_initial_conditions(
                   const char* filename,
		   boundary_conditions_t* bconditions, /* Output */
		   material_t* mat,                    /* Output */
		   char* enable_plane_stress_analysis, /* Output */
		   double *thickness)                  /* Output */
{      
  uint i, j; /* Iterative variables */

  /* Initialize custom format to read file */
  customized_format_reader_t* cfreader = cfreader_create(filename, "#");
  if(cfreader == NULL) return NULL;
  /* Read modele vertices */
  uint N_model_vtx = 0;
  if(cfreader_read_uint(cfreader, &N_model_vtx) != 0){
    cfreader_destroy(cfreader);
    printf("Error: The 'number of vertices' in the modeletry can't be readed.\n");
    return NULL;
  }
  if(N_model_vtx < 1){
    cfreader_destroy(cfreader);
    printf("Error: There are not vertices to read in the input file.\n");
    return NULL;
  }
  double *model_vtx = (double*)malloc(2 * N_model_vtx * sizeof(double));
  for(i=0; i < 2 * N_model_vtx; i++){
    if(cfreader_read_double(cfreader, &(model_vtx[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      printf("Error: The 'vertices' of the modeletry can't be readed.\n");
      return NULL;
    }
  }
  /* Read model segments */
  uint N_model_sgm = 0;
  if(cfreader_read_uint(cfreader, &N_model_sgm) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    printf("Error: The 'number of segments' in the modeletry can't be readed.\n");
    return NULL;
  }
  if(N_model_sgm < 1){
    cfreader_destroy(cfreader);
    free(model_vtx);
    printf("Error: There are not segments to read in the input file\n");
    return NULL;    
  }
  uint *model_sgm = (uint*)malloc(2 * N_model_sgm * sizeof(uint));
  for(i=0; i < 2 * N_model_sgm; i++){
    if(cfreader_read_uint(cfreader, &(model_sgm[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      free(model_sgm);
      printf("Error: The 'segments' of the modeletry can't be readed.\n");
      return NULL;
    }
  }
  /* Read modele holes */
  uint N_model_holes;
  if(cfreader_read_uint(cfreader, &N_model_holes) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    printf("Error: The 'number of holes' in the modeletry can't be readed.\n");
    return NULL;
  }
  double *model_holes = NULL;
  if(N_model_holes > 0){
    model_holes = (double*)malloc(2 * N_model_holes * sizeof(double));
    for(i=0; i < 2 * N_model_holes; i++){
      if(cfreader_read_double(cfreader, &(model_holes[i])) != 0){
	cfreader_destroy(cfreader);
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
  if(cfreader_read_uint(cfreader, &(bconditions->N_Dirichlet_on_vtx)) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'number of Dirichlet conditions on \n");
    printf("       vertices' can not be readed.\n");
    return NULL;
  }
  if(bconditions->N_Dirichlet_on_vtx > 0){
    bconditions->Dirichlet_on_vtx_idx = 
      (uint*)malloc(bconditions->N_Dirichlet_on_vtx * sizeof(uint));
    bconditions->Dirichlet_on_vtx_dof_mask = (bool*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Dirichlet_on_vtx, sizeof(bool));
    bconditions->Dirichlet_on_vtx_val = (double*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Dirichlet_on_vtx, 
	     sizeof(double));
  }
  for(i = 0; i < bconditions->N_Dirichlet_on_vtx; i++){
    /* Read vertex id of Dirichlet condition */
    if(cfreader_read_uint(cfreader,
			  &(bconditions->Dirichlet_on_vtx_idx[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      free(model_sgm);
      if(N_model_holes > 0)
	free(model_holes);
      printf("Error: Can't read 'vertex index' of Dirichlet conditions. \n");
      return NULL;
    }
    /* Read mask of Dirichlet conditions */
    for(j=0; j < bconditions->N_dof; j++){
      int mask;
      if(cfreader_read_int(cfreader, &mask) != 0){
	cfreader_destroy(cfreader);
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
    for(j = 0; j < bconditions->N_dof; j++){
      if(bconditions->Dirichlet_on_vtx_dof_mask
	 [i * bconditions->N_dof + j]){
	if(cfreader_read_double(cfreader,
				&(bconditions->Dirichlet_on_vtx_val
				  [i * bconditions->N_dof + j])) != 0){
	  cfreader_destroy(cfreader);
	  free(model_vtx);
	  free(model_sgm);
	  if(N_model_holes > 0)
	    free(model_holes);
	  printf("Error: Can't read 'values' of Dirichlet conditions. \n");
	  return NULL;
	}
      }
    }
  }

  /* Read Neuman conditions upon vertices */
  if(cfreader_read_uint(cfreader, &(bconditions->N_Neuman_on_vtx)) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'number of Neuman conditions on \n");
    printf("       vertices' can not be readed.\n");
    return NULL;
  }
  if(bconditions->N_Neuman_on_vtx > 0){
    bconditions->Neuman_on_vtx_idx = 
      (uint*)malloc(bconditions->N_Neuman_on_vtx * sizeof(uint));
    bconditions->Neuman_on_vtx_dof_mask = (bool*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Neuman_on_vtx, sizeof(bool));
    bconditions->Neuman_on_vtx_val = (double*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Neuman_on_vtx, 
	     sizeof(double));
  }
  for(i = 0; i < bconditions->N_Neuman_on_vtx; i++){
    /* Read vertex id of Neuman condition */
    if(cfreader_read_uint(cfreader,
			  &(bconditions->Neuman_on_vtx_idx[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      free(model_sgm);
      if(N_model_holes > 0)
	free(model_holes);
      printf("Error: Can't read 'vertex index' of Neuman conditions. \n");
      return NULL;
    }
    /* Read mask of Neuman conditions */
    for(j=0; j < bconditions->N_dof; j++){
      int mask;
      if(cfreader_read_int(cfreader, &mask) != 0){
	cfreader_destroy(cfreader);
	free(model_vtx);
	free(model_sgm);
	if(N_model_holes > 0)
	  free(model_holes);
	printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
	return NULL;
      }
      bconditions->Neuman_on_vtx_dof_mask
	[i * bconditions->N_dof + j] = (mask==1)?true:false;
    }
    /* Read Neuman condition components */
    for(j = 0; j < bconditions->N_dof; j++){
      if(bconditions->Neuman_on_vtx_dof_mask
	 [i * bconditions->N_dof + j]){
	if(cfreader_read_double(cfreader,
				&(bconditions->Neuman_on_vtx_val
				  [i * bconditions->N_dof + j])) != 0){
	  cfreader_destroy(cfreader);
	  free(model_vtx);
	  free(model_sgm);
	  if(N_model_holes > 0)
	    free(model_holes);
	  printf("Error: Can't read 'values' of Neuman conditions. \n");
	  return NULL;
	}
      }
    }
  }

  /* Read Dirichlet conditions upon segments */
  if(cfreader_read_uint(cfreader, &(bconditions->N_Dirichlet_on_sgm)) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'number of Dirichlet conditions on \n");
    printf("       segments' can not be readed.\n");
    return NULL;
  }
  if(bconditions->N_Dirichlet_on_sgm > 0){
    bconditions->Dirichlet_on_sgm_idx = 
      (uint*)malloc(bconditions->N_Dirichlet_on_sgm * sizeof(uint));
    bconditions->Dirichlet_on_sgm_dof_mask = (bool*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Dirichlet_on_sgm, sizeof(bool));
    bconditions->Dirichlet_on_sgm_val = (double*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Dirichlet_on_sgm, 
	     sizeof(double));
  }
  for(i = 0; i < bconditions->N_Dirichlet_on_sgm; i++){
    /* Read vertex id of Dirichlet condition */
    if(cfreader_read_uint(cfreader,
			  &(bconditions->Dirichlet_on_sgm_idx[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      free(model_sgm);
      if(N_model_holes > 0)
	free(model_holes);
      printf("Error: Can't read 'segment index' of Dirichlet conditions. \n");
      return NULL;
    }
    /* Read mask of Dirichlet conditions */
    for(j=0; j < bconditions->N_dof; j++){
      int mask;
      if(cfreader_read_int(cfreader, &mask) != 0){
	cfreader_destroy(cfreader);
	free(model_vtx);
	free(model_sgm);
	if(N_model_holes > 0)
	  free(model_holes);
	printf("Error: Can't read 'DoF mask' of Dirichlet conditions. \n");
	return NULL;
      }
      bconditions->Dirichlet_on_sgm_dof_mask
	[i * bconditions->N_dof + j] = (mask==1)?true:false;
    }
    /* Read Dirichlet condition components */
    for(j = 0; j < bconditions->N_dof; j++){
      if(bconditions->Dirichlet_on_sgm_dof_mask
	 [i * bconditions->N_dof + j]){
	if(cfreader_read_double(cfreader,
				&(bconditions->Dirichlet_on_sgm_val
				  [i * bconditions->N_dof + j])) != 0){
	  cfreader_destroy(cfreader);
	  free(model_vtx);
	  free(model_sgm);
	  if(N_model_holes > 0)
	    free(model_holes);
	  printf("Error: Can't read 'values' of Dirichlet conditions. \n");
	  return NULL;
	}
      }
    }
  }
  
  /* Read Neuman conditions upon segments */
  if(cfreader_read_uint(cfreader, &(bconditions->N_Neuman_on_sgm)) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'number of Neuman conditions on \n");
    printf("       segments' can not be readed.\n");
    return NULL;
  }
  if(bconditions->N_Neuman_on_sgm > 0){
    bconditions->Neuman_on_sgm_idx = 
      (uint*)malloc(bconditions->N_Neuman_on_sgm * sizeof(uint));
    bconditions->Neuman_on_sgm_dof_mask = (bool*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Neuman_on_sgm, sizeof(bool));
    bconditions->Neuman_on_sgm_val = (double*)
      calloc(bconditions->N_dof * 
	     bconditions->N_Neuman_on_sgm, 
	     sizeof(double));
  }
  for(i = 0; i < bconditions->N_Neuman_on_sgm; i++){
    /* Read vertex id of Neuman condition */
    if(cfreader_read_uint(cfreader,
			  &(bconditions->Neuman_on_sgm_idx[i])) != 0){
      cfreader_destroy(cfreader);
      free(model_vtx);
      free(model_sgm);
      if(N_model_holes > 0)
	free(model_holes);
      printf("Error: Can't read 'segment index' of Neuman conditions. \n");
      return NULL;
    }
    /* Read mask of Neuman conditions */
    for(j=0; j < bconditions->N_dof; j++){
      int mask;
      if(cfreader_read_int(cfreader, &mask) != 0){
	cfreader_destroy(cfreader);
	free(model_vtx);
	free(model_sgm);
	if(N_model_holes > 0)
	  free(model_holes);
	printf("Error: Can't read 'DoF mask' of Neuman conditions. \n");
	return NULL;
      }
      bconditions->Neuman_on_sgm_dof_mask
	[i * bconditions->N_dof + j] = (mask==1)?true:false;
    }
    /* Read Neuman condition components */
    for(j = 0; j < bconditions->N_dof; j++){
      if(bconditions->Neuman_on_sgm_dof_mask
	 [i * bconditions->N_dof + j]){
	if(cfreader_read_double(cfreader,
				&(bconditions->Neuman_on_sgm_val
				  [i * bconditions->N_dof + j])) != 0){
	  cfreader_destroy(cfreader);
	  free(model_vtx);
	  free(model_sgm);
	  if(N_model_holes > 0)
	    free(model_holes);
	  printf("Error: Can't read 'values' of Neuman conditions. \n");
	  return NULL;
	}
      }
    }
  }
  
  /* Read materials properties */
  double poisson_module;
  if(cfreader_read_double(cfreader, &poisson_module) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'poisson module' of the material can't be readed.\n");
    return NULL;
  }
  material_set_poisson_module(mat, poisson_module);

  double elasticity_module;
  if(cfreader_read_double(cfreader, &elasticity_module) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'elasticity module' of the material can't be readed.\n");
    return NULL;
  }
  material_set_elasticity_module(mat, elasticity_module);

  double fracture_energy;
  if(cfreader_read_double(cfreader, &fracture_energy) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'fracture energy' of the material can't be readed.\n");
    return NULL;
  }
  material_set_fracture_energy(mat, fracture_energy);

  double compression_limit_stress;
  if(cfreader_read_double(cfreader, &compression_limit_stress) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'max strength' of the material can't be readed.\n");
    return NULL;
  }
  material_set_compression_limit_stress(mat, compression_limit_stress);

  double traction_limit_stress;
  if(cfreader_read_double(cfreader, &traction_limit_stress) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'max strength' of the material can't be readed.\n");
    return NULL;
  }
  material_set_traction_limit_stress(mat, traction_limit_stress);

  /* Read analysis params and output flags */
  int iaux;
  if(cfreader_read_int(cfreader, &iaux) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'plane stress enabling flag' \
            can't be readed.\n");
    return NULL;
  }
  *enable_plane_stress_analysis = iaux;

  if(cfreader_read_double(cfreader, thickness) != 0){
    cfreader_destroy(cfreader);
    free(model_vtx);
    free(model_sgm);
    if(N_model_holes > 0)
      free(model_holes);
    printf("Error: The 'thickness' of the material can't be readed.\n");
    return NULL;
  }

  /* Create domain's modeletry structure */
  model_t* model = model_create(model_vtx, N_model_vtx,
					      model_sgm, N_model_sgm,
					      model_holes, N_model_holes);
  /* Free memory */
  cfreader_destroy(cfreader);
  free(model_vtx);
  free(model_sgm);
  if(N_model_holes > 0)
    free(model_holes);

  /* Return domain's modeletry structure */
  return model;
}

void output_save_dma(const char* filename,
		     const char* author,
		     const char* project_name,
		     bool enable_double_precision,
		     uint N_vertices, 
		     double* vertices,
		     uint N_elements,
		     uint* elements_connectivity_matrix,
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
  uchar version = 1;
  fwrite(&version, 1, 1, fp);
  /* Write author's name */
  uchar author_name_length = strlen(author);
  fwrite(&author_name_length, 1, 1, fp);
  fwrite(author, 1, author_name_length, fp);
  /* Write mesh name */
  uchar project_name_length = strlen(project_name);
  fwrite(&project_name_length, 1, 1, fp);
  fwrite(project_name, 1, project_name_length, fp);
  /* Write configuration number */
  uchar configuration = 1;
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
  uchar N_attributes = 2;
  fwrite(&N_attributes, 1, 1, fp);
  
  /* Write displacement data */
  {
    uchar attribute_name_length = 12;
    fwrite(&attribute_name_length, 1, 1, fp);
    char attribute_name[12] = "Displacement";
    fwrite(attribute_name, 1, 12, fp);
    char attribute_type = 'p'; /* Dynamic vector of vertices */
    fwrite(&attribute_type, 1, 1, fp);
    uchar vec_length = 2;
    fwrite(&vec_length, 1, 1, fp);
  }

  /* Write strain data */
  {
    uchar attribute_name_length = 6;
    fwrite(&attribute_name_length, 1, 1, fp);
    char attribute_name[6] = "Strain";
    fwrite(attribute_name, 1, 6, fp);
    char attribute_type = 'P'; /* Dynamic vector on vertices */
    fwrite(&attribute_type, 1, 1, fp);
    uchar vec_length = 3;
    fwrite(&vec_length, 1, 1, fp);
  }

  /****************************************************************/
  /******** Write vertices coordinates and displacements **********/
  /****************************************************************/
  if(enable_double_precision){
    fwrite(vertices, 8, N_vertices * 2, fp);
    fwrite(displacement, 8, N_vertices * 2, fp);
  }else{
    for(register uint i=0; i < N_vertices * 2; i++){
      float v = vertices[i]; /* Cast to simple precision */
      fwrite(&v, 4, 1, fp);
    }
    for(register uint i=0; i < N_vertices * 2; i++){
      float v = displacement[i]; /* Cast to simple precision */
      fwrite(&v, 4, 1, fp);
    }
  }

  
  /****************************************************************/
  /****************** Write mesh connectivity *********************/
  /****************************************************************/
  aux = N_elements; /* Number of triangular elements */
  fwrite(&aux, 4, 1, fp);
  for(register uint i=0; i < N_elements; i++){
    aux = i;
    fwrite(&aux, 4, 1, fp); /* Write element ID */
    /* Write vertices ids conforming the element */
    for(register int j = 0; j < 3; j++){
      aux = elements_connectivity_matrix[i*3+j];
      fwrite(&aux, 4, 1, fp);
    }
  }
  aux = 0;          /* Number of cuadrilateral elements */
  fwrite(&aux, 4, 1, fp);

  /* Write strains */
  if(enable_double_precision){
    fwrite(strain, 8, N_elements * 3, fp);
  }else{
    for(register uint i=0; i < N_elements * 3; i++){
      float v = strain[i]; /* Cast to simple precision */
      fwrite(&v, 4, 1, fp);
    }
  }
  /******************* End Writing DMA format *********************/
  fclose(fp);
}
