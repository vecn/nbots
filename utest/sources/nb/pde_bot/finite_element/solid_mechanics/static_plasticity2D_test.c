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
/*#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "nb/pde_bot/finite_element/solid_mechanics/pipeline.h"
#define INPUTS_DIR "/home/david/Desktop/nbots/utest_plastic/static_plasticity2D_inputs"
//#define INPUTS_DIR "C:/Users/David/Desktop/nbots/utest_plastic/static_plasticity2D_inputs"*/

#define INPUTS_DIR "../utest/sources/nb/pde_bot/static_plasticity2D_inputs"

#define POW2(a) ((a)*(a))

typedef struct {
	double *disp;
	double *strain;
	double *stress;
	bool *plastic_elements;
	double *damage;
    	uint32_t N_vtx;
	uint32_t N_trg;
	double *nodal_strain;
	double *nodal_stress;
	double *nodal_damage;
} plastic_results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_beam_cantilever(void);
static void test_continuous_beam(void);
static void check_beam_cantilever(const void *part,
				  const plastic_results_t *results);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     	void (*check_results)(const void*,
			const plastic_results_t*),
			char *printable_name);
static int simulate(const char *problem_data,
		    nb_mesh2D_t *part, plastic_results_t *results,
		    uint32_t N_vtx, char *printable_name);
static void get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx);
static void print_plastic_results(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                    double *total_strain, const nb_mesh2D_t *const part,
                    bool *plastic_elements, const char *problem_data, char *printable_name);
static void results_init_plastic(plastic_results_t *results, uint32_t N_vtx, uint32_t N_trg);
static void results_finish(plastic_results_t *results);
static int read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(nb_cfreader_t *cfreader, nb_model_t *model);
static int read_plastic_material(nb_cfreader_t *cfreader, nb_material_t *mat);
static int read_plasticity2D_params(nb_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static void plastic_check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D);


void cunit_nb_pde_bot_fem_sm_static_plasticity(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/"\
			     "static_plasticity2D.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Beam cantilever", test_beam_cantilever);
	CU_add_test(suite, "Continuous beam", test_continuous_beam);
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
        run_test("%s/plastic_beam_cantilever.txt", 500,
        check_beam_cantilever, "cantilever");
}
static void check_beam_cantilever(const void *part,
				  const plastic_results_t *results)
{
    /*
	double max_X_stress = 0;
	double max_Y_stress = 0;
	double max_XY_stress = 0;
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	for (uint32_t i = 0; i < N_elem; i++) {
        double X_stress = results->stress[3*i];
        double Y_stress = results->stress[3*i+1];
        double XY_stress = results->stress[3*i+2];
		if (abs(max_X_stress) < abs(X_stress))
			max_X_stress = X_stress;
		if (abs(max_Y_stress) < abs(Y_stress))
			max_Y_stress = Y_stress;
		if (abs(max_XY_stress) < abs(XY_stress))
			max_XY_stress = XY_stress;
	}
    if (abs(max_X_stress - 1.0836e9) < 1e-2) {
        printf("The maximum Sx is inside tolerance");
    }
    else {
        printf("Error in Sx! It's value is: %lf \n", max_X_stress);
    }
	if (abs(max_Y_stress - 7.2167e8) < 1e-2) {
        printf("The maximum Sy is inside tolerance");
    }
    else {
        printf("Error in Sy! It's value is: %lf \n", max_Y_stress);
    }
    if (abs(max_XY_stress - 3.6373e7) < 1e-2) {
        printf("The maximum Sxy is inside tolerance");
    }
    else {
        printf("Error in Sxy! It's value is: %lf \n", max_XY_stress);
    }
    */
    ;
}

static void test_continuous_beam(void)
{
        run_test("%s/plastic_continuous_beam.txt", 100,
        check_beam_cantilever, "continuous");
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		    	void (*check_results)(const void*,
			const plastic_results_t*), char *printable_name)
{
    	plastic_results_t results;

	nb_mesh2D_t *part = nb_allocate_mem(nb_mesh2D_get_memsize(NB_TRIAN));
	nb_mesh2D_init(part, NB_TRIAN);

	int status = simulate(problem_data, part, &results, N_vtx, printable_name);

    	printf("Simulation status: %i\n", status); /* TEMPORAL */
	if (status == 0) {
        printf("The simulation ran properly. \n");
	}
	//check_results(part, &results);

	nb_mesh2D_finish(part);
	nb_free_mem(part);
    	results_finish(&results);
}

static int simulate(const char *problem_data,
		    nb_mesh2D_t *part, plastic_results_t *results,
		    uint32_t N_vtx, char *printable_name)
{
	int status = 1;
	nb_model_t* model = nb_allocate_mem(nb_model_get_memsize());
	nb_model_init(model);
	uint32_t bcond_size = nb_bcond_get_memsize(2);
	nb_bcond_t *bcond = nb_allocate_mem(bcond_size);
	nb_bcond_init(bcond, 2);
	nb_material_t* material = nb_material_create();
	nb_analysis2D_t analysis2D;
	nb_analysis2D_params params2D;

	char input[255];
	sprintf(input, problem_data, INPUTS_DIR);
	int read_status =
		read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);
    	printf("Read status: %i\n", read_status); /* TEMPORAL */
    	plastic_check_input_values(material, analysis2D, &params2D);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	get_mesh(model, part, N_vtx);

	nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(part);
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	results_init_plastic(results, N_nodes, N_elems);
	uint32_t N_force_steps = 100;
	double accepted_tol = 1.0;
	double gravity[2] = {0.0 , -9.81};
    	int status_fem = fem_compute_plastic_2D_Solid_Mechanics(part, elem, material,
                                                             bcond, true, gravity,
                                                             analysis2D, &params2D, NULL, results->strain, results->stress,
                                                             results->disp, N_force_steps,
                                                             accepted_tol, results->plastic_elements,results->nodal_strain,
						 	     results->nodal_stress);

    	print_plastic_results(N_nodes, N_elems, results->disp, results->stress, results->strain,
                          part, results->plastic_elements, problem_data, printable_name);

	if (0 != status_fem)
		goto CLEANUP_FEM;
	status = 0;

	CLEANUP_FEM:
		nb_fem_elem_destroy(elem);
	CLEANUP_INPUT:
		nb_model_finish(model);
		nb_free_mem(model);
		nb_bcond_finish(bcond);
		nb_free_mem(bcond);
		nb_material_destroy(material);

	return status;
}

static void get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx)
{
	uint32_t mesh_memsize = nb_tessellator2D_get_memsize();
	nb_tessellator2D_t* mesh = nb_allocate_mem(mesh_memsize);
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     N_vtx);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	nb_tessellator2D_generate_from_model(mesh, model);

	nb_mesh2D_load_from_tessellator2D(part, mesh);
	nb_tessellator2D_finish(mesh);
	nb_free_mem(mesh);
}

static void results_init_plastic(plastic_results_t *results, uint32_t N_vtx, uint32_t N_trg)
{
	uint32_t size_disp = N_vtx * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_trg * 3 * sizeof(*(results->strain));
	uint32_t size_plastic_elements = N_trg * sizeof(*(results->plastic_elements));
	uint32_t size_nodal_strain = N_vtx * 3 * sizeof(*(results->nodal_strain));
	uint32_t size_nodal_stress = N_vtx * 3 * sizeof(*(results->nodal_stress));
	uint32_t total_size = size_disp + 2 * size_strain + size_plastic_elements + size_nodal_strain + size_nodal_stress;// + size_damage;
	char *memblock = nb_allocate_mem(total_size);

	results->N_vtx = N_vtx;
	results->N_trg = N_trg;
	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
	results->plastic_elements = (void*)(memblock + size_disp + 2 * size_strain);
	results->nodal_strain = (void*)(memblock + size_disp + 2 * size_strain + size_plastic_elements);
	results->nodal_stress = (void*)(memblock + size_disp + 2 * size_strain + size_plastic_elements + size_nodal_strain);
	results->damage = NULL;//void*)(memblock + size_disp + 2 * size_strain + size_plastic_elements);
	results->nodal_damage = NULL;
}	

static inline void results_finish(plastic_results_t *results)
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
	nb_cfreader_t* cfreader = nb_cfreader_create();
	nb_cfreader_add_line_comment_token(cfreader, "#");
	status = nb_cfreader_open_file(cfreader, filename);
	if (0 != status) {
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
	if (0 != read_plastic_material(cfreader, mat)) {
		printf("\nERROR: Material contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_plasticity2D_params(cfreader, analysis2D, params2D)) {
		printf("\nERROR: Reading problem type in %s.\n",
		       filename);
		goto EXIT;
	}
	status = 0;

	EXIT:
		nb_cfreader_close_file(cfreader);
    		nb_cfreader_destroy(cfreader);

	return status;
}

static int read_geometry(nb_cfreader_t *cfreader, nb_model_t *model)
{
	int status = 1;
	/* Read model vertices */
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

	if (1 > N) {
		goto CLEANUP_VERTICES;
	}
	model->M = N;
	model->edge = nb_allocate_mem(2 * model->M * sizeof(*model->edge));
	for (uint32_t i = 0; i < 2 * model->M; i++) {
		if (0 != nb_cfreader_read_uint(cfreader, &(model->edge[i])))
			goto CLEANUP_SEGMENTS;
	}

	status = 0;
	goto EXIT;

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

static int read_plastic_material(nb_cfreader_t *cfreader, nb_material_t *mat)
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

	double plasticity_module;
	if (0 != nb_cfreader_read_double(cfreader, &plasticity_module))
		goto EXIT;
	nb_material_set_plasticity_module(mat, plasticity_module);

	double yield_stress;
	if (0 != nb_cfreader_read_double(cfreader, &yield_stress))
		goto EXIT;
	nb_material_set_yield_stress(mat, yield_stress);
	double density;
	if (0 != nb_cfreader_read_double(cfreader, &density));
	nb_material_set_density(mat, density);

	status = 0;

	EXIT:
	return status;
}


static int read_plasticity2D_params(nb_cfreader_t *cfreader,
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

static void plastic_check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D) {
	printf("Elastic Modulus: %lf\n", nb_material_get_elasticity_module(material)); /* TEMPORAL */
    	printf("Plastic Modulus: %lf\n", nb_material_get_plasticity_module(material)); /* TEMPORAL */
    	printf("Poisson module: %f\n", nb_material_get_poisson_module(material)); /* TEMPORAL */
    	printf("Density: %lf\n", nb_material_get_density(material)); /* TEMPORAL */
   	 printf("Thickness: %lf\n", params2D->thickness);
    	switch (analysis2D) {
    		case 0:
    			printf("Problem type: NB_PLANE_STRESS\n");
   		break;
    		case 1:
    			printf("Problem type: NB_PLANE_STRAIN\n");
    		break;
    		case 2:
    			printf("Problem type: NB_SOLID_OF_REVOLUTION\n");
    		break;
    		default:
    			printf("Problem type: NB_PLANE_STRESS\n");
    	}
}

static void print_plastic_results(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                            double *total_strain, const nb_mesh2D_t *const part,
                            bool *plastic_elements, const char *problem_data, char *printable_name)
{
    	uint32_t nodesize = N_nod * sizeof(double);
    	uint32_t elemsize = N_elem * sizeof(double);

    	uint32_t avg_displacement_size = nodesize;
    	uint32_t total_displacement_x_size = nodesize;
    	uint32_t total_displacement_y_size = nodesize;
    	uint32_t Sx_size = elemsize;
    	uint32_t Sy_size = elemsize;
    	uint32_t Sxy_size = elemsize;
    	uint32_t Ex_size = elemsize;
   	uint32_t Ey_size = elemsize;
    	uint32_t Exy_size = elemsize;

   	uint64_t memsize = avg_displacement_size + total_displacement_x_size + total_displacement_y_size +
                        Sx_size + Sy_size + Sxy_size + Ex_size + Ey_size + Exy_size;

    	char* memblock = malloc(memsize);

    	double *avg_displacement = (void*)memblock;
    	for(int i = 1; i < N_nod; i++){
        	avg_displacement[i] = sqrt(POW2(total_displacement[2*i])+ POW2(total_displacement[2*i +1]));
    	}

    	double *total_displacement_x = (void*)(memblock + avg_displacement_size);
    	double *total_displacement_y = (void*)(memblock + avg_displacement_size + total_displacement_x_size);

    	for(int j = 0; j < N_nod; j++){
        	total_displacement_x[j] = total_displacement[2*j];
        	total_displacement_y[j] = total_displacement[2*j +1];
    	}
    	double *Sx = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size);
    	double *Sy =(void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size);
    	double *Sxy =(void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size);
   	for(int j = 0; j < N_elem; j++) {
        	Sx[j] = stress[3*j];
        	Sy[j] = stress[3*j + 1];
        	Sxy[j] = stress[3*j +2];
    	}

    	double *Ex = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size);
    	double *Ey = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size + Ex_size);
    	double *Exy = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size + Ex_size + Ey_size);

    	for(int j = 0; j < N_elem; j++) {
        	Ex[j] = total_strain[3*j];
        	Ey[j] = total_strain[3*j + 1];
        	Exy[j] = total_strain[3*j +2];
    	}

    	/* Position of the palette is controlled by float label_width in static void add_palette() of draw.c*/
    	//nb_mesh2D_distort_with_field(part, NB_NODE, total_displacement, 0.5);
	
	char* pltc_elems = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.elems.png"));
	nb_mesh2D_export_draw(part, strg_concat(pltc_elems, printable_name, "_pltc.elems.png"), 
				1200, 800, NB_ELEMENT, NB_CLASS, plastic_elements, true);
	nb_free_mem(pltc_elems);

	char* pltc_avDisp = nb_allocate_mem(strlen(printable_name) + strlen("_pltc.avg.disp.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_avDisp, printable_name, "_pltc.avg.disp.png"), 
				1200, 800, NB_NODE, NB_FIELD, avg_displacement, true);
	nb_free_mem(pltc_avDisp);

	char* pltc_dispX = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.disp.X.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_dispX, printable_name, "_pltc.disp.X.png"), 
				1200, 800, NB_NODE, NB_FIELD, total_displacement_x, true);
	nb_free_mem(pltc_dispX);
	
	char* pltc_dispY = nb_allocate_mem(strlen(printable_name) + sizeof("pltc.disp.Y.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_dispY, printable_name, "_pltc.disp.Y.png"), 
				1200, 800, NB_NODE, NB_FIELD, total_displacement_y, true);
	nb_free_mem(pltc_dispY);
	
	char* pltc_Sx = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Sx.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Sx, printable_name, "_pltc.Sx.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Sx, true);
	nb_free_mem(pltc_Sx);
	
	char* pltc_Sy = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Sy.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Sy, printable_name, "_pltc.Sy.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Sy, true);
	nb_free_mem(pltc_Sy);
	
	char* pltc_Sxy = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Sxy.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Sxy, printable_name, "_pltc.Sxy.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Sxy, true);
	nb_free_mem(pltc_Sxy);
	
	char* pltc_Ex = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Ex.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Ex, printable_name, "_pltc.Ex.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Ex, true);
	nb_free_mem(pltc_Ex);
	
	char* pltc_Ey = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Ey.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Ey, printable_name, "_pltc.Ey.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Ey, true);
	nb_free_mem(pltc_Ey);
	
	char* pltc_Exy = nb_allocate_mem(strlen(printable_name) + sizeof("_pltc.Exy.png"));
    	nb_mesh2D_export_draw(part, strg_concat(pltc_Exy, printable_name, "_pltc.Exy.png"), 
				1200, 800, NB_ELEMENT, NB_FIELD, Exy, true);
	nb_free_mem(pltc_Exy);

    	nb_free_mem(memblock);
}

char* strg_concat(char* print_name, char* problem_name, char* extension) /* Function's referenece declared in static_plasticity2D.h*/
{
	strcpy(print_name, problem_name);
	strcat(print_name, extension);
	
	return print_name;
}
