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
/*#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_damage2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "nb/pde_bot/finite_element/solid_mechanics/pipeline.h"*/

/*#define INPUTS_DIR "/home/david/Desktop/nbots/utest_plastic/static_plasticity2D_inputs"*/
#define INPUTS_DIR "../utest/sources/nb/pde_bot/damage_inputs"

#define POW2(a) ((a)*(a))

typedef struct {
	double *disp;
	double *strain;
	double *stress;
	bool *plastic_elements;
	double *damage;
    uint32_t N_vtx;
	uint32_t N_trg;
} damage_results_t;

static int suite_init(void);
static int suite_clean(void);


static void test_damage_f1(void);
static void test_damage_f2(void);

static void check_problem_results(const void *part,
				  const damage_results_t *results);
static void damage_run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const void*,
					   const damage_results_t*));
static int damage_simulate(const char *problem_data,
		    nb_mesh2D_t *part, damage_results_t *results,
		    uint32_t N_vtx);
static void damage_get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx);
static void damage_results_init(damage_results_t *results, uint32_t N_vtx, uint32_t N_trg);
static void damage_results_finish(damage_results_t *results);
static int damage_read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int damage_read_geometry(nb_cfreader_t *cfreader, nb_model_t *model);
static int damage_read_material(nb_cfreader_t *cfreader, nb_material_t *mat);
static int damage_read_plasticity2D_params(nb_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
void damage_check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D);
void print_damage_results(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                            double *total_strain, const nb_mesh2D_t *const part,
                            double *damage, const char *problem_data, const nb_fem_elem_t *const elem);
uint8_t nb_fem_compute_2D_Damage_Solid_Mechanics
			(const nb_mesh2D_t *const part,
			 const nb_fem_elem_t *const elem,
			 const nb_material_t *const material,
			 const nb_bcond_t *const bcond,
			 bool enable_self_weight,
			 double gravity[2],
			 bool enable_Cholesky_solver,
			 nb_analysis2D_t analysis2D,
			 nb_analysis2D_params *params2D,
			 nb_fem_implicit_t* params,
			 const char* logfile,
			 double *damage,
			 double *displacement,
			 double *strain,
			 double *stress);

void cunit_nb_pde_bot_fem_sm_static_elasticity(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/"\
			     "static_damage2D.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Damage failure 1", test_damage_f1);
	CU_add_test(suite, "Damage failure 2", test_damage_f2);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_damage_f1(void)
{
        damage_run_test("%s/damage_beam_failure_one.txt", 300,
        check_problem_results);
}
static void check_problem_results(const void *part,
				  const damage_results_t *results)
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

static void test_damage_f2(void)
{
        damage_run_test("%s/damage_beam_failure_one.txt", 300,
        check_problem_results);
}

static void damage_run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const void*,
					   const damage_results_t*))
{
	damage_results_t results;
	nb_mesh2D_t *part = nb_allocate_on_stack(nb_mesh2D_get_memsize(NB_TRIAN));
	nb_mesh2D_init(part, NB_TRIAN);

	int status = damage_simulate(problem_data, part, &results, N_vtx);

    printf("Simulation status: %i\n", status); /* TEMPORAL */
	if (status == 0) {
        printf("The simulation ran properly. \n");
	}
	check_results(part, &results);

	nb_mesh2D_finish(part);
	damage_results_finish(&results);
}

static int damage_simulate(const char *problem_data,
		    nb_mesh2D_t *part, damage_results_t *results,
		    uint32_t N_vtx)
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
		damage_read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);
    printf("Read status: %i\n", read_status); /* TEMPORAL */
    damage_check_input_values(material, analysis2D, &params2D);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	damage_get_mesh(model, part, N_vtx);

	nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(part);
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	damage_results_init(results, N_nodes, N_elems);
	uint32_t N_force_steps = 100;
	double accepted_tol = 1;
	double gravity[2] = {0, 0}; /*Antes era {0, -9.81}*/
	bool *elements_enabled = nb_allocate_mem(N_elems*sizeof(elements_enabled));
	nb_fem_implicit_t *params = nb_fem_implicit_create();
	nb_fem_implicit_set_N_max_iter(params, 300);
	nb_fem_implicit_set_N_max_iter_without_enhance(params, 300);
	nb_fem_implicit_set_N_steps(params, 100);
	nb_fem_implicit_set_residual_tolerance(params, 1);

    status = nb_fem_compute_2D_Damage_Solid_Mechanics
                            (part, elem, material, bcond, true, gravity,
                             false, analysis2D, &params2D, params, problem_data, results->damage,
                             results->disp, results->strain, results->stress);
    print_damage_results(N_nodes, N_elems, results->disp, results->stress, results->strain,
                         part, results->damage, problem_data, elem);

CLEANUP_FEM:
	nb_fem_elem_destroy(elem);
CLEANUP_INPUT:
	nb_model_finish(model);
	nb_bcond_finish(bcond);
	nb_material_destroy(material);

	return status;
}

static void damage_get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx)
{
	uint32_t mesh_memsize = nb_tessellator2D_get_memsize();
	nb_tessellator2D_t* mesh = nb_allocate_on_stack(mesh_memsize);
	nb_tessellator2D_init(mesh);
	nb_tessellator2D_set_size_constraint(mesh,
				     NB_MESH_SIZE_CONSTRAINT_MAX_VTX,
				     N_vtx);
	nb_tessellator2D_set_geometric_constraint(mesh,
					  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
					  NB_GEOMETRIC_TOL);
	nb_tessellator2D_generate_from_model(mesh, model);

	nb_mesh2D_load_from_mesh(part, mesh);
	nb_tessellator2D_finish(mesh);
}

static void damage_results_init(damage_results_t *results, uint32_t N_vtx, uint32_t N_trg)
{
	uint32_t size_disp = N_vtx * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_trg * 3 * sizeof(*(results->strain));
	uint32_t size_damage = N_trg * sizeof(*(results->damage));
	uint32_t total_size = size_disp + 2 * size_strain + size_damage;
	char *memblock = nb_allocate_mem(total_size);

	results->N_vtx = N_vtx;
	results->N_trg = N_trg;
	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
	results->damage = (void*)(memblock + size_disp + 2 * size_strain);
	results->plastic_elements = NULL;
}

static inline void damage_results_finish(damage_results_t *results)
{
	nb_free_mem(results->disp);
}

static int damage_read_problem_data
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
	if (0 != damage_read_geometry(cfreader, model)) {
		printf("\nERROR: Geometry contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != nb_bcond_read(bcond, cfreader)) {
		printf("\nERROR: Boundary C. contain errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != damage_read_material(cfreader, mat)) {
		printf("\nERROR: Material contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != damage_read_plasticity2D_params(cfreader, analysis2D, params2D)) {
		printf("\nERROR: Reading problem type in %s.\n",
		       filename);
		goto EXIT;
	}
	status = 0;
EXIT:
    nb_cfreader_destroy(cfreader);
	return status;
}

static int damage_read_geometry(nb_cfreader_t *cfreader, nb_model_t *model)
{
	int status = 1;
	/* Read model vertices */
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

	if (1 > N) {
		goto CLEANUP_VERTICES;}
	model->M = N;
	model->edge = malloc(2 * model->M * sizeof(*model->edge));
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

static int damage_read_material(nb_cfreader_t *cfreader, nb_material_t *mat)
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

	double traction_limit_stress;
	if (0 != nb_cfreader_read_double(cfreader, &traction_limit_stress))
		goto EXIT;
	nb_material_set_traction_limit_stress(mat, traction_limit_stress);

	double fracture_energy;
	if (0 != nb_cfreader_read_double(cfreader, &fracture_energy))
		goto EXIT;
	nb_material_set_fracture_energy(mat, fracture_energy);
	double density;
	if (0 != nb_cfreader_read_double(cfreader, &density));
	nb_material_set_density(mat, density);

	status = 0;
EXIT:
	return status;
}


static int damage_read_plasticity2D_params(nb_cfreader_t *cfreader,
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

void damage_check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D) {
	printf("Elastic Modulus: %lf\n", nb_material_get_elasticity_module(material)); /* TEMPORAL */
    printf("Traction limit stress: %lf\n", nb_material_get_traction_limit_stress(material)); /* TEMPORAL */
    printf("Fracture energy: %lf\n", nb_material_get_fracture_energy(material));
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

void print_damage_results(uint32_t N_nod, uint32_t N_elem, double *total_displacement, double *stress,
                            double *total_strain, const nb_mesh2D_t *const part,
                            double *damage, const char *problem_data, const nb_fem_elem_t *const elem)
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

    /*****************************************************************/
	/******************** > Print Results ****************************/
	/*****************************************************************/
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
    uint8_t N_gp = nb_fem_elem_get_N_gpoints(elem);
    for(int j = 0; j < N_elem; j++) {
        for(int k = 0; k < N_gp; k++) {
        Sx[j + k] = stress[3*j + k];
        Sy[j + k] = stress[3*j + k + 1];
        Sxy[j + k] = stress[3*j + k + 2];
        }
    }
    double *strainX = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size);
    double *strainY = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size + Ex_size);
    double *strainXY = (void*)(memblock + avg_displacement_size + total_displacement_x_size + total_displacement_y_size + Sx_size +
                          Sy_size + Sxy_size + Ex_size + Ey_size);
    for(int i = 0; i < N_elem; i++){
        for(int j = 0; j < N_gp; j++){
            strainX[i + j] = total_strain[3*i + j];
            strainY[i + j] = total_strain[3*i + j + 1];
            strainXY[i + j] = total_strain[3*i + j + 2];
        }
    }
    /* Position of the palette is controlled by float label_width in static void add_palette() of draw.c*/
    nb_mesh2D_distort_with_field(part, NB_NODE, total_displacement, 0.01);
    nb_mesh2D_export_draw(part, "DMG_damaged_elements.png", 1200, 800, NB_ELEMENT, NB_FIELD, damage, true);
    nb_mesh2D_export_draw(part,"DMG_avg_disp.png", 1200, 800, NB_NODE, NB_FIELD, avg_displacement, true);
    nb_mesh2D_export_draw(part,"DMG_disp_x.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_x, true);
    nb_mesh2D_export_draw(part,"DMG_disp_y.png", 1200, 800, NB_NODE, NB_FIELD, total_displacement_y, true);
    nb_mesh2D_export_draw(part, "DMG_strainX.png", 1200, 800, NB_ELEMENT, NB_FIELD, strainX, true);
    nb_mesh2D_export_draw(part, "DMG_strainY.png", 1200, 800, NB_ELEMENT, NB_FIELD, strainY, true);
    nb_mesh2D_export_draw(part, "DMG_strainXY.png", 1200, 800, NB_ELEMENT, NB_FIELD, strainXY, true);
    nb_mesh2D_export_draw(part,"DMG_StressX.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sx, true);
    nb_mesh2D_export_draw(part,"DMG_StressY.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sy, true);
    nb_mesh2D_export_draw(part,"DMG_StressXY.png", 1200, 800, NB_ELEMENT, NB_FIELD, Sxy, true);

    nb_free_mem(memblock);
}