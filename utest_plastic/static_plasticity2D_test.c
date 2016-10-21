#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "nb/memory_bot.h"
#include "nb/container_bot.h"
#include "nb/cfreader_cat.h"
#include "nb/geometric_bot.h"
#include "nb/pde_bot.h"
#include "nb/pde_bot/common_solid_mechanics/analysis2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/static_plasticity2D.h"
#include "nb/pde_bot/finite_element/solid_mechanics/set_bconditions.h"
#include "nb/pde_bot/finite_element/solid_mechanics/pipeline.h"

#define INPUTS_DIR "C:/Users/David/Desktop/nbots/utest_plastic/static_plasticity2D_inputs"

#define POW2(a) ((a)*(a))

typedef struct {
	double *disp;
	double *strain;
	double *stress;
    uint32_t N_vtx;
	uint32_t N_trg;
} plastic_results_t;

static void test_beam_cantilever(void);
static void check_beam_cantilever(const void *part,
				  const plastic_results_t *results);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const void*,
					   const plastic_results_t*));
static int simulate(const char *problem_data,
		    nb_mesh2D_t *part, plastic_results_t *results,
		    uint32_t N_vtx);
static void get_mesh(const nb_model_t *model, void *part,
		     uint32_t N_vtx);
static void results_init(plastic_results_t *results, uint32_t N_vtx, uint32_t N_trg);
static void results_finish(plastic_results_t *results);
static int read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(nb_cfreader_t *cfreader, nb_model_t *model);
static int read_material(nb_cfreader_t *cfreader, nb_material_t *mat);
static int read_plasticity2D_params(nb_cfreader_t *cfreader,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static int read_2Dthickness(nb_cfreader_t *cfreader, nb_analysis2D_params *params2D);
void check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D);

int main() {

    test_beam_cantilever();

    return 0;
}

static void test_beam_cantilever(void)
{
        run_test("%s/plastic_beam_cantilever.txt", 1000,
        check_beam_cantilever);
}
static void check_beam_cantilever(const void *part,
				  const plastic_results_t *results)
{
    /* Hacer el check con tensiones porque aún no he conseguido los desplazamientos. Tensiones en X e Y */
	double max_X_stress = 0;
	double max_Y_stress = 0;
	double max_XY_stress = 0;
	uint32_t N_elem = nb_mesh2D_get_N_elems(part);
	for (uint32_t i = 0; i < N_elem; i++) {
        double X_stress = results->stress[3*i];
        double Y_stress = results->stress[3*i+1];
        double XY_stress = results->stress[3*i+2];
		if (max_X_stress < X_stress)
			max_X_stress = X_stress;
		if (max_Y_stress < Y_stress)
			max_Y_stress = Y_stress;
		if (max_XY_stress < XY_stress)
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
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     void (*check_results)(const void*,
					   const plastic_results_t*))
{
	plastic_results_t results;
	nb_mesh2D_t *part = nb_allocate_on_stack(nb_mesh2D_get_memsize(NB_TRIAN));
	nb_mesh2D_init(part, NB_TRIAN);

	int status = simulate(problem_data, part, &results, N_vtx);
    printf("Simulation status: %i\n", status); /* TEMPORAL */
	if (status = 0) {
        printf("The simulation ran properly. \n");
	}
	check_results(part, &results);

	nb_mesh2D_finish(part);
	results_finish(&results);
}

static int simulate(const char *problem_data,
		    nb_mesh2D_t *part, plastic_results_t *results,
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
		read_problem_data(input, model, bcond, material,
				  &analysis2D, &params2D);
    printf("Read status: %i\n", read_status);
    check_input_values(material, analysis2D, &params2D);

	if (0 != read_status)
		goto CLEANUP_INPUT;

	get_mesh(model, part, N_vtx);

	nb_fem_elem_t* elem = nb_fem_elem_create(NB_TRG_LINEAR);

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(part);
	uint32_t N_elems = nb_mesh2D_get_N_elems(part);
	results_init(results, N_nodes, N_elems);
	uint32_t N_force_steps = 100;
	double accepted_tol = 1e-3;
	double gravity[2] = {0, -9.81};
	bool *elements_enabled = nb_allocate_mem(N_elems*sizeof(elements_enabled));
    int status_fem = fem_compute_plastic_2D_Solid_Mechanics (part, elem, material,
                                                             bcond, true, gravity,
                                                             analysis2D, &params2D, elements_enabled, results->strain, results->stress,
                                                             results->disp, N_force_steps,
                                                             accepted_tol);
	if (0 != status_fem)
		goto CLEANUP_FEM;
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

static void results_init(plastic_results_t *results, uint32_t N_vtx, uint32_t N_trg)
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
    if (0 != read_2Dthickness(cfreader, params2D)) {
        printf("\nError: Reading thickness of the problem in %s.\n",
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
	return status;

nb_cfreader_destroy(cfreader);
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
	/*
	if (0 != nb_cfreader_read_double(cfreader, &(params2D->thickness)))
		goto EXIT;
	*/
	status = 0;
EXIT:
	return status;
}

static int read_2Dthickness(nb_cfreader_t *cfreader, nb_analysis2D_params *params2D) {

    int status = 1;
    if (0 != nb_cfreader_read_double(cfreader, &(params2D->thickness)))
    goto EXIT;

    status = 0;
    EXIT:
        return status;
}

void check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D) {
	printf("Elastic Modulus: %lf\n", nb_material_get_elasticity_module(material)); /* TEMPORAL */
    printf("Plastic Modulus: %lf\n", nb_material_get_plasticity_module(material)); /* TEMPORAL */
    printf("Poisson module: %f\n", nb_material_get_poisson_module(material)); /* TEMPORAL */
    printf("Density: %lf\n", nb_material_get_density(material)); /* TEMPORAL */
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

