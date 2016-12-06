#ifndef STATIC_PLASTIC2D_TEST_H_INCLUDED
#define STATIC_PLASTIC2D_TEST_H_INCLUDED

typedef struct  damage_results_t plastic_results_t;

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
void check_input_values(nb_material_t *material, nb_analysis2D_t analysis2D,
                    nb_analysis2D_params *params2D);

#endif // STATIC_PLASTIC2D_TEST_H_INCLUDED
