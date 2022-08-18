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

#define INPUTS_DIR "../utest/sources/nb/pde_bot/static_elasticity2D_inputs"

#define POW2(a) ((a)*(a))
#define CHECK_ZERO(a) ((fabs(a)<1e-25)?1:(a))
#define MIN(a,b) (((a)<(b))?(a):(b))

typedef struct {
	uint32_t N_faces;
	uint32_t N_elems;
	double *disp;
	double *strain;
	double *stress;
	char *boundary_mask;
	bool allocated;
} results_t;

static int suite_init(void);
static int suite_clean(void);

static void test_beam_cantilever(void);
static void check_beam_cantilever(const void *mesh,
				  const results_t *results);
static void test_plate_with_hole(void);
static void check_plate_with_hole(const void *mesh,
				  const results_t *results);
static double get_error_avg_pwh(const void *mesh, const double *stress,
				const char *boundary_mask,
				double error[3]);
static double get_face_error_avg_pwh(const void *mesh, const double *stress,
				     uint32_t face_id, double error[3]);
static void get_face_avg_of_analytic_stress_pwh(const nb_mesh2D_t *mesh,
						uint32_t face_id,
						double stress_avg[3]);
static void get_analytic_stress_pwh(double x, double y, double stress[3]);
static void modify_bcond_pwh(const void *mesh, nb_bcond_t *bcond);
static void pwh_DC_SGM_cond(const double *x, double t, double *out);
static void pwh_BC_SGM_cond(const double *x, double t, double *out);
static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_mesh2D_type mesh_type,
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */);
static int simulate(const char *problem_data, nb_mesh2D_t *mesh,
		    results_t *results, uint32_t N_vtx,
		    void (*modify_bcond)(const void*,
					     nb_bcond_t*)/* Can be NULL */);
static void get_mesh(const nb_model_t *model, void *mesh,
		     uint32_t N_vtx);
static void results_init(results_t *results, uint32_t N_faces,
			 uint32_t N_elems);
static void results_finish(results_t *results);
static int read_problem_data
		(const char* filename,
		 nb_model_t *model,
		 nb_bcond_t* bcond,
		 nb_material_t* mat,
		 nb_analysis2D_t *analysis2D,
		 nb_analysis2D_params *params2D);
static int read_geometry(nb_cfreader_t *cfr, nb_model_t *model);
static int read_material(nb_cfreader_t *cfr, nb_material_t *mat);
static int read_elasticity2D_params(nb_cfreader_t *cfr,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D);
static void show_cvfa_error_msg(int cvfa_status);

void cunit_nb_pde_bot_cvfa_sm_static_elasticity(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_element/solid_mechanics/" \
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
	run_test("%s/beam_cantilever.txt", 4000, NB_QUAD,
		 check_beam_cantilever, NULL);
}

static void check_beam_cantilever(const void *mesh,
				  const results_t *results)
{
	double max_disp = 0;
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	for (uint32_t i = 0; i < N_elems; i++) {
		double disp2 = POW2(results->disp[i * 2]) +
			POW2(results->disp[i*2+1]);
		if (max_disp < disp2)
			max_disp = disp2;
	}
	max_disp = sqrt(max_disp);
	printf("--- MAX DISP: %e\n", max_disp); /* TEMPORAL */
	CU_ASSERT(fabs(max_disp - 1.000109e-1) < 1e-3);
}


static void test_plate_with_hole(void)
{
	run_test("%s/plate_with_hole.txt", 4000, NB_POLY,
		 check_plate_with_hole,
		 modify_bcond_pwh);
}

static void check_plate_with_hole(const void *mesh,
				  const results_t *results)
{
	double error[3];
	double avg_error = get_error_avg_pwh(mesh, results->stress,
					     results->boundary_mask,
					     error);

	printf("-- AVG XX ERROR: %e\n", error[0]); /* TEMPORAL */
	printf("-- AVG YY ERROR: %e\n", error[1]); /* TEMPORAL */
	printf("-- AVG XY ERROR: %e\n", error[2]); /* TEMPORAL */
	printf("-- AVG VM ERROR: %e\n", avg_error); /* TEMPORAL */
	printf("--        ELEMS: %i\n", nb_mesh2D_get_N_elems(mesh)); /* TEMP */
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	double length = 0;
	for (uint32_t i = 0; i < N_faces; i++)
		length += nb_mesh2D_edge_get_length(mesh, i);
	length /= N_faces;
	printf("--      Delta X: %e\n", length); /* TEMPORAL */
	CU_ASSERT(avg_error < 9.7e-3);
}

static double get_error_avg_pwh(const void *mesh,
				const double *stress,
				const char *boundary_mask,
				double error[3])
{
	FILE *fp = fopen("./stress.txt", "w");             /* TEMPORAL */
	fprintf(fp,                                               /* TEMPORAL */
		"# Ex Ey Exy Enx Eny Etx Ety |En| |Et| Enn Ent Etn Ett E1 E2 " \
		"Sx Sy Sxy Snx Sny Stx Sty |Sn| |St| Snn Snt Stn Stt S1 S2 " \
		"Ax Ay Axy Anx Any Atx Aty |An| |At| Ann Ant Atn Att A1 A2\n");/**/
	fclose(fp);                                               /* TEMPORAL */
	double avg = 0.0;
	memset(error, 0, 3 * sizeof(*error));
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		double face_err[3];
		if (!(boundary_mask[i])) {
			avg += get_face_error_avg_pwh(mesh, stress,
						      i, face_err);
			error[0] += face_err[0];
			error[1] += face_err[1];
			error[2] += face_err[2];
		}
		
	}
	error[0] /= N_faces;
	error[1] /= N_faces;
	error[2] /= N_faces;
	return avg / N_faces;
}

static void TEMPORAL3(const void *mesh, const double analytic_stress[3],
		      const double *stress, uint32_t face_id)
{
	FILE *fp = fopen("./stress.txt", "a");  /* TEMPORAL */
	double nf[2];                                  /* TEMPORAL */
	nb_mesh2D_edge_get_normal(mesh, face_id, nf);  /* TEMPORAL */
	uint32_t id = face_id;                         /* TEMPORAL */
	double error[15];                              /* TEMPORAL */
	error[0] = fabs((analytic_stress[0] - stress[id * 3]) /      /**/
			CHECK_ZERO(analytic_stress[0]));         /**/
	error[1] = fabs((analytic_stress[1] - stress[id*3+1]) /      /**/
			CHECK_ZERO(analytic_stress[1]));         /**/
	error[2] = fabs((analytic_stress[2] - stress[id*3+2]) /      /**/
			CHECK_ZERO(analytic_stress[2]));         /**/
	double Sn[2];                                  /* TEMPORAL */
	Sn[0] = stress[id * 3]*nf[0] + 0.5*stress[id*3+2]*nf[1]; /**/
	Sn[1] = 0.5*stress[id*3+2]*nf[0] + stress[id*3+1]*nf[1]; /**/
	double St[2];                                  /* TEMPORAL */
	St[0] = stress[id * 3]*nf[1] - 0.5*stress[id*3+2]*nf[0]; /**/
	St[1] = 0.5*stress[id*3+2]*nf[1] - stress[id*3+1]*nf[0]; /**/
	double mSn = sqrt(POW2(Sn[0]) + POW2(Sn[1]));  /* TEMPORAL */
	double mSt = sqrt(POW2(St[0]) + POW2(St[1]));  /* TEMPORAL */
	double Snn = Sn[0] * nf[0] + Sn[1] * nf[1];    /* TEMPORAL */
	double Snt = Sn[0] * nf[1] - Sn[1] * nf[0];    /* TEMPORAL */
	double Stn = St[0] * nf[0] + St[1] * nf[1];    /* TEMPORAL */
	double Stt = St[0] * nf[1] - St[1] * nf[0];    /* TEMPORAL */
	double main_s[2];                                        /**/
	nb_pde_get_main_stress(stress[id*3], stress[id*3+1],     /**/
			       stress[id*3+2], main_s);          /**/
	double An[2];                                  /* TEMPORAL */
	An[0] = analytic_stress[0]*nf[0] +                       /**/
		0.5 * analytic_stress[2]*nf[1];                  /**/
	An[1] = 0.5 * analytic_stress[2]*nf[0] +                 /**/
		analytic_stress[1]*nf[1];                        /**/
	double At[2];                                  /* TEMPORAL */
	At[0] = analytic_stress[0]*nf[1] -                       /**/
		0.5 * analytic_stress[2]*nf[0];                  /**/
	At[1] = 0.5 * analytic_stress[2]*nf[1] -                 /**/
		analytic_stress[1]*nf[0];                        /**/
	double mAn = sqrt(POW2(An[0]) + POW2(An[1]));            /**/
	double mAt = sqrt(POW2(At[0]) + POW2(At[1]));            /**/
	double Ann = An[0] * nf[0] + An[1] * nf[1];    /* TEMPORAL */
	double Ant = An[0] * nf[1] - An[1] * nf[0];    /* TEMPORAL */
	double Atn = At[0] * nf[0] + At[1] * nf[1];    /* TEMPORAL */
	double Att = At[0] * nf[1] - At[1] * nf[0];    /* TEMPORAL */
	double main_a[2];                                        /**/
	nb_pde_get_main_stress(analytic_stress[0],		 /**/
			       analytic_stress[1],		 /**/
			       analytic_stress[2], main_a);	 /**/
	error[3] = fabs((An[0] - Sn[0]) / CHECK_ZERO(An[0]));    /**/
	error[4] = fabs((An[1] - Sn[1]) / CHECK_ZERO(An[1]));    /**/
	error[5] = fabs((At[0] - St[0]) / CHECK_ZERO(At[0]));    /**/
	error[6] = fabs((At[1] - St[1]) / CHECK_ZERO(At[1]));    /**/
	error[7] = fabs((mAn - mSn)/CHECK_ZERO(mAn));		 /**/
	error[8] = fabs((mAt - mSt)/CHECK_ZERO(mAt));		 /**/
	error[9] = fabs((Ann - Snn)/CHECK_ZERO(Ann));		 /**/
	error[10] = fabs((Ant - Snt)/CHECK_ZERO(Ant));		 /**/
	error[11] = fabs((Atn - Stn)/CHECK_ZERO(Atn));		 /**/
	error[12] = fabs((Att - Stt)/CHECK_ZERO(Att));           /**/
	error[13] = fabs((main_a[0] - main_s[0])/CHECK_ZERO(main_a[0]));
	error[14] = fabs((main_a[1] - main_s[1])/CHECK_ZERO(main_a[1]));
	fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e " \
		"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e " \
		"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", /**/
		error[0], error[1], error[2], error[3], error[4], /**/
		error[5], error[6], error[7], error[8], /* TEMPORAL */
		error[9], error[10], error[11], error[12],        /**/
		error[13], error[14],				  /**/
		stress[id*3], stress[id*3+1], stress[id*3+2],     /**/
		Sn[0], Sn[1], St[0], St[1], mSn, mSt, Snn, Snt,	  /**/
		Stn, Stt, main_s[0], main_s[1],                   /**/
		analytic_stress[0], analytic_stress[1], /* TEMPORAL */
		analytic_stress[2],                     /* TEMPORAL */
		An[0], An[1], At[0], At[1], mAn, mAt, Ann, Ant,   /**/
		Atn, Att, main_a[0], main_a[1]);                  /**/
	fclose(fp);                                     /* TEMPORAL */
}

static double get_face_error_avg_pwh(const void *mesh, const double *stress,
				     uint32_t face_id, double error[3])
{
	
	double analytic_stress[3];
	get_face_avg_of_analytic_stress_pwh(mesh, face_id,
					    analytic_stress);
	
	error[0] = fabs(1.0 - stress[face_id * 3] / analytic_stress[0]);
	error[1] = fabs(1.0 - stress[face_id*3+1] / analytic_stress[1]);
	error[2] = fabs(1.0 - stress[face_id*3+2] / analytic_stress[2]);

	TEMPORAL3(mesh, analytic_stress, stress, face_id);

	double vm_stress = 
		nb_pde_get_vm_stress(stress[face_id * 3],
				     stress[face_id*3+1],
				     stress[face_id*3+2]);
	double analytic_vm_stress =
		nb_pde_get_vm_stress(analytic_stress[0],
				     analytic_stress[1],
				     analytic_stress[2]);
	double vm_error = fabs(1.0 - vm_stress / analytic_vm_stress);
	return vm_error;
}

static void get_face_avg_of_analytic_stress_pwh(const nb_mesh2D_t *mesh,
						uint32_t face_id,
						double stress_avg[3])
{
	int Np = 100;
	double lf = nb_mesh2D_edge_get_length(mesh, face_id);
	memset(stress_avg, 0, 3 * sizeof(*stress_avg));
	
	/* Get numerical integral over face */
	double step = lf / Np;
	double x1[2], x2[2];
	nb_mesh2D_edge_get_midpoint(mesh, face_id, 0, x1);
	for (int i = 0; i < Np; i++) {
		double w = step * (i+1);
		nb_mesh2D_edge_get_midpoint(mesh, face_id, w, x2);
		double s1[3], s2[3];
		get_analytic_stress_pwh(x1[0], x1[1], s1);
		get_analytic_stress_pwh(x2[0], x2[1], s2);
		stress_avg[0] += ((s1[0] + s2[0]) * step) / 2;
		stress_avg[1] += ((s1[1] + s2[1]) * step) / 2;
		stress_avg[2] += ((s1[2] + s2[2]) * step) / 2;
		memcpy(x1, x2, 2 * sizeof(*x1));
	}
	
	stress_avg[0] /= lf;
	stress_avg[1] /= lf;
	stress_avg[2] /= lf;
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

static void modify_bcond_pwh(const void *mesh,
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

static void pwh_DC_SGM_cond(const double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[0];
	out[1] = stress[2];
}

static void pwh_BC_SGM_cond(const double *x, double t, double *out)
{
	double stress[3];
	get_analytic_stress_pwh(x[0], x[1], stress);
	out[0] = stress[2];
	out[1] = stress[1];
}

static void TEMPORAL1(nb_mesh2D_t *mesh, results_t *results)
{
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	double *disp = nb_allocate_mem(N_elems * sizeof(*disp));

	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	double *disp_nodes = nb_allocate_mem(N_nodes * sizeof(*disp_nodes));

	nb_mesh2D_distort_with_field(mesh, NB_ELEMENT, results->disp, 0.5);

	for (uint32_t i = 0; i < N_elems; i++)
		disp[i] = results->disp[i*2];

	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, disp, disp_nodes);
	nb_mesh2D_export_draw(mesh, "./CVFA_dx.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      disp_nodes, 0, 0, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_elems; i++)
		disp[i] = results->disp[i*2+1];

	nb_mesh2D_extrapolate_elems_to_nodes(mesh, 1, disp, disp_nodes);
	nb_mesh2D_export_draw(mesh, "./CVFA_dy.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      disp_nodes, 0, 0, true);/* TEMPORAL */

	nb_free_mem(disp);
	nb_free_mem(disp_nodes);
}

static void TEMPORAL2(nb_mesh2D_t *mesh, results_t *results)
{
	uint32_t N_nodes = nb_mesh2D_get_N_nodes(mesh);
	uint32_t memsize = N_nodes * (4 * sizeof(double) + sizeof(uint16_t));
	char *memblock = nb_soft_allocate_mem(memsize);
	double *stress = (void*) memblock;
	double *vm_stress = (void*) (memblock + N_nodes* 3 * sizeof(double));
	uint16_t *counter = (void*) (memblock + N_nodes* 4 * sizeof(double));

	memset(stress, 0, 3 * N_nodes * sizeof(*stress));
	memset(counter, 0, N_nodes * sizeof(*counter));
	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	for (uint32_t i = 0; i < N_faces; i++) {
		uint32_t v1 = nb_mesh2D_edge_get_1n(mesh, i);
		uint32_t v2 = nb_mesh2D_edge_get_2n(mesh, i);
		if (!(results->boundary_mask[i])) {
			stress[v1 * 3] += results->stress[i * 3];
			stress[v1*3+1] += results->stress[i*3+1];
			stress[v1*3+2] += results->stress[i*3+2];
			counter[v1] += 1;

			stress[v2 * 3] += results->stress[i * 3];
			stress[v2*3+1] += results->stress[i*3+1];
			stress[v2*3+2] += results->stress[i*3+2];
			counter[v2] += 1;
		} else {
			double x = nb_mesh2D_node_get_x(mesh, v1);
			double y = nb_mesh2D_node_get_y(mesh, v1);
			double s[3];
			get_analytic_stress_pwh(x, y, s);

			stress[v1 * 3] += s[0];
			stress[v1*3+1] += s[1];
			stress[v1*3+2] += s[2];
			counter[v1] += 1;
			
			x = nb_mesh2D_node_get_x(mesh, v2);
			y = nb_mesh2D_node_get_y(mesh, v2);
			get_analytic_stress_pwh(x, y, s);

			stress[v2 * 3] += s[0];
			stress[v2*3+1] += s[1];
			stress[v2*3+2] += s[2];
			counter[v2] += 1;
		}
	}
	for (uint32_t i = 0; i < N_nodes; i++) {
		stress[i * 3] /= counter[i];
		stress[i*3+1] /= counter[i];
		stress[i*3+2] /= counter[i];
	}

	/*
	for (uint32_t i = 0; i < N_nodes; i++) {
		vm_stress[i] = nb_pde_get_vm_stress(stress[i * 3],
						    stress[i*3+1],
						    stress[i*3+2]);
	}
	*/

	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_mesh2D_node_get_x(mesh, i);
		double y = nb_mesh2D_node_get_y(mesh, i);
		double s[3];
		get_analytic_stress_pwh(x, y, s);
		vm_stress[i] = MIN(1,fabs((s[0]-stress[i*3])/CHECK_ZERO(s[0])));
	}

	nb_mesh2D_export_draw(mesh, "./CVFA_xErr.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_mesh2D_node_get_x(mesh, i);
		double y = nb_mesh2D_node_get_y(mesh, i);
		double s[3];
		get_analytic_stress_pwh(x, y, s);
		vm_stress[i] = MIN(1,fabs((s[1]-stress[i*3+1])/CHECK_ZERO(s[1])));
	}

	nb_mesh2D_export_draw(mesh, "./CVFA_yErr.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++) {
		double x = nb_mesh2D_node_get_x(mesh, i);
		double y = nb_mesh2D_node_get_y(mesh, i);
		double s[3];
		get_analytic_stress_pwh(x, y, s);
		vm_stress[i] = MIN(1,fabs((s[2]-stress[i*3+2])/CHECK_ZERO(s[2])));
	}

	nb_mesh2D_export_draw(mesh, "./CVFA_xyErr.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0, true);/* TEMPORAL */

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3];

	nb_mesh2D_export_draw(mesh, "./CVFA_Sxx.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0,  true);/* TEMPORAL */
	nb_mesh2D_export_level_sets(mesh, "./CVFA_ls_Sxx.eps",
				    1000, 800, vm_stress, 20, false);

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3+1];

	nb_mesh2D_export_draw(mesh, "./CVFA_Syy.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0, true);/* TEMPORAL */
	nb_mesh2D_export_level_sets(mesh, "./CVFA_ls_Syy.eps",
				    1000, 800, vm_stress, 20, false);

	for (uint32_t i = 0; i < N_nodes; i++)
		vm_stress[i] = stress[i*3+2];

	nb_mesh2D_export_draw(mesh, "./CVFA_Sxy.png", 1000, 800,
			      NB_NODE, NB_FIELD,
			      vm_stress, 0, 0, true);/* TEMPORAL */
	nb_mesh2D_export_level_sets(mesh, "./CVFA_ls_Sxy.eps",
				    1000, 800, vm_stress, 20, false);

	nb_soft_free_mem(memsize, memblock);
}

static void run_test(const char *problem_data, uint32_t N_vtx,
		     nb_mesh2D_type mesh_type,
		     void (*check_results)(const void*,
					   const results_t*),
		     void (*modify_bcond)(const void*,
					  nb_bcond_t*)/* Can be NULL */)
{
	results_t results;
	results.allocated = false;
	nb_mesh2D_t *mesh = nb_allocate_on_stack(nb_mesh2D_get_memsize(mesh_type));
	nb_mesh2D_init(mesh, mesh_type);

	int status = simulate(problem_data, mesh, &results,
			      N_vtx, modify_bcond);
	
	CU_ASSERT(0 == status);

	check_results(mesh, &results);

	TEMPORAL2(mesh, &results);
	TEMPORAL1(mesh, &results);

	nb_mesh2D_finish(mesh);
	results_finish(&results);
}

static int simulate(const char *problem_data, nb_mesh2D_t *mesh,
		    results_t *results, uint32_t N_vtx,
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
		modify_bcond(mesh, bcond);

	if (0 != read_status)
		goto CLEANUP;

	get_mesh(model, mesh, N_vtx);

	uint32_t N_faces = nb_mesh2D_get_N_edges(mesh);
	uint32_t N_elems = nb_mesh2D_get_N_elems(mesh);
	results_init(results, N_faces, N_elems);

	int status_cvfa = 
		nb_cvfa_compute_2D_Solid_Mechanics(mesh, material,
						   bcond,
						   false, NULL,
						   analysis2D,
						   &params2D,
						   results->disp,
						   results->strain,
						   results->boundary_mask);

	if (0 != status_cvfa) {
		show_cvfa_error_msg(status_cvfa);
		goto CLEANUP;
	}

	nb_cvfa_compute_stress_from_strain(mesh,  material, analysis2D,
					   results->strain,
					   results->stress);

	status = 0;
CLEANUP:
	nb_model_finish(model);
	nb_bcond_finish(bcond);
	nb_material_destroy(material);

	return status;
}

static void get_mesh(const nb_model_t *model, void *mesh,
		     uint32_t N_vtx)
{
	uint32_t t2d_memsize = nb_tessellator2D_get_memsize();
	nb_tessellator2D_t* t2d = nb_allocate_on_stack(t2d_memsize);
	nb_tessellator2D_init(t2d);
	nb_tessellator2D_set_size_constraint(t2d,
					     NB_MESH_SIZE_CONSTRAINT_MAX_TRG,
					     N_vtx);
	nb_tessellator2D_set_geometric_constraint(t2d,
						  NB_MESH_GEOM_CONSTRAINT_MAX_EDGE_LENGTH,
						  NB_GEOMETRIC_TOL);
	nb_tessellator2D_generate_from_model(t2d, model);

	nb_mesh2D_load_from_tessellator2D(mesh, t2d);
	nb_tessellator2D_finish(t2d);
	nb_mesh2D_export_draw(mesh, "./mesh.png", 1000, 800,
			      NB_NULL, NB_NULL, NULL, 0, 0, true);/* TEMPORAL */
	nb_cvfa_draw_integration_mesh(mesh, "./CVFA_alpha_x.eps",/*T*/
				      1000, 800);              /* TEMPORAL */
}

static void results_init(results_t *results, uint32_t N_faces,
			 uint32_t N_elems)
{
	uint32_t size_disp = N_elems * 2 * sizeof(*(results->disp));
	uint32_t size_strain = N_faces * 3 * sizeof(*(results->strain));
	uint32_t size_mask = N_faces * sizeof(*(results->boundary_mask));
	uint32_t total_size = size_disp + 2 * size_strain + size_mask;
	char *memblock = nb_allocate_mem(total_size);

	results->N_faces = N_faces;
	results->N_elems = N_elems;
	results->disp = (void*) memblock;
	results->strain = (void*)(memblock + size_disp);
	results->stress = (void*)(memblock + size_disp + size_strain);
	results->boundary_mask = (void*)(memblock + size_disp +
					 2 * size_strain);
	results->allocated = true;
}

static inline void results_finish(results_t *results)
{
	if (results->allocated)
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
	nb_cfreader_t* cfr = nb_cfreader_create();
	nb_cfreader_add_line_comment_token(cfr, "#");
	int status = nb_cfreader_open_file(cfr, filename);
	if (0 != status) {
		printf("\nERROR: Can not open file %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_geometry(cfr, model)) {
		printf("\nERROR: Geometry contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != nb_bcond_read(bcond, cfr)) {
		printf("\nERROR: Boundary C. contain errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_material(cfr, mat)) {
		printf("\nERROR: Material contains errors in %s.\n",
		       filename);
		goto EXIT;
	}
	if (0 != read_elasticity2D_params(cfr, analysis2D, params2D)) {
		printf("\nERROR: Reading numerical params in %s.\n",
		       filename);
		goto EXIT;
	}
	status = 0;
EXIT:
	nb_cfreader_close_file(cfr);
	nb_cfreader_destroy(cfr);
	return status;
}

static int read_geometry(nb_cfreader_t *cfr, nb_model_t *model)
{
	int status = 1;
	/* Read modele vertices */
	uint32_t N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto EXIT;

	if (1 > N)
		goto EXIT;
	model->N = N;

	model->vertex = nb_allocate_mem(2 * model->N * sizeof(*(model->vertex)));
	for (uint32_t i = 0; i < 2 * model->N; i++) {
		if (0 != nb_cfreader_read_double(cfr, &(model->vertex[i])))
			goto CLEANUP_VERTICES;
	}
	/* Read model segments */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto CLEANUP_VERTICES;

	if (1 > N)
		goto CLEANUP_VERTICES;
	model->M = N;

	model->edge = nb_allocate_mem(2 * model->M * sizeof(*model->edge));
	for (uint32_t i = 0; i < 2 * model->M; i++) {
		if (0 != nb_cfreader_read_uint(cfr, &(model->edge[i])))
			goto CLEANUP_SEGMENTS;
	}
	/* Read model holes */
	N = 0;
	if (0 != nb_cfreader_read_uint(cfr, &N))
		goto CLEANUP_SEGMENTS;
	model->H = N;

	model->holes = NULL;
	if (0 < model->H) {
		model->holes = nb_allocate_mem(2 * model->H *
					       sizeof(*(model->holes)));
		for (uint32_t i = 0; i < 2 * model->H; i++) {
			if (0 != nb_cfreader_read_double(cfr,
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


static int read_material(nb_cfreader_t *cfr, nb_material_t *mat)
{
	int status = 1;
	double poisson_module;
	if (0 != nb_cfreader_read_double(cfr, &poisson_module))
		goto EXIT;
	nb_material_set_poisson_module(mat, poisson_module);

	double elasticity_module;
	if (0 != nb_cfreader_read_double(cfr, &elasticity_module))
		goto EXIT;
	nb_material_set_elasticity_module(mat, elasticity_module);

	double fracture_energy;
	if (0 != nb_cfreader_read_double(cfr, &fracture_energy))
		goto EXIT;
	nb_material_set_fracture_energy(mat, fracture_energy);

	double compression_limit_stress;
	if (0 != nb_cfreader_read_double(cfr, &compression_limit_stress))
		goto EXIT;
	nb_material_set_compression_limit_stress(mat,
						      compression_limit_stress);

	double traction_limit_stress;
	if (0 != nb_cfreader_read_double(cfr, &traction_limit_stress))
		goto EXIT;
	nb_material_set_traction_limit_stress(mat, traction_limit_stress);
	status = 0;
EXIT:
	return status;
}


static int read_elasticity2D_params(nb_cfreader_t *cfr,
				    nb_analysis2D_t *analysis2D,
				    nb_analysis2D_params *params2D)
{
	int status = 1;
	int iaux;
	if (0 != nb_cfreader_read_int(cfr, &iaux))
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
	if (0 != nb_cfreader_read_double(cfr, &(params2D->thickness)))
		goto EXIT;
	status = 0;
EXIT:
	return status;
}

static void show_cvfa_error_msg(int cvfa_status)
{
	switch (cvfa_status) {
	case 1:
		printf("-- NB_ERROR: CVFA assembly fails.\n");
		break;
	case 2:
		printf("-- NB_ERROR: CVFA solver fails.\n");
		break;
	default:
		printf("-- NB_ERROR: CVFA fails.\n");
		break;
	}
}
