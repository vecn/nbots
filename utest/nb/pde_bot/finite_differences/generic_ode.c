#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <CUnit/Basic.h>

#include "nb/memory_bot.h"
#include "nb/pde_bot.h"

#define L_DAMAGE 0.02

#define POW2(a) ((a)*(a))

static int suite_init(void);
static int suite_clean(void);

static void test_Laplace_eq(void);
static double unity(double x);
static void test_heat_eq(void);
static void test_damage_eq(void);
static double damage_diff(double x);
static double damage_source(double x);
static void test_convection_eq(void);
static double convection(double x);

void cunit_nb_pde_bot_fd_generic_ode(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/finite_differences/"\
			     "generic_ode.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Laplace equation", test_Laplace_eq);
	CU_add_test(suite, "Heat equation", test_heat_eq);
	CU_add_test(suite, "Damage equation", test_damage_eq);
	CU_add_test(suite, "Diff-Convection equation", test_convection_eq);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_Laplace_eq(void)
{
	int N = 11;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, unity, NULL, NULL, NULL,
			   0.0, 1.0, sol);
	
	CU_ASSERT(fabs(sol[5]-0.5) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double unity(double x)
{
	return 1.0;
}

static void test_heat_eq(void)
{
	int N = 11;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, unity, NULL, NULL, unity,
			   0.0, 0.0, sol);

	CU_ASSERT(fabs(sol[5]-0.125) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static void test_damage_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, damage_diff, NULL, unity,
			   damage_source,
			   0.0, 0.0, sol);
	
	CU_ASSERT(fabs(sol[15]-6.401844e-1) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double damage_diff(double x)
{
	return POW2(L_DAMAGE);
}

static double damage_source(double x)
{
	return (fabs(x-0.5) < 1e-2)?1.0:0.0;
}

static void test_convection_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, unity, convection, NULL, unity,
			   0.0, 0.0, sol);
	FILE *fp = fopen("../../../TEMPORAL_ode.txt", "w"); /**/
	/* Gnuplot 
	   f(x,A,D,k,a,c)=A + c*(k/a**2) - B(A,D,k,a,c)*(k/a) + (c/a)*x + (B(A,D,k,a,c)-(k/a)-c*(k/a**2))*exp((a/k)*x)
	   B(A,D,k,a,c)= (D-A-(c/a)*(1+k/a)+(k/a)*(1+c/a)*exp(a/k)) / (exp(a/k) - k/a)

	 */
	for (int i = 0; i < N; i++) {
		double x = i * 1.0/(N-1);
		fprintf(fp, "%lf %lf\n", x, sol[i]);
	}                                                   /**/
	fclose(fp);                                         /**/
	
	CU_ASSERT(fabs(sol[15] - 4.936130e-02) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double convection(double x)
{
	return 10.0;
}
