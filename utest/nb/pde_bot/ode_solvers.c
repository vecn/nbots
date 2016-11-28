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

static void test_fd_laplace_eq(void);
static double unity(double x);
static void test_fd_heat_eq(void);
static void test_fd_diffusion_reaction_eq(void);
static double diff_reaction_diff(double x);
static double diff_reaction_source(double x);
static void test_fem_diffusion_reaction_eq(void);
static void test_fd_convection_eq(void);
static double fd_convection(double x);
static void test_afd_convection_eq(void);
static double afd_convection(double x);
static void test_afem_convection_eq(void);

void cunit_nb_pde_bot_ode_solver(void)
{
	CU_pSuite suite =
		CU_add_suite("nb/pde_bot/ode_solvers.c",
			     suite_init, suite_clean);
	CU_add_test(suite, "Laplace equation FD", test_fd_laplace_eq);
	CU_add_test(suite, "Heat equation FD", test_fd_heat_eq);
	CU_add_test(suite, "Diffusion + Reaction equation FD",
		    test_fd_diffusion_reaction_eq);
	CU_add_test(suite, "Diffusion + Reaction equation FEM",
		    test_fem_diffusion_reaction_eq);
	CU_add_test(suite, "Diffusion + Convection equation FD",
		    test_fd_convection_eq);
	CU_add_test(suite, "Diffusion + Convection equation (Analytic FD)",
		    test_afd_convection_eq);
	CU_add_test(suite, "Diffusion + Convection equation (Analytic FEM)",
		    test_afem_convection_eq);
}

static int suite_init(void)
{
	return 0;
}

static int suite_clean(void)
{
	return 0;
}

static void test_fd_laplace_eq(void)
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

static void test_fd_heat_eq(void)
{
	int N = 11;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, unity, NULL, NULL, unity,
			   0.0, 0.0, sol);

	CU_ASSERT(fabs(sol[5]-0.125) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static void test_fd_diffusion_reaction_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, diff_reaction_diff, NULL, unity,
			   diff_reaction_source,
			   0.0, 0.0, sol);
	
	CU_ASSERT(fabs(sol[15]-6.401844e-1) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double diff_reaction_diff(double x)
{
	return POW2(L_DAMAGE);
}

static double diff_reaction_source(double x)
{
	return (fabs(x-0.5) < 1e-2)?1.0:0.0;
}

static void test_fem_diffusion_reaction_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, diff_reaction_diff, NULL, unity,
			   diff_reaction_source,
			   0.0, 0.0, sol);
	
	CU_ASSERT(fabs(sol[15]-6.401844e-1) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static void test_fd_convection_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_fd_solve_ode_h2(N, 1, unity, fd_convection, NULL, NULL,
			   0.0, 1.0, sol);
	
	CU_ASSERT(fabs(sol[25] - 1.859008e-01) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double fd_convection(double x)
{
	return 10.0;
}

static void test_afd_convection_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_analytic_fd_solve_diffusion_convection(N, 1, unity,
						   afd_convection, NULL,
						   0.0, 1.0, sol);
	
	CU_ASSERT(fabs(sol[25] - 0.0) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}

static double afd_convection(double x)
{
	return 300.0;
}

static void test_afem_convection_eq(void)
{
	int N = 31;
	uint32_t memsize = N * sizeof(double);
	double *sol = nb_soft_allocate_mem(memsize);

	nb_analytic_fem_solve_diffusion_convection(N, 1, unity,
						    afd_convection, unity,
						    0.0, 0.0, sol);
	
	CU_ASSERT(fabs(sol[25] - 0.0) < 1e-8);

	nb_soft_free_mem(memsize, sol);
}
