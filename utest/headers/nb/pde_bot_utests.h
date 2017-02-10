#ifndef __UTEST_NB_PDE_BOT_H__
#define __UTEST_NB_PDE_BOT_H__

void cunit_nb_pde_bot_ode_solver(void);
void cunit_nb_pde_bot_fem_sm_static_elasticity(void);
void cunit_nb_pde_bot_fem_sm_static_plasticity(void);
void cunit_nb_pde_bot_fem_sm_static_damage(void);
void cunit_nb_pde_bot_cvfa_sm_static_elasticity(void);
void cunit_nb_pde_bot_cvfa_sm_static_damage_phase_field(void);

static void cunit_suites_nb_pde_bot(void)
{
	//cunit_nb_pde_bot_ode_solver();
	//cunit_nb_pde_bot_fem_sm_static_elasticity();
	cunit_nb_pde_bot_fem_sm_static_plasticity();
	//cunit_nb_pde_bot_fem_sm_static_damage();
	//cunit_nb_pde_bot_cvfa_sm_static_elasticity();
	//cunit_nb_pde_bot_cvfa_sm_static_damage_phase_field();
}

#endif
