#ifndef __UTEST_NB_PDE_BOT_H__
#define __UTEST_NB_PDE_BOT_H__

void cunit_nb_pde_bot_fem_sm_static_elasticity(void);

static void cunit_suites_nb_pde_bot(void)
{
	cunit_nb_pde_bot_fem_sm_static_elasticity();
	cunit_nb_pde_bot_cvfa_sm_static_elasticity();
}

#endif
