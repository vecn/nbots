#ifndef __UTEST_NB_GEOMETRIC_BOT_H__
#define __UTEST_NB_GEOMETRIC_BOT_H__

void cunit_nb_geometric_bot_utils2D(void);
void cunit_nb_geometric_bot_point2D(void);
void cunit_nb_geometric_bot_bins2D(void);
void cunit_nb_geometric_bot_bins2D_iterator(void);
void cunit_nb_geometric_bot_dewall(void);
void cunit_nb_geometric_bot_constrained_delaunay(void);
void cunit_nb_geometric_bot_mesh2D(void);
void cunit_nb_geometric_bot_mesh2D_image_density(void);
void cunit_nb_msh3trg(void);
void cunit_nb_mshquad(void);
void cunit_nb_mshpoly(void);
void cunit_nb_geometric_bot_model2D_clipper(void);
void cunit_nb_geometric_bot_model2D_verifier(void);

static void cunit_suites_nb_geometric_bot(void)
{
	cunit_nb_geometric_bot_utils2D();
	cunit_nb_geometric_bot_point2D();
	cunit_nb_geometric_bot_bins2D();
	cunit_nb_geometric_bot_bins2D_iterator();
	cunit_nb_geometric_bot_dewall();
	cunit_nb_geometric_bot_constrained_delaunay();
	cunit_nb_geometric_bot_tessellator2D();
	cunit_nb_geometric_bot_tessellator2D_image_density();
	cunit_nb_msh3trg();
	cunit_nb_mshquad();
	cunit_nb_mshpoly();
	cunit_nb_geometric_bot_model2D_clipper();
	cunit_nb_geometric_bot_model2D_verifier();
}

#endif
