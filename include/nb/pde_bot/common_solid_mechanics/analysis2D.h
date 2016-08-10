#ifndef __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_ANALYSIS2D_H__
#define __NB_PDE_BOT_FINITE_ELEMENT_SOLID_MECHANICS_ANALYSIS2D_H__

typedef enum {
	NB_PLANE_STRESS,
	NB_PLANE_STRAIN,
	NB_SOLID_OF_REVOLUTION
} nb_analysis2D_t;

typedef struct {
	double thickness;
	double revolution_axe;
	char revolution_coordinate; /* 0: x, 1: y */
} nb_analysis2D_params;

#endif
