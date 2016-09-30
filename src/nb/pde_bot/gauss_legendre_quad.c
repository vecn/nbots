#include <stdint.h>
#include <string.h>

#include "nb/pde_bot/gauss_legendre_quad.h"

static void set_1_pnt(nb_glquadrature_t *glq);
static void set_2_pnt(nb_glquadrature_t *glq);
static void set_3_pnt(nb_glquadrature_t *glq);
static void set_4_pnt(nb_glquadrature_t *glq);
static void set_5_pnt(nb_glquadrature_t *glq);
static void set_6_pnt(nb_glquadrature_t *glq);
static void set_7_pnt(nb_glquadrature_t *glq);
static void set_8_pnt(nb_glquadrature_t *glq);
static void set_9_pnt(nb_glquadrature_t *glq);
static void set_10_pnt(nb_glquadrature_t *glq);

void nb_glquadrature_load(nb_glquadrature_t *glq, uint8_t N_points)
{  
	glq->N = N_points;

	switch (N_points) {
	case 1:
		set_1_pnt(glq);
		break;
	case 2:
		set_2_pnt(glq);
		break;
	case 3:
		set_3_pnt(glq);
		break;
	case 4:
		set_4_pnt(glq);
		break;
	case 5:
		set_5_pnt(glq);
		break;
	case 6:
		set_6_pnt(glq);
		break;
	case 7:
		set_7_pnt(glq);
		break;
	case 8:
		set_8_pnt(glq);
		break;
	case 9:
		set_9_pnt(glq);
		break;
	case 10:
		set_10_pnt(glq);
		break;
	default:
		glq->N = 0;
	}
}

static void set_1_pnt(nb_glquadrature_t *glq)
{
	glq->w[0] = 2.0;
	glq->x[0] = 0.0;
}

static void set_2_pnt(nb_glquadrature_t *glq)
{
	double w[2] = {1.0, 1.0};
	memcpy(glq->w, w, 2 * sizeof(*(glq->w)));

	double x[2] = {-0.5773502691896257,
		       0.5773502691896257};
	memcpy(glq->x, x, 2 * sizeof(*(glq->x)));
}

static void set_3_pnt(nb_glquadrature_t *glq)
{
	double w[3] = {0.5555555555555556,
		       0.8888888888888888,
		       0.5555555555555556};
	memcpy(glq->w, w, 3 * sizeof(*(glq->w)));

	double x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
	memcpy(glq->x, x, 3 * sizeof(*(glq->x)));
}

static void set_4_pnt(nb_glquadrature_t *glq)
{
	double w[4] = {0.6521451548625461, 0.6521451548625461,
		       0.3478548451374538, 0.3478548451374538};
	memcpy(glq->w, w, 4 * sizeof(*(glq->w)));

	double x[4] = {-0.3399810435848563, 0.3399810435848563,
		       -0.8611363115940526, 0.8611363115940526};
	memcpy(glq->x, x, 4 * sizeof(*(glq->x)));
}

static void set_5_pnt(nb_glquadrature_t *glq)
{
	double w[5] = {0.5688888888888889, 0.4786286704993665,
		       0.4786286704993665, 0.2369268850561891,
		       0.2369268850561891};
	memcpy(glq->w, w, 5 * sizeof(*(glq->w)));

	double x[5] = {0.0,
		       -0.5384693101056831, 0.5384693101056831,
		       -0.9061798459386640, 0.9061798459386640};
	memcpy(glq->x, x, 5 * sizeof(*(glq->x)));
}

static void set_6_pnt(nb_glquadrature_t *glq)
{
	double w[6] = {0.3607615730481386, 0.3607615730481386,
		       0.4679139345726910, 0.4679139345726910,
		       0.1713244923791704, 0.1713244923791704};
	memcpy(glq->w, w, 6 * sizeof(*(glq->w)));

	double x[6] = {0.6612093864662645, -0.6612093864662645,
		       -0.2386191860831969, 0.2386191860831969,
		       -0.9324695142031521, 0.9324695142031521};
	memcpy(glq->x, x, 6 * sizeof(*(glq->x)));
}

static void set_7_pnt(nb_glquadrature_t *glq)
{
	double w[7] = {0.4179591836734694, 0.3818300505051189,
		       0.3818300505051189, 0.2797053914892766,
		       0.2797053914892766, 0.1294849661688697,
		       0.1294849661688697};
	memcpy(glq->w, w, 7 * sizeof(*(glq->w)));

	double x[7] = {0.0,
		       0.4058451513773972, -0.4058451513773972,
		       -0.7415311855993945, 0.7415311855993945,
		       -0.9491079123427585, 0.9491079123427585};
	memcpy(glq->x, x, 7 * sizeof(*(glq->x)));
}

static void set_8_pnt(nb_glquadrature_t *glq)
{
	double w[8] = {0.3626837833783620, 0.3626837833783620,
		       0.3137066458778873, 0.3137066458778873,
		       0.2223810344533745, 0.2223810344533745,
		       0.1012285362903763, 0.1012285362903763};
	memcpy(glq->w, w, 8 * sizeof(*(glq->w)));

	double x[8] = {-0.1834346424956498, 0.1834346424956498,
		       -0.5255324099163290, 0.5255324099163290,
		       -0.7966664774136267, 0.7966664774136267,
		       -0.9602898564975363, 0.9602898564975363};
	memcpy(glq->x, x, 8 * sizeof(*(glq->x)));
}

static void set_9_pnt(nb_glquadrature_t *glq)
{
	double w[9] = {0.3302393550012598, 0.1806481606948574,
		       0.1806481606948574, 0.0812743883615744,
		       0.0812743883615744, 0.3123470770400029,
		       0.3123470770400029, 0.2606106964029354,
		       0.2606106964029354};
	memcpy(glq->w, w, 9 * sizeof(*(glq->w)));

	double x[9] = {0.0,
		       -0.8360311073266358, 0.8360311073266358,
		       -0.9681602395076261, 0.9681602395076261,
		       -0.3242534234038089, 0.3242534234038089,
		       -0.6133714327005904, 0.6133714327005904};
	memcpy(glq->x, x, 9 * sizeof(*(glq->x)));
}

static void set_10_pnt(nb_glquadrature_t *glq)
{
	double w[10] = {0.2955242247147529, 0.2955242247147529,
			0.2692667193099963, 0.2692667193099963,
			0.2190863625159820, 0.2190863625159820,
			0.1494513491505806, 0.1494513491505806,
			0.0666713443086881, 0.0666713443086881};
	memcpy(glq->w, w, 10 * sizeof(*(glq->w)));

	double x[10] = {-0.1488743389816312, 0.1488743389816312,
			-0.4333953941292472, 0.4333953941292472,
			-0.6794095682990244, 0.6794095682990244,
			-0.8650633666889845, 0.8650633666889845,
			-0.9739065285171717, 0.9739065285171717};
	memcpy(glq->x, x, 10 * sizeof(*(glq->x)));
}
