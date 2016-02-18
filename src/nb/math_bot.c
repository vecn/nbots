/******************************************************************************
 *   Math Cat: Math utilities.                                                *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "nb/math_bot.h"

inline int vcn_math_pow2i(int a)
{
	return a*a;
}

inline uint32_t vcn_math_pow2u(uint32_t a)
{
	return a*a;
}

inline double vcn_math_pow2(double a)
{
	return a*a;
}

inline double vcn_math_pow3(double a)
{
	return a*a*a;
}

inline double vcn_math_pow4(double a)
{
	double b = vcn_math_pow2(a);
	return vcn_math_pow2(b);
}

inline double vcn_math_pow5(double a)
{
	double b = vcn_math_pow2(a);
	double c = vcn_math_pow3(a);
	return b*c;
}
  
inline double vcn_math_pow6(double a)
{
	double b = vcn_math_pow3(a);
	return vcn_math_pow2(b);
}

inline double vcn_math_pow7(double a)
{
	double b = vcn_math_pow2(a);
	double c = vcn_math_pow3(a);
	return vcn_math_pow2(b)*c;
}

inline double vcn_math_pow8(double a)
{
	double b = vcn_math_pow2(a);
	double c = vcn_math_pow2(b);
	return vcn_math_pow2(c);
}

inline double vcn_math_pow9(double a)
{
	double b = vcn_math_pow3(a);
	return vcn_math_pow3(b);
}

inline int vcn_math_min(int a, int b)
{
	return (a < b)? a : b;
}

inline int vcn_math_max(int a, int b)
{
	return (a > b)? a : b;
}

inline uint32_t vcn_math_minu(uint32_t a, uint32_t b)
{
	return (a < b)? a : b;
}

inline uint32_t vcn_math_maxu(uint32_t a, uint32_t b)
{
	return (a > b)? a : b;
}

inline double vcn_math_mind(double a, double b)
{
	return (a < b)? a : b;
}

inline double vcn_math_maxd(double a, double b)
{
	return (a > b)? a : b;
}

int vcn_math_powk(int a, uint32_t k)
{
	int out = 1;
	for(uint32_t i = 0; i < k; i++)
		out *= a;
	return out;
}


double vcn_math_hypo(double a, double b)
{
	/* Compute the hypotenuse with non-destructive overflow or underflow */
	a = fabs(a);
	b = fabs(b);
	double hypo;
	if (a > b) {
		hypo = a * sqrt(1.0 + vcn_math_pow2(b/a));
	} else {
		double val = b * sqrt(1.0 + vcn_math_pow2(a/b));
		hypo = (b < 1e-9)? 0.0 : val;
	}
	return hypo;
}

inline double vcn_math_harmonic_avg(double a, double b)
{
	return (2 * a * b) / (a + b);
}

inline double vcn_math_log2(double arg)
{
	return log(arg) / NB_MATH_LOG2;
}

vcn_Gauss_Legendre_table_t* vcn_GLtable_create(uint32_t N_points)
{
	if(N_points == 0)
		return NULL;
	if(N_points > 10)
		return NULL;

	vcn_Gauss_Legendre_table_t *table = calloc(1, sizeof(*table));
  
	table->N = N_points;
	table->w = (double*) calloc(table->N, sizeof(double));
	table->x = (double*) calloc(table->N, sizeof(double));

	switch(N_points){
	case 1:
		/* Set weights */
		table->w[0] = 2.0;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.0;
		break;
	case 2:
		/* Set weights */
		table->w[0] = 1.0;
		table->w[1] = 1.0;
		/* Set roots of Legendre polynomials */
		table->x[0] = -0.5773502691896257;
		table->x[1] = 0.5773502691896257;
		break;
	case 3:
		/* Set weights */
		table->w[0] = 0.8888888888888888;
		table->w[1] = 0.5555555555555556;
		table->w[2] = 0.5555555555555556;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.0;
		table->x[1] = -0.7745966692414834;
		table->x[2] = 0.7745966692414834;
		break;
	case 4:
		/* Set weights */
		table->w[0] = 0.6521451548625461;
		table->w[1] = 0.6521451548625461;
		table->w[2] = 0.3478548451374538;
		table->w[3] = 0.3478548451374538;
		/* Set roots of Legendre polynomials */
		table->x[0] = -0.3399810435848563;
		table->x[1] = 0.3399810435848563;
		table->x[2] = -0.8611363115940526;
		table->x[3] = 0.8611363115940526;
		break;
	case 5:
		/* Set weights */
		table->w[0] = 0.5688888888888889;
		table->w[1] = 0.4786286704993665;
		table->w[2] = 0.4786286704993665;
		table->w[3] = 0.2369268850561891;
		table->w[4] = 0.2369268850561891;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.0;
		table->x[1] = -0.5384693101056831;
		table->x[2] = 0.5384693101056831;
		table->x[3] = -0.9061798459386640; 
		table->x[4] = 0.9061798459386640;
		break;
	case 6:
		/* Set weights */
		table->w[0] = 0.3607615730481386;
		table->w[1] = 0.3607615730481386;
		table->w[2] = 0.4679139345726910;
		table->w[3] = 0.4679139345726910;
		table->w[4] = 0.1713244923791704;
		table->w[5] = 0.1713244923791704;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.6612093864662645;
		table->x[1] = -0.6612093864662645;
		table->x[2] = -0.2386191860831969;
		table->x[3] = 0.2386191860831969;
		table->x[4] = -0.9324695142031521; 
		table->x[5] = 0.9324695142031521;
		break;
	case 7:
		/* Set weights */
		table->w[0] = 0.4179591836734694;
		table->w[1] = 0.3818300505051189;
		table->w[2] = 0.3818300505051189;
		table->w[3] = 0.2797053914892766;
		table->w[4] = 0.2797053914892766;
		table->w[5] = 0.1294849661688697;
		table->w[6] = 0.1294849661688697;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.0;
		table->x[1] = 0.4058451513773972;
		table->x[2] = -0.4058451513773972;
		table->x[3] = -0.7415311855993945;
		table->x[4] = 0.7415311855993945;
		table->x[5] = -0.9491079123427585; 
		table->x[6] = 0.9491079123427585;
		break;
	case 8:
		/* Set weights */
		table->w[0] = 0.3626837833783620;
		table->w[1] = 0.3626837833783620;
		table->w[2] = 0.3137066458778873;
		table->w[3] = 0.3137066458778873;
		table->w[4] = 0.2223810344533745;
		table->w[5] = 0.2223810344533745;
		table->w[6] = 0.1012285362903763;
		table->w[7] = 0.1012285362903763;
		/* Set roots of Legendre polynomials */
		table->x[0] = -0.1834346424956498;
		table->x[1] = 0.1834346424956498;
		table->x[2] = -0.5255324099163290;
		table->x[3] = 0.5255324099163290;
		table->x[4] = -0.7966664774136267; 
		table->x[5] = 0.7966664774136267;
		table->x[6] = -0.9602898564975363; 
		table->x[7] = 0.9602898564975363;
		break;
	case 9:
		/* Set weights */
		table->w[0] = 0.3302393550012598;
		table->w[1] = 0.1806481606948574;
		table->w[2] = 0.1806481606948574;
		table->w[3] = 0.0812743883615744;
		table->w[4] = 0.0812743883615744;
		table->w[5] = 0.3123470770400029;
		table->w[6] = 0.3123470770400029;
		table->w[7] = 0.2606106964029354;
		table->w[8] = 0.2606106964029354;
		/* Set roots of Legendre polynomials */
		table->x[0] = 0.0;
		table->x[1] = -0.8360311073266358;
		table->x[2] = 0.8360311073266358;
		table->x[3] = -0.9681602395076261;
		table->x[4] = 0.9681602395076261;
		table->x[5] = -0.3242534234038089;
		table->x[6] = 0.3242534234038089;
		table->x[7] = -0.6133714327005904;
		table->x[8] = 0.6133714327005904;
		break;
	case 10:
		/* Set weights */
		table->w[0] = 0.2955242247147529;
		table->w[1] = 0.2955242247147529;
		table->w[2] = 0.2692667193099963;
		table->w[3] = 0.2692667193099963;
		table->w[4] = 0.2190863625159820;
		table->w[5] = 0.2190863625159820;
		table->w[6] = 0.1494513491505806;
		table->w[7] = 0.1494513491505806;
		table->w[8] = 0.0666713443086881;
		table->w[9] = 0.0666713443086881;
		/* Set roots of Legendre polynomials */
		table->x[0] = -0.1488743389816312;
		table->x[1] = 0.1488743389816312;
		table->x[2] = -0.4333953941292472;
		table->x[3] = 0.4333953941292472;
		table->x[4] = -0.6794095682990244;
		table->x[5] = 0.6794095682990244;
		table->x[6] = -0.8650633666889845;
		table->x[7] = 0.8650633666889845;
		table->x[8] = -0.9739065285171717;
		table->x[9] = 0.9739065285171717;
		break;
	}
	return table;
}

void vcn_GLtable_destroy(vcn_Gauss_Legendre_table_t* GLtable)
{
	free(GLtable->w);
	free(GLtable->x);
	free(GLtable);
}
