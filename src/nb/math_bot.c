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
