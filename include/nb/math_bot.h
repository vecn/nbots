/******************************************************************************
 *   Math Cat: Math utilities.                                                *
 *   2011-2015 Victor Eduardo Cardoso Nungaray                                *
 *   Twitter: @victore_cardoso                                                *
 *   email: victorc@cimat.mx                                                  *
 ******************************************************************************/

/**
 * @file math_cat.h
 * @brief Simple math routines.
 * @author Victor Eduardo Cardoso Nungaray
 * @n victorc@@cimat.mx
 * @n @@victore_cardoso
 * @date 10 August 2015
 */

#ifndef __VCN_MATH_CAT_H__
#define __VCN_MATH_CAT_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define VCN_MATH_PI (3.14159265358979323846264338327)
#define VCN_MATH_MAX_UINT32_T ((unsigned int)(-1))
#define VCN_MATH_SQRT2 (1.41421356237309504880)
#define VCN_MATH_SQRT3 (1.73205080756887729352)
#define VCN_MATH_INV_SQRT3 (0.57735026919)  /* 1/sqrt(3) */
#define VCN_MATH_INV_SQRT6 (0.40824829046)  /* 1/sqrt(6) */
#define VCN_MATH_LOG2 (0.69314718056)

	typedef struct {
		uint32_t N;
		double* w;
		double* x;
	} vcn_Gauss_Legendre_table_t;

  
	int vcn_math_pow2i(int a);
	int vcn_math_powk(int a, uint32_t k);
	double vcn_math_pow2(double a);
	double vcn_math_pow3(double a);
	double vcn_math_pow4(double a);
	double vcn_math_pow5(double a);
	double vcn_math_pow6(double a);
	double vcn_math_pow7(double a);
	double vcn_math_pow8(double a);
	double vcn_math_pow9(double a);
	int vcn_math_min(int a, int b);
	int vcn_math_max(int a, int b);
	uint32_t vcn_math_minu(uint32_t a, uint32_t b);
	uint32_t vcn_math_maxu(uint32_t a, uint32_t b);
	double vcn_math_mind(double a, double b);
	double vcn_math_maxd(double a, double b);
	double vcn_math_hypo(double a, double b);
	double vcn_math_harmonic_avg(double a, double b);
	double vcn_math_log2(double arg);  

	vcn_Gauss_Legendre_table_t* vcn_GLtable_create(uint32_t N_points);
	void vcn_GLtable_destroy(vcn_Gauss_Legendre_table_t* GLtable);

#ifdef __cplusplus
}
#endif

#endif
