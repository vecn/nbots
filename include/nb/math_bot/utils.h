#ifndef __NB_MATH_BOT_UTILS_H__
#define __NB_MATH_BOT_UTILS_H__

#include <stdint.h>

#define NB_PI (3.14159265359)
#define NB_PHI (1.61803398875)

#define NB_MATH_PI (3.14159265358979323846264338327)
#define NB_MATH_MAX_UINT32_T ((unsigned int)(-1))
#define NB_MATH_SQRT2 (1.41421356237309504880)
#define NB_MATH_SQRT3 (1.73205080756887729352)
#define NB_MATH_INV_SQRT3 (0.57735026919)  /* 1/sqrt(3) */
#define NB_MATH_INV_SQRT6 (0.40824829046)  /* 1/sqrt(6) */
#define NB_MATH_LOG2 (0.69314718056)
  
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

#endif
