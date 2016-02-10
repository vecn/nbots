#include "math.h"
#include "nb/math_bot.h"

double eval_1(double x, double coeff)
{
  return coeff;
}

double eval_logn(double x, double coeff)
{
  return coeff * vcn_math_log2(x);
}

double eval_logn2(double x, double coeff)
{
  return coeff * vcn_math_log2(x*x);
}

double eval_logn3(double x, double coeff)
{
  return coeff * vcn_math_log2(x*x*x);
}

double eval_n(double x, double coeff)
{
  return coeff * x;
}

double eval_nlogn(double x, double coeff)
{
  return coeff * x * vcn_math_log2(x);
}

double eval_nlogn2(double x, double coeff)
{
  return coeff * x * vcn_math_log2(x*x);
}

double eval_nlogn3(double x, double coeff)
{
  return coeff * x * vcn_math_log2(x*x*x);
}

double eval_n2(double x, double coeff)
{
  return coeff * vcn_math_pow2(x);
}

double eval_n2logn(double x, double coeff)
{
  return coeff * vcn_math_pow2(x) * vcn_math_log2(x);
}

double eval_n2logn2(double x, double coeff)
{
  return coeff * vcn_math_pow2(x) * vcn_math_log2(x*x);
}

double eval_n2logn3(double x, double coeff)
{
  return coeff * vcn_math_pow2(x) * vcn_math_log2(x*x*x);
}

double eval_n3(double x, double coeff)
{
  return coeff * vcn_math_pow3(x);
}

double eval_n3logn(double x, double coeff)
{
  return coeff * vcn_math_pow3(x) * vcn_math_log2(x);
}

double eval_n3logn2(double x, double coeff)
{
  return coeff * vcn_math_pow3(x) * vcn_math_log2(x*x);
}

double eval_n3logn3(double x, double coeff)
{
  return coeff * vcn_math_pow3(x) * vcn_math_log2(x*x*x);
}

double eval_n4(double x, double coeff)
{
  return coeff * vcn_math_pow4(x);
}

double eval_expn(double x, double coeff)
{
  return coeff * exp(x);
}
