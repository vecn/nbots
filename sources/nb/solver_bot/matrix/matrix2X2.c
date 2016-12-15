#include <math.h>

#include "nb/math_bot.h"
#include "nb/solver_bot/matrix/matrix2X2.h"

#define POW2(a) ((a)*(a))

double nb_matrix_2X2_inverse(const double *const A, double* A_inv)
{
	double det = (A[0] * A[3] - A[1] * A[2]);
	A_inv[0] = A[3] / det;
	A_inv[1] = -A[1] / det;
	A_inv[2] = -A[2] / det;
	A_inv[3] = A[0] / det;
	return det;
}

void nb_matrix_2X2_eigen(const double *const A,
			  double* Lambda, /* Output (diag) */
			  double* P,      /* Output */
			  double tolerance)
/* Returns A decomposed into P Lambda P' */
{
	if (fabs(A[1]) < tolerance && fabs(A[2]) < tolerance) {
		Lambda[0] = A[0];
		Lambda[1] = A[3];
		P[0] = 1.0;
		P[2] = 0.0;
		P[1] = 0.0;
		P[3] = 1.0;
	} else {
		double T = A[0] + A[3];
		double D = A[0] * A[3] - A[1] * A[2];
		double root = sqrt(POW2(T)/4.0 - D);
		Lambda[0] = T/2.0 + root;
		Lambda[1] = T/2.0 - root;
		if (fabs(A[2]) > tolerance) {
			P[0] = Lambda[0] - A[3];
			P[2] = A[2];
			P[1] = Lambda[1] - A[3];
			P[3] = A[2];
		} else {
			P[0] = A[1];
			P[2] = Lambda[0] - A[0];
			P[1] = A[1];
			P[3] = Lambda[1] - A[0];
		}
	}
	if (fabs(Lambda[0]) < fabs(Lambda[1])) {
		double aux = Lambda[0];
		Lambda[0] = Lambda[1];
		Lambda[1] = aux;

		aux = P[0];
		P[0] = P[1];
		P[1] = aux;

		aux = P[2];
		P[2] = P[3];
		P[3] = aux;
	}
	/* Normalize eigenvectors */
	double normalizer = sqrt(POW2(P[0]) + POW2(P[2]));
	P[0] /= normalizer;
	P[2] /= normalizer;
	normalizer = sqrt(POW2(P[1]) + POW2(P[3]));
	P[1] /= normalizer;
	P[3] /= normalizer;
}
double nb_matrix_2X2_det(double *A)
{
	return A[0] * A[3] - A[1] * A[2];
}

double nb_matrix_2X2_inverse_destructive(double *A)
{
	double det = A[0] * A[3] - A[1] * A[2];
	double tmp = A[0];
	A[0] = A[3]/det;
	A[3] = tmp/det;
	A[1] = -A[1]/det;
	A[2] = -A[2]/det;
	return det;
}
