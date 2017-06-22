#include <math.h>

#include "nb/math_bot.h"
#include "nb/solver_bot/matrix/mtxsym2X2.h"

#define POW2(a) ((a)*(a))

double nb_mtxsym_2X2_inverse(const double A[3], double A_inv[3])
{
	double det = nb_mtxsym_2X2_det(A);
	A_inv[0] = A[1] / det;
	A_inv[1] = A[0] / det;
	A_inv[2] = -A[2] / det;
	return det;
}

void nb_mtxsym_2X2_eigen(const double A[3], double Lambda[2],
			 double P[4], double tol)
/* Returns A decomposed into P Lambda P' */
{
	/* Get first eigenvectors */
	if (fabs(A[2]) < tol) {
		Lambda[0] = A[0];
		Lambda[1] = A[1];
		P[0] = 1.0;
		P[2] = 0.0;
		P[1] = 0.0;
		P[3] = 1.0;
	} else {
		/* Calculate eigenvalues with Mohr's circle*/
		double avg = 0.5 * (A[0] + A[1]);
		double r = sqrt(POW2(0.5 * (A[0] - A[1])) + POW2(A[2]));
		Lambda[0] = avg + r;
		Lambda[1] = avg - r;

		P[0] = Lambda[0] - A[1];
		P[2] = A[2];
		P[1] = Lambda[1] - A[1];
		P[3] = A[2];

		/* Normalize eigenvectors */
		double normalizer = sqrt(POW2(P[0]) + POW2(P[2]));
		P[0] /= normalizer;
		P[2] /= normalizer;
		normalizer = sqrt(POW2(P[1]) + POW2(P[3]));
		P[1] /= normalizer;
		P[3] /= normalizer;
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
}

double nb_mtxsym_2X2_det(const double A[3])
{
	return A[0] * A[1] - POW2(A[2]);
}

double nb_mtxsym_2X2_inverse_destructive(double A[3])
{
	double det = nb_mtxsym_2X2_det(A);
	double tmp = A[0];
	A[0] = A[1]/det;
	A[1] = tmp/det;
	A[2] = -A[2]/det;
	return det;
}
