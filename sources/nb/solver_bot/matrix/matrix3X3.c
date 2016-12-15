#include "nb/solver_bot/matrix/matrix3X3.h"

double nb_matrix_3X3_det(double *A)
{
	return 
		A[0] * A[4] * A[8] + 
		A[3] * A[7] * A[2] +
		A[6] * A[1] * A[5] -
		A[2] * A[4] * A[6] -
		A[5] * A[7] * A[0] -
		A[8] * A[1] * A[3];
}

double nb_matrix_3X3_inverse_destructive(double *A)
{
	double det = nb_matrix_3X3_det(A);
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	a11 = A[0];
	a12 = A[1];
	a13 = A[2];
	a21 = A[3];
	a22 = A[4];
	a23 = A[5];
	a31 = A[6];
	a32 = A[7];
	a33 = A[8];
	A[0] = (a33*a22-a32*a23)/det;
	A[1] = -(a33*a12-a32*a13)/det;
	A[2] = (a23*a12-a22*a13)/det;
	A[3] = -(a33*a21-a31*a23)/det;
	A[4] = (a33*a11-a31*a13)/det;
	A[5] = -(a23*a11-a21*a13)/det;
	A[6] = (a32*a21-a31*a22)/det;
	A[7] = -(a32*a11*a31*a12)/det;
	A[8] = (a22*a11-a21*a12)/det;
	return det;
}
