#ifndef GRAVITY_H
#define GRAVITY_H

#include <Rcpp.h>
using namespace Rcpp;

#define NUM_TELESCOPES 4
#define NUM_BASELINES 6
#define MAX_PHASE_SHIFTS 4
#define SCI_NUM_IO_OUTPUT 23
#define SCI_NUM_PHASE_SHIFTS {4, 4, 3, 3, 4, 4}
#define SCI_IDX_PHASE_SHIFTS {0, 2, 1, 3, 0, 2, 1, 3, 0, 2, 1, -1, -1, 2, 1, 3, 0, 2, 1, 3, 0, 2, 1, 3}

#define R_IDX(I,J,M,N) (I)+(J)*(M)
#define C_IDX(I,J,M,N) (I)*(N)+(J)

List gv_const();
List gv_solvels(NumericMatrix rA, NumericVector rb, NumericVector rw);
List gv_vis2gd(ComplexVector vis);
void abcd2vis(double *ii, double *kx, double *sx, double *vv, double fx);
void gv_print_fits_err(int err);

//' Transpose a Vector/Matrix
template <typename T>
T gv_transpose(T x, int nrow, int ncol)
{
  T y(nrow*ncol);
	for (int i=0; i<nrow; i++) for (int j=0; j<ncol; j++)
		y[i*ncol+j] = x[j*nrow+i];
	return y;
}

#endif