#include <Rcpp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include "gRavity.h"
using namespace Rcpp;

#define NUMERIC_SINGULAR_VALUE 1e-12

void sv_inverse(gsl_matrix *U, gsl_matrix *V, gsl_vector *s, gsl_matrix *Z, int m, int n)
{
	if ((U == NULL) || (V == NULL) || (s == NULL))
		stop("unexpected null matrices/vectors\n");

	double si;
	gsl_matrix *S = gsl_matrix_calloc(n, n);
	gsl_matrix *X = gsl_matrix_alloc(n, n);

	// diag(S) <- s
	for (int i=0; i<n; i++)
	{
		si = gsl_vector_get(s,i);
		gsl_matrix_set(S, i, i, (si == 0) ? 0 : 1/si);
	}
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, S, 0.0, X);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, X, U, 0.0, Z);
	

	gsl_matrix_free(X);
	gsl_matrix_free(S);
}

void sv_threshold(gsl_vector *s, int n)
{
	for (int i=0; i<n; i++)
		if ((gsl_vector_get(s,i) < NUMERIC_SINGULAR_VALUE) &&
				(gsl_vector_get(s,i) > -NUMERIC_SINGULAR_VALUE))
			gsl_vector_set(s,i,0);
}


//' Solve for x in an Ax=b linear system
//'
//' The solution for x is inv(A)%*%b where inv(A) is computed by means of the
//' singular value decomposition method.
//'
// [[Rcpp::export]]
List gv_solvels(NumericMatrix rA, NumericVector rb, NumericVector rw = NumericVector())
{
  if (rA.nrow() != rb.length())
  {
    REprintf("A and b sizes mismatch\n");
    return R_NilValue;
  }

	int nrow = rA.nrow();
	int ncol = rA.ncol();
	int i, j;
	
	// setup the main structure, i.e. Ax = b
	gsl_matrix *A = gsl_matrix_alloc(nrow, ncol);
	gsl_vector *b = gsl_vector_alloc(nrow);
	gsl_matrix *W = (rw.length() == 0) ? NULL : gsl_matrix_calloc(nrow, nrow);
	NumericVector x(ncol), var_x(ncol);
  
	for (i=0; i<nrow; i++)
	{
		gsl_vector_set(b, i, rb[i]);
		if (W != NULL) gsl_matrix_set(W, i, i, (rw[i] == 0) ? 0 : 1/rw[i]);
		for (j=0; j<ncol; j++)
		{
			gsl_matrix_set(A, i, j, rA(i,j));
		}
	}


	// copied from cppgvspc
	gsl_matrix *C  = gsl_matrix_alloc(ncol, ncol);
	gsl_matrix *U  = gsl_matrix_alloc(ncol, ncol);
	gsl_matrix *V  = gsl_matrix_alloc(ncol, ncol);
	gsl_vector *s  = gsl_vector_alloc(ncol);
	gsl_vector *a  = gsl_vector_alloc(ncol);
	gsl_vector *r  = gsl_vector_alloc(nrow);
	gsl_vector *Mb = gsl_vector_alloc(ncol);
	gsl_matrix *WA = NULL;
	gsl_vector *Wb = NULL;
	
	if (W != NULL)
	{
		WA = gsl_matrix_alloc(nrow, ncol);
		Wb = gsl_vector_alloc(nrow);
		// first, compute the inverse covariance matrix of x
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W, A, 0.0, WA);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, WA, 0.0, C);
		// then, compute t(A)*W*b -> Mb
		gsl_blas_dgemv(CblasNoTrans, 1.0, W, b, 0.0, Wb);
		gsl_blas_dgemv(CblasTrans, 1.0, A, Wb, 0.0, Mb);
	}
	else
	{
		// first, compute the inverse covariance matrix of x
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, C);
		// then, compute t(A)*b -> Mb
		gsl_blas_dgemv(CblasTrans, 1.0, A, b, 0.0, Mb);
	}
	
	// then, decompose C
	gsl_matrix_memcpy(U, C);
	gsl_linalg_SV_decomp(U, V, s, a);

	// apply threshold for singular values
	sv_threshold(s, ncol);
	
	// solve it!
	gsl_linalg_SV_solve(U, V, s, Mb, a);

	// save x
	for (int i=0; i<ncol; i++) x[i] = gsl_vector_get(a, i);

	// compute covariance of x by inverting C
	// C <- Cov(x), var_x <- sigma(residual)^2*diag(C)
	sv_inverse(U, V, s, C, ncol, ncol);
	
	if (W != NULL)
	{
		for (int i=0; i<ncol; i++) var_x[i] = gsl_matrix_get(C, i, i);
	}
	else
	{
		// compute residuals
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, a, 0.0, r);
		std::vector<double> resd(nrow,0);
		for (int i=0; i<nrow; i++)
			resd[i] = gsl_vector_get(r, i) - gsl_vector_get(b, i);
		double var_resd = gsl_stats_variance(resd.data(), 1, nrow);
		for (int i=0; i<ncol; i++) var_x[i] = var_resd*gsl_matrix_get(C, i, i);
	}

	
	if (W != NULL)
	{
		gsl_vector_free(Wb);
		gsl_matrix_free(WA);
		gsl_matrix_free(W);
	}
	gsl_vector_free(Mb);
	gsl_vector_free(r);
	gsl_vector_free(a);
	gsl_vector_free(s);
	gsl_matrix_free(V);
	gsl_matrix_free(U);
	gsl_matrix_free(C);
	
	gsl_vector_free(b);
	gsl_matrix_free(A);

  return List::create(
			Named("x")=x,
			Named("var_x")=var_x);
}
