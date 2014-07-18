#include <Rcpp.h>
#include "gRavity.h"
using namespace Rcpp;

//' Strip images to relevant pixels
//' 
//' gv_sci_img2pix(img, cal$idx, cal$cnr, dim(cal$idx))
// [[Rcpp::export]]
NumericVector gv_sci_img2pix(NumericMatrix img,
    NumericVector idx, IntegerVector cnr, IntegerVector dim,
		IntegerVector n_wd = 3)
{
  // TODO: get dimensions from idx directly to save 1 argument
	int nch = dim[0];
	int nio = dim[1];
	int npl = dim[2];
	int nwd = as<int>(n_wd);
	
  // TODO: return only the 
	NumericVector pix;

	if (idx.length() != nch*nio*npl) stop("size of idx matrix mismatch\n");
  if (cnr.length() != 2) stop("unexpected length of corner vector\n");

//	pix = NumericVector(Dimension(nio, nch, npl));
	pix = NumericVector(Dimension(nch, nio, npl));

	int i, j, k, l, p, cen;
	double val;
  
  // copied from gvspcImageProc.c
  /* the format of idx is determined from the txt file it was read in */
  /* number of channel x number of IO output x number of polarization states */
  for (i=0; i<nio; i++) for (j=0; j<nch; j++) for (p=0; p<npl; p++)
	{
		/* get the center of slits */
    cen = idx[p*nio*nch+C_IDX(i,j,nio,nch)]+cnr[1];
		val = 0.0;
		for (k=0; k<nwd; k++)
		{
			l = cen + k - nwd/2;
      val += img[R_IDX(j+cnr[0],l-1,img.nrow(),img.ncol())];
		} /* k */
		pix[p*nio*nch+C_IDX(i,j,nio,nch)] = val;
  } /* p,j,i */

	return pix;
}
