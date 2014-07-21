#include <Rcpp.h>
#include "gRavity.h"
using namespace Rcpp;

// copied from gvspcPix class
long convert_img_indices(int l, int j, int p, int isYjunc, int n_ch, int n_pl)
{
	int n_ph = MAX_PHASE_SHIFTS, n_bl = NUM_BASELINES;
  // l - io output, j - spectral channel, p - polarization
	int i, k; // k - phase, i - baseline
	if (isYjunc == 1) // at the Y-junction
	{
		return p*n_ch + j;
	}
	else
	{
		i = (l < 11) ? l/4 : ((l > 11) ? (l+1)/4 : -1); // translate l-th IO to i-th baseline
		k = (l < 11) ? l%4 : ((l > 11) ? (l+1)%4 : -1); // translate l-th IO to k-th phase shift
		k = (k == 1) ? 2 : ((k == 2) ? 1 : k); // swap phase C and B
		return p*n_ch*n_bl*n_ph + j*n_bl*n_ph + i*n_ph + k;
	}
}

//' Strip images to relevant pixels
//' 
//' gv_sci_img2pix(img, cal$idx, cal$cnr)
//' 
// [[Rcpp::export]]
NumericVector gv_sci_img2pix(NumericMatrix img,
    NumericVector idx, IntegerVector cnr, IntegerVector isYjunc = 0, IntegerVector n_wd = 3)
{
//  if (idx.length() < 10)
//	{
//		REprintf("unexpected size of idx array\n");
//		return wrap(0);
//	}
//  if (cnr.length() != 2)
//	{
//		REprintf("unexpected length of corner vector\n");
//		return wrap(0);
//	}
//  if (img.length() < 10)
//  {
//    REprintf("unexpected size of image");
//    return wrap(0);
//  }
//  if (dim.length() != 3)
//	{
//		REprintf("unexpected length of idx vector\n");
//		return wrap(0);
//	}
  
  IntegerVector dim = idx.attr("dim");
	int n_ph = MAX_PHASE_SHIFTS, n_bl = NUM_BASELINES;
	int n_ch = dim[0];
	int n_io = dim[1];
	int n_pl = dim[2];
	int nwd = as<int>(n_wd);
	int isY = as<int>(isYjunc);
	
	std::vector<double> v(n_ph*n_bl*n_ch*n_pl, 0);
	std::vector<double> y(n_ch*n_pl, 0);

	if ((idx.length() != n_ch*n_io*n_pl) && (n_io != SCI_NUM_IO_OUTPUT))
	{
		REprintf("size of idx matrix mismatch\n");
		return wrap(0);
	}

	int i, j, k, l, p, cen;
	double sum;
  
  // copied from gvspcImageProc.c
  /* the format of idx is determined from the txt file it was read in */
  /* number of channel x number of IO output x number of polarization states */
  for (i=0; i<n_io; i++) for (j=0; j<n_ch; j++) for (p=0; p<n_pl; p++)
	{
		/* get the center of slits */
    cen = idx[p*n_io*n_ch+C_IDX(i,j,n_io,n_ch)] + cnr[1];
		sum = 0.0;
		for (k=0; k<nwd; k++)
		{
			l = cen + k - nwd/2;
      sum += img[R_IDX(j+cnr[0],l-1,img.nrow(),img.ncol())];
		} /* k */
		//pix[p*nio*nch+C_IDX(i,j,nio,nch)] = sum;
		if (i != 11)
			v[convert_img_indices(i, j, p, 0, n_ch, n_pl)] = sum;
		else
			y[convert_img_indices(i, j, p, 1, n_ch, n_pl)] = sum;
  } /* p,j,i */
  
	NumericVector pix;
  if (isY == 0)
  {
    pix = NumericVector(v.begin(), v.end());
    pix.attr("dim") = IntegerVector::create(n_ph*n_bl, n_ch, n_pl);
  }
  else
  {
    pix = NumericVector(y.begin(), y.end());
    pix.attr("dim") = IntegerVector::create(1, n_ch, n_pl);
  }

	return pix;
}
