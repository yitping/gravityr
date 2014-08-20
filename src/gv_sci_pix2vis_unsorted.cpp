#include <Rcpp.h>
#include "gRavity.h"
using namespace Rcpp;

// copied from cppgvspc
double get_var_v2(double x, double y, double ex, double ey)
{
  double e;
	e  = 4*x*x*ex;
	e += 4*y*y*ey;
	return e;
}

// copied from cppgvspc
double get_var_pd(double x, double y, double ex, double ey)
{
	double e;
	e  = x*x*ex;
	e += y*y*ey;
	e *= (x*x + y*y)*(x*x + y*y);
	return e;
}

//' Compute complex visibility of SCI fringes
//' 
//' Copied from gvspc codes.
//' @param pixels Matrix of pixel values.
//' @param v2pms List of v2pms.
//' @param rdnoiz Optional. Matrix of variance of each pixels (obtained from dark frames). 
//' @param bad Optional. Vector indicating bad spectral channels.
//' 
//' @examples \dontrun{
//'  img  <- gv_readfits(fits$dark, hdu='imaging_data_spe')
//'  dark_img  <- apply(img[,,1:dim(img)[3]], c(1,2), mean)
//'  dark_pixs <- lapply(1:dim(img)[3], function (i) gv_sci_img2pix(img[,,i], cal$idx, cal$corner))
//'  dark_pixs <- simplify2array(dark_pixs, higher=T)
//'  rdnoiz <- apply(dark_pixs, c(1,2,3), var)
//'  img <- gv_readfits(fits$flux, hdu='imaging_data_spe')
//'  flux_pixs <- lapply(1:dim(img)[3], function (i) gv_sci_img2pix(img[,,i]-dark_img, cal$idx, cal$corner))
//'  bad <- which(sapply(v2pms, function (v) any(is.na(v))))
//'  coh <- gv_sci_pix2vis_unsorted(flux_pixs[[1]][,,1], v2pms, rdnoiz, bad)
//' }
//' 
// [[Rcpp::export]]
List gv_sci_pix2vis_unsorted(NumericVector pixels, List v2pms,
  NumericVector rdnoiz = NumericVector(), NumericVector bad = NumericVector())
{
  IntegerVector dim = pixels.attr("dim");
  if ((dim.length() == 2) && (dim[0] != MAX_PHASE_SHIFTS*NUM_BASELINES))
  {
    REprintf("unexpected pixel dimension");
    return List::create();
  }
  if (dim[1] != v2pms.length())
  {
    REprintf("pixel dimension and v2pms length mismatch");
    return List::create();
  }
  
  // for codes copied from cppgvspc
  int n_tel=4, n_bl=6;
  NumericMatrix tels(2,6);
  // 2, 0, 0, 1, 1, 0 -> last bit slightly different from gv_telmat()
  // 3, 3, 2, 3, 2, 1 -> last bit slightly different from gv_telmat()
  tels(0,0)=2; tels(0,1)=0; tels(0,2)=0; tels(0,3)=1; tels(0,4)=1; tels(0,5)=0;
  tels(1,0)=3; tels(1,1)=3; tels(1,2)=2; tels(1,3)=3; tels(1,4)=2; tels(1,5)=1;
  
  int n_row = dim[0], n_ch = dim[1];
	int n_col = as<NumericVector>(v2pms[0]).length()/n_row;
  NumericVector pix, rdn, w, coh, var_coh;
  NumericVector flux(n_tel), var_flux(n_tel);
  ComplexMatrix vis(n_bl,n_ch);
  NumericMatrix var_v_cos_pd(n_bl,n_ch), var_v_sin_pd(n_bl,n_ch);
  NumericMatrix var_v2(n_bl,n_ch), var_pd(n_bl,n_ch);
  List xx;
  int t, i, j;
  double x, y, f, e_v, e_f;
  
  // bad channels, NB: given index starts with 1 (R style)
  NumericVector badch = NumericVector(n_ch,0.0);
  if (bad.length() != 0) for (j=0; j<bad.length(); j++) badch[bad[j]-1] = 1;
  
  for (j=0; j<n_ch; j++)
  {
    int pix_offset = j*n_row;
    pix = NumericVector(pixels.begin()+pix_offset, pixels.begin()+pix_offset+n_row);
    w   = abs(pix);
    if (rdnoiz.length() != 0)
    {
      rdn = NumericVector(rdnoiz.begin()+pix_offset, rdnoiz.begin()+pix_offset+n_row);
      w  += rdn;
    }
    // if this is a bad channel, don't bother solving the LS
    if (badch[j] != 1)
    {
      xx  = gv_solvels(v2pms[j], pix, w);
      coh = as<NumericVector>(xx["x"]);
			var_coh = as<NumericVector>(xx["var_x"]);
			for (i=0; i<n_col; i++)
			{
				// should not expect any print out, else more coding to take care of
				// this error.
				if (isnan(coh[i])) Rcout << "j = " << j << " is bad too (i=" << i << ")" << std::endl;
			}
    }
    else
		{
			coh = NumericVector(n_col, NAN);
			var_coh = NumericVector(n_col, NAN);
		}
		// compute flux
    if (j == 0)
    {
      flux      = (badch[j] == 1) ? NumericVector(n_tel,0.0) : NumericVector(coh.begin(), coh.begin()+4);
      var_flux  = (badch[j] == 1) ? NumericVector(n_tel,0.0) : NumericVector(var_coh.begin(), var_coh.begin()+4);
    }
    else
    {
      flux     += (badch[j] == 1) ? NumericVector(n_tel,0.0) : NumericVector(coh.begin(), coh.begin()+4);
      var_flux += (badch[j] == 1) ? NumericVector(n_tel,0.0) : NumericVector(var_coh.begin(), var_coh.begin()+4);
    }
		// compute complex visibility
		// copied from cppgvspc but modified [] -> () for NumericMatrix
    for (i=0; i<n_bl; i++)
		{
			if (badch[j] == 1)
			{
				vis(i,j).r = NAN;
				vis(i,j).i = NAN;
				var_v2(i,j) = NAN;
				var_pd(i,j) = NAN;
			}
			else
			{
				f = coh[tels(0,i)] + coh[tels(1,i)];
				x = coh[n_tel+0*n_bl+i];
				y = coh[n_tel+1*n_bl+i];
				vis(i,j).r = x/f;
				vis(i,j).i = y/f;
				e_f  = var_coh[tels(0,i)] + var_coh[tels(1,i)];
				e_f /= f*f;
				e_v  = var_coh[n_tel+0*n_bl+i];
				e_v /= x*x;
				var_v_cos_pd(i,j)  = vis(i,j).r*vis(i,j).r;
				var_v_cos_pd(i,j) *= e_v + e_f;
				e_v  = var_coh[n_tel+1*n_bl+i];
				e_v /= y*y;
				var_v_sin_pd(i,j)  = vis(i,j).i*vis(i,j).i;
				var_v_sin_pd(i,j) *= e_v + e_f;
				// copied from another section of cppgvspc
				var_v2(i,j) = get_var_v2(vis(i,j).r, vis(i,j).i, var_v_cos_pd(i,j), var_v_sin_pd(i,j));
				var_pd(i,j) = get_var_pd(vis(i,j).r, vis(i,j).i, var_v_cos_pd(i,j), var_v_sin_pd(i,j));
			}
		}
  }
  
  return List::create(
    Named("flux")=flux,
    Named("var_flux")=var_flux,
    Named("vis")=vis,
    Named("var_v2")=var_v2,
    Named("var_pd")=var_pd
    );
}

