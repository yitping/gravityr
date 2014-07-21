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

// [[Rcpp::export]]
List gv_sci_pix2vis_unsorted(NumericVector pixels, NumericVector rdnoiz, List v2pms)
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
  
  int nrow = dim[0], n_ch = dim[1];
  NumericVector pix, rdn, w, coh, var_coh;
  NumericVector flux(n_tel), var_flux(n_tel);
//  NumericVector gd(n_bl), var_gd(n_bl);
  ComplexMatrix vis(n_bl,n_ch);
  NumericMatrix var_v_cos_pd(n_bl,n_ch), var_v_sin_pd(n_bl,n_ch);
  NumericMatrix var_v2(n_bl,n_ch), var_pd(n_bl,n_ch);
  List xx;
  int t, i, j;
  double x, y, f, e_v, e_f;
  
  for (int j=0; j<n_ch; j++)
  {
    int pix_offset = j*nrow;
    pix = NumericVector(pixels.begin()+pix_offset, pixels.begin()+pix_offset+nrow);
    rdn = NumericVector(rdnoiz.begin()+pix_offset, rdnoiz.begin()+pix_offset+nrow);
    w   = abs(pix) + rdn;
    xx  = gv_solvels(v2pms[j], pix, w);
    coh = as<NumericVector>(xx["x"]); var_coh = as<NumericVector>(xx["var_x"]);
    if (j == 0)
    {
      flux      = NumericVector(coh.begin(), coh.begin()+4);
      var_flux  = NumericVector(var_coh.begin(), var_coh.begin()+4);
    }
    else
    {
      flux     += NumericVector(coh.begin(), coh.begin()+4);
      var_flux += NumericVector(var_coh.begin(), var_coh.begin()+4);
    }
		// copied from cppgvspc but modified [] -> () for NumericMatrix
    for (i=0; i<n_bl; i++)
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
  
//  for (i=0; i<n_bl; i++)
//  {
//    xx = gv_vis2gd(vis.row(i));
//    gd[i] = as<double>(xx["val"]); var_gd[i] = as<double>(xx["var"]);
//  }
  
  return List::create(
    Named("flux")=flux,
    Named("var_flux")=var_flux,
//    Named("gd")=gd,
//    Named("var_gd")=var_gd,
    Named("vis")=vis,
    Named("var_v2")=var_v2,
    Named("var_pd")=var_pd
    );
}

