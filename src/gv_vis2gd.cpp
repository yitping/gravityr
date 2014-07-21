#include <Rcpp.h>
#include <complex>
using namespace Rcpp;

//' Compute group delay from complex visibilities
//' 
//' Copied from cppgvspc
//' 
// [[Rcpp::export]]
List gv_vis2gd(ComplexVector vis)
{
  int n_ch=vis.length();
  double gd, var_gd;
  
  int i, j;
	std::complex<double> d, m;
	std::vector<std::complex<double> > dd(n_ch-1);
	double e, s;
  
	// compute average of phasor difference
	for (j=0; j<(n_ch-1); j++)
	{
	  d  = std::complex<double>(vis(j).r,    vis(j).i);
	  d *= std::complex<double>(vis(j+1).r, -vis(j+1).i);
	  m  = (j == 0) ? d : m+d;
	  dd[j] = d;
	}
	gd = arg(m)*(n_ch-1);
	// compute the deviation of
	m = conj(m);
	for (j=0; j<(n_ch-1); j++)
	{
	  e  = arg(dd[j]*m);
	  e *= e;
	  s  = (j == 0) ? e : s+e;
	}
	var_gd = s*(n_ch-1)/(n_ch-2);
  
	return List::create(Named("val")=gd, Named("var")=var_gd);
  
}
