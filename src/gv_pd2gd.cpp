#include <Rcpp.h>
using namespace Rcpp;

ComplexVector Conj(ComplexVector x)
{
  ComplexVector xt = clone(x);
  for (ComplexVector::iterator it = xt.begin(); it != xt.end(); it++)
    it->i = -it->i;
  return xt;
}

// [[Rcpp::export]]
ComplexVector gv_pd2gd(ComplexVector pd)
{
  int i, j, n=pd.length();
  ComplexVector cj_pd = Conj(pd);
  ComplexVector dd = ComplexVector();
  for (i=0; i<n-1; i++) for (j=i; j<n; j++)
  {
    dd.push_back(pd[i]*cj_pd[j]);
  }
  return dd;
}
