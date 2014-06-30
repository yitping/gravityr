#include <Rcpp.h>
using namespace Rcpp;

void abcd2vis(double *ii, double *kx, double *sx,  double *vv, double fx)
{
  double v_cos_ph, v_sin_ph, v_scale;
  v_sin_ph = kx[3]*(sx[0]*ii[0]-sx[2]*ii[2]) - kx[1]*(sx[1]*ii[1]-sx[3]*ii[3]);
  v_cos_ph = kx[2]*(sx[0]*ii[0]-sx[2]*ii[2]) - kx[0]*(sx[1]*ii[1]-sx[3]*ii[3]);
  v_scale  = fx*(kx[1]*kx[2] - kx[0]*kx[3]);
  v_scale  = (v_scale < 0) ? -1/v_scale : 1/v_scale;
//  v_scale  = 1/abs(v_scale);
  vv[0] = v_cos_ph*v_scale;
  vv[1] = v_sin_ph*v_scale;
}

//' Compute complex visibility based on ABCD principle
//'
//' @param ii ABCD
//' @param kx P2VM. By default, c(0,2,2,0)
//' @param sx P2VM. By default, c(4,4,4,4)
//' @param fx Total flux in ABCD
//' 
// [[Rcpp::export]]
ComplexVector gv_abcd2vis(NumericVector ii,
  NumericVector kx = NumericVector::create(4),
  NumericVector sx = NumericVector::create(4),
  NumericVector fx = NumericVector::create(1))
{
  double vv[2] = {0,0};
  ComplexVector vis = ComplexVector(1);
  abcd2vis(ii.begin(), kx.begin(), sx.begin(), &vv[0], as<double>(fx));
  vis[0].r = vv[0];
  vis[0].i = vv[1];

  return vis;
}
