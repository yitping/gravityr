#include <Rcpp.h>
#include "gRavity.h"
using namespace Rcpp;

//' Compute PD
//' 
//' Copied from gvspc codes.
//' 
// [[Rcpp::export]]
ComplexMatrix gv_sci_pix2vis(NumericMatrix flux, NumericMatrix dark, List p2vm)
{
  int ps_i[] = SCI_IDX_PHASE_SHIFTS;
  int n_ps[] = SCI_NUM_PHASE_SHIFTS;
  int n_ch = flux.nrow();
  NumericVector p2vm_k = as<NumericVector>(p2vm["k"]);
  NumericVector p2vm_s = as<NumericVector>(p2vm["sx"]);
  
  int i, j, k, l, m, n, phi, missing_phi;
  double px[MAX_PHASE_SHIFTS], kx[MAX_PHASE_SHIFTS], sx[MAX_PHASE_SHIFTS], vv[2];
  double sum_ch=0.0;
  
  // Sanity check
  if (flux.ncol() != SCI_NUM_IO_OUTPUT) stop("unexpected number of io outputs\n");
  
  ComplexMatrix vis = ComplexMatrix(Dimension(n_ch, NUM_BASELINES));
  
  /* compute only 1 selected polarization */
  m = 0; /* cummulative sum of the number of phase shifts */
  for (i=0; i<NUM_BASELINES; i++)
  {
    for (j=0; j<n_ch; j++)
    {
      /* See comments in gvspc_calinfo_compute_v2pm() for details */
      l = i*MAX_PHASE_SHIFTS*n_ch + j;
      if (i > 2) l -= n_ch; /* skip the Y-junction */
      
      missing_phi = 0;
      if (n_ps[i] != MAX_PHASE_SHIFTS)
      {
        for (k=0; k<MAX_PHASE_SHIFTS; k++)
        {
          missing_phi += k;
          if ((ps_i[i*MAX_PHASE_SHIFTS+k] >= 0) &&
          (ps_i[i*MAX_PHASE_SHIFTS+k] <= MAX_PHASE_SHIFTS))
          missing_phi -= ps_i[i*MAX_PHASE_SHIFTS+k];
        }
      }
      
      sum_ch = 0.0;
      for (k=0; k<MAX_PHASE_SHIFTS; k++)
      {
        /* After this loop, px, sx and kx are arranged so that the index
        * 0,1,2,3 == A,B,C,D */
        phi = ps_i[i*MAX_PHASE_SHIFTS+k];
        if (phi < 0) continue;
        n = l+k*n_ch;
        px[phi] = flux[n] - dark[n];
        sx[phi] = p2vm_s[n];
        kx[phi] = p2vm_k[n];
        sum_ch += px[phi];
      }
      /* this corrects for baselines with less than MAX_PHASE_SHIFTS */
      if (n_ps[i] != MAX_PHASE_SHIFTS)
      {
        px[missing_phi] = sum_ch;
        kx[missing_phi]  = 0.0;
        sx[missing_phi]  = 1.0;
      }
      abcd2vis(&px[0], &kx[0], &sx[0], &vv[0], sum_ch);
      vis[i*n_ch+j].r = vv[0];
      vis[i*n_ch+j].i = vv[1];
    }
    m += n_ps[i];
  }

  return vis;
}
