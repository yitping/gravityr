#include <Rcpp.h>
#include "gRavity.h"
using namespace Rcpp;

//' Constants of the spectrometers
//' 
//' Constants which are used in several other functions.
//' 
//' @return A named list of constants
//'
// [[Rcpp::export]]
List gv_const()
{
  int default_sci_n_ps[] = SCI_NUM_PHASE_SHIFTS;
  int default_sci_ps_i[] = SCI_IDX_PHASE_SHIFTS;
  std::vector<int> sci_n_ps(std::begin(default_sci_n_ps), std::end(default_sci_n_ps));
  std::vector<int> sci_ps_i(std::begin(default_sci_ps_i), std::end(default_sci_ps_i));
  
  return List::create(
    Named("ft_io_out")=24,
    Named("ft_telarr")=NumericVector::create(2,3,1,4),
    Named("sci_io_out")=SCI_NUM_IO_OUTPUT,
    Named("sci_telarr")=NumericVector::create(2,1,3,4),
    Named("sci_io_phases")=wrap(sci_n_ps),
    Named("sci_io_phase")=wrap(sci_ps_i)
    );
}
