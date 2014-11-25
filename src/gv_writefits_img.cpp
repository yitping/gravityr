#include <Rcpp.h>
#include <fitsio.h>
#include "gRavity.h"
using namespace Rcpp;

//' FITS image writer
//' 
//' Writes a vector, matrix or 3D array to a FITS file as an image.
//' The data is written to the primary HDU.
//' 
// [[Rcpp::export]]
int gv_writefits_img(NumericVector img, CharacterVector fits_name, CharacterVector hdu_name = "")
{
  IntegerVector dim;
  if (!img.hasAttribute("dim"))
  {
    REprintf("ERROR: image has not been dimensioned.\n");
    return 1;
  }
  dim = img.attr("dim");
  if (dim.length() > 3)
	{
		REprintf("ERROR: dimension of more than 3 unsupported.\n");
		return 1;
	}

  fitsfile *pfits=NULL;
	int err=0;

	std::string fname = as<std::string>(fits_name[0]);
	fits_create_file(&pfits, (char *) fname.c_str(), &err);
	if (err)
	{
		gv_print_fits_err(err);
		return err;
	}

#ifdef GV_DEBUG
  Rcout << "Number of dim: " << dim.length() << std::endl;
	for (int i=0; i<dim.length(); i++)
	{
		Rcout << "Dim[" << i << "]: " << dim[i] << std::endl;
	}
	Rcout << "Number of elements: " << img.length() << std::endl;
  double *p = &(*img.begin());
  for (int i=0; i<img.length(); i++)
  {
    Rcout << "*(p+" << i << ") = " << *(p+i) << std::endl;
  }
#endif
  
  long longdim[3], startpix[3] = {1,1,1}; // default start
  for (int i=0; i<dim.length(); i++) longdim[i] = (long) dim[i];
  
	// start writing to file
	fits_create_img(pfits, DOUBLE_IMG, dim.length(), longdim, &err);
	fits_write_pix(pfits, TDOUBLE, startpix, img.length(), &(*img.begin()), &err);

	fits_close_file(pfits, &err);
  return err;
}

