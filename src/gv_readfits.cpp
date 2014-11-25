#include <Rcpp.h>
#include <fitsio.h>
#include "gRavity.h"
using namespace Rcpp;

#define MAX_IMG_DIM 3

void gv_print_fits_err(int err)
{
  char errmsg[256];
  fits_get_errstatus(err, errmsg);
  REprintf("WARNING: [%d] %s\n", err, errmsg);
}

int gv_readfits_tab(fitsfile *pfits, int colnum, long nrows, long ncols, double *vec)
{
  long i;
  int  j, anynull=0, err=0;
  double nullval=0.0, data[ncols];
  
#ifdef GV_DEBUG
  Rcout << "column id = " << colnum << std::endl;
  Rcout << "m = " << nrows << std::endl;
  Rcout << "n = " << ncols << std::endl;
  Rcout << "length(data) = " << sizeof(data)/sizeof(double) << std::endl;
#endif
  
  // read a FITS column the old school style
  for (i=0; i<nrows; i++)
  {
    fits_read_col(pfits, TDOUBLE, colnum, i+1, 1, ncols,
				&nullval, data, &anynull, &err);
    for (j=0; j<ncols; j++) vec[j*nrows+i] = data[j];
  }
  return err;
}

int gv_readfits_img(fitsfile *pfits, long n_elem, double *vec)
{
  int i, anynull=0, err=0;
  double nullval=0;
 
  // read a FITS image the old school style
  fits_read_img(pfits, TDOUBLE, 1, n_elem, &nullval, vec, &anynull, &err);
  
  return err;
}

//' FITS file reader
//' 
//' Reads a specific image HDU or extended binary table of a FITS file.
//' This reader function is tailor made for GRAVITY FITS format but can be used
//' for generic FITS files. However, not all data types are supported.
//' 
//' Extended ASCII table and complex data types are currently not supported.
//' 
//' @param fits_name FITS file. If a vector is given, only the first is read.
//' @param hdu_name Name of the HDU to be read. Case-insensitive. If not provided, the Primary HDU is read.
//' @param col_name Name of the column of an extended HDU to be read. This is ignored if the HDU is an image.
//' 
//' @return A numeric array.
//' 
// [[Rcpp::export]]
NumericVector gv_readfits(CharacterVector fits_name,
  CharacterVector hdu_name = "", CharacterVector col_name = "")
{
  fitsfile *pfits=NULL;
  int i, err=0, hdutype=0, hdunum=0, coltype=0, colnum=0, naxes=0;
  long nrows=0, ncols=0, axis[MAX_IMG_DIM], n_elem=1;
  NumericVector d;
  
#ifdef GV_DEBUG
  Rcout << "size of fits_name is " << fits_name.size() << std::endl;
  Rcout << fits_name[0] << std::endl;
#endif

  if (fits_name.size() > 1)
  {
    REprintf("WARNING: only the first file will be read.\n");
  }
  
  std::string fname   = as<std::string>(fits_name[0]);
  std::string hduname = as<std::string>(hdu_name[0]);
  std::string colname = as<std::string>(col_name[0]);
  
  if (fits_open_file(&pfits, (char *) fname.c_str(), READONLY, &err) != 0)
  {
    gv_print_fits_err(err);
    return wrap(0);
  }
  
	// By default, HDU to read is the Primary HDU
	// otherwise move to the extended HDU given by the user in hduname
  if (hduname.length() != 0)
  {
    if (fits_movnam_hdu(pfits, ANY_HDU, (char *) hduname.c_str(), 0, &err) != 0)
    {
      gv_print_fits_err(err);
      return wrap(0);
    }
    fits_get_hdu_num(pfits, &hdunum);
  }

  // Check HDU type and handle them accordingly
  if (fits_get_hdu_type(pfits, &hdutype, &err) != 0)
  {
    gv_print_fits_err(err);
    return wrap(0);
  }
  switch (hdutype)
  {
    case IMAGE_HDU:
      fits_get_img_dim(pfits, &naxes, &err);
      if (naxes > MAX_IMG_DIM) naxes = MAX_IMG_DIM;
      fits_get_img_size(pfits, MAX_IMG_DIM, &axis[0], &err);
			for (i=0; i<naxes; i++) n_elem *= axis[i];
#ifdef GV_DEBUG
      Rcout << "IMAGE_HDU" << std::endl;
      Rcout << "hdunum: " << hdunum << std::endl;
      Rcout << "Number of axes: " << naxes << std::endl;
      for (i=0; i<naxes; i++)
				Rcout << "Size of axis[" << i << "]: " << axis[i] << std::endl;
			Rcout << "Number of pixels: " << n_elem << std::endl;
#endif
			// assume we don't have to deal with an image hypercube!
      switch (naxes)
      {
        case 3:
          d = NumericVector(Dimension(axis[0], axis[1], axis[2]));
          break;
        case 2:
          d = NumericVector(Dimension(axis[0], axis[1]));
          break;
        default:
          d = NumericVector(Dimension(axis[0]));
          break;
      }
#ifdef GV_DEBUG
      Rcout << "Length of d: " << d.length() << std::endl;
      Rcout << "Size of d: " << d.size() << std::endl;
#endif
      if ((err=gv_readfits_img(pfits, n_elem, d.begin())) != 0)
        gv_print_fits_err(err);
      break;
      // end of IMAGE_HDU

    case BINARY_TBL:
      if (colname.length() == 0)
      {
        REprintf("column name must be specified!\n");
        return wrap(0);
      }
      fits_get_num_rows(pfits, &nrows, &err);
      fits_get_colnum(pfits, CASEINSEN, (char *) colname.c_str(), &colnum, &err);
      fits_get_coltype(pfits, colnum, &coltype, &ncols, NULL, &err);
      d = NumericVector(Dimension(nrows, ncols));
#ifdef GV_DEBUG
      Rcout << "BINARY_TBL" << std::endl;
      Rcout << "hdunum: " << hdunum << std::endl;
      Rcout << "Number of rows: " << nrows << std::endl;
      Rcout << "Number of cols: " << ncols << std::endl;
      Rcout << "Index of *" << colname << "*: " << colnum << std::endl;
      Rcout << "Length of d: " << d.length() << std::endl;
      Rcout << "Size of d: " << d.size() << std::endl;
#endif
      if ((err=gv_readfits_tab(pfits, colnum, nrows, ncols, d.begin())) != 0)
        gv_print_fits_err(err);
      break;
      // end of BINARY_TBL

    case ASCII_TBL:
      // should never reach here!
      Rcout << "ASCII_TBL" << std::endl;
    default:
      REprintf("unsupported foreign table type");
      return wrap(0);
  }
  
  if (pfits != NULL)
    if (fits_close_file(pfits, &err) != 0) gv_print_fits_err(err);
  
  return d;
}
