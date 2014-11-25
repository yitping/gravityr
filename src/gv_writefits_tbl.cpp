#include <Rcpp.h>
#include <fitsio.h>
#include "gRavity.h"
using namespace Rcpp;

//' Transpose a Vector/Matrix
template <typename T>
T gv_transpose(T x, int nrow, int ncol)
{
	T y(nrow*ncol);
	for (int i=0; i<nrow; i++) for (int j=0; j<ncol; j++)
		y[i*ncol+j] = x[j*nrow+i];
	return y;
}

//' FITS table writer
//' 
//' Writes a named list to a FITS file as a binary table.
//' The data is written to an extended HDU. Each member of the list is written
//' as a separate column in the table.
//'
//' Supports only numeric and integer data type for now.
//' 
// [[Rcpp::export]]
int gv_writefits_tbl(List dat, CharacterVector fits_name, CharacterVector hdu_name = "gRavity", CharacterVector units = "")
{

  fitsfile *pfits=NULL;
  int err=0;

	std::string fname = as<std::string>(fits_name[0]);
	fits_create_file(&pfits, (char *) fname.c_str(), &err);
	if (err)
	{
		gv_print_fits_err(err);
		return err;
	}
  
  // define name of each table column
  std::vector<char *> ttype;
  std::vector<std::string> ttype_str, names;
  if (dat.hasAttribute("names")) names = dat.attr("names");
  
  // define the format of each table column
  std::vector<char *> tform;
  std::vector<std::string> tform_str;

  // define the units of each table column
  std::vector<char *> tunit;
  std::vector<std::string> tunit_str;

	// determine the number of columns to be written out to the FITS file
	int nrow=0;
  std::vector<int> colsel, colsize;
	std::vector<char> coltyp;
  for (int i=0; i<dat.length(); i++)
  {
		// check data type
    switch (TYPEOF(dat[i]))
    {
      case REALSXP:
				coltyp.push_back('D');
        break;
      case INTSXP:
				coltyp.push_back('I');
        break;
      default:
				Rcout << "WARNING: unsupported data type at [" << i << "]" << std::endl;
				continue;
    }
		
    // .hasAttribute("dim") does not work
    char lastdim[10];
    if (!Rf_isNull(Rf_getAttrib(dat[i], R_DimSymbol)))
    {
      // multi-dimension array will be converted into 2D with last dimension used as
      // vector size in FITS
      IntegerVector dim = as<IntegerVector>(Rf_getAttrib(dat[i], R_DimSymbol));
			if (dim.length() > 2) Rf_error("too many dimensions");
			if (i == 0) nrow = dim[0];
			else if (nrow != dim[0]) Rf_error("inconsistent number of rows");
      sprintf(lastdim, "%d%c", dim[1], coltyp[i]);
			colsize.push_back(dim[1]);
    }
    else
		{
			if (i == 0) nrow = Rf_length(dat[i]);
			else if (nrow != Rf_length(dat[i])) Rf_error("inconsistent number of rows");
			sprintf(lastdim, "%c", coltyp[i]);
			colsize.push_back(1);
		}

		// add this column
    colsel.push_back(i);

		// define name of each table column
		ttype_str.push_back(names[i]);

		// define the format of each table column
    tform_str.push_back(lastdim);
		
		// define the units of each table column
		if (units.length() == dat.length())
			tunit_str.push_back(as<std::string>(units[i]));
  }

	ttype.resize(colsel.size(), NULL);
	tform.resize(colsel.size(), NULL);
	tunit.resize(colsel.size(), NULL);

  for (int i=0; i<colsel.size(); i++)
	{
		ttype[i] = (char *) ttype_str[i].data();
		tform[i] = (char *) tform_str[i].data();
		if (units.length() == dat.length())
			tunit[i] = (char *) tunit_str[i].data();
	}

  fits_create_tbl(pfits, BINARY_TBL, 0, colsel.size(), ttype.data(),
    tform.data(), tunit.data(), as<std::string>(hdu_name).c_str(), &err);

	for (int j=0; j<colsel.size(); j++)
	{
		switch (coltyp[j])
		{
			case 'D':
				{
					NumericVector tmp = gv_transpose<NumericVector>(as<NumericVector>(dat[colsel[j]]), nrow, colsize[j]);
					for (int i=0; i<nrow; i++)
						fits_write_col(pfits, TDOUBLE, j+1, i+1, 1, colsize[j], &(tmp[i*colsize[j]]), &err);
					break;
				}
			case 'I':
				{
					IntegerVector tmp = gv_transpose<IntegerVector>(as<IntegerVector>(dat[colsel[j]]), nrow, colsize[j]);
					for (int i=0; i<nrow; i++)
						fits_write_col(pfits, TINT, j+1, i+1, 1, colsize[j], &(tmp[i*colsize[j]]), &err);
					break;
				}
			default:
				Rf_error("unexpected column type");
		}
		if (err)
		{
			gv_print_fits_err(err);
			break;
		}
	}

	fits_close_file(pfits, &err);
  return err;
}
