#' Computes the complex visibility of the fringes
#' 
#' Computes the complex visibilty of the SCI fringes in a given FITS file. 
#'
#' @param fits A GRAVITY FITS file
#' @param cal A list of calibration coefficients
#' @param p2vm A list of P2VM coefficients
#' @param dark Optional (temp). Dark pixels from gv_sci_img2pix().
#' @param p Optional. Polarization states to process.
#' @param hdu Optional. The hdu name that contains the SCI inteferograms.
#' @param blseq Optional. Baseline order and OPD polarity. See gv_blorder().
#' 
#' @return A 3D cube of complex visibility vs. spectral channel vs. baseline vs. number of image frames.
#'
gv_sci_fits2vis <- function (fits, p2vm, cal,
                             dark=0, p=1, hdu='imaging_data_spe',
                             blseq=NULL)
{
  img <- gv_readfits(fits, hdu)
  nd  <- length(dim(img))
  cv  <- lapply(1:ifelse(nd == 3, dim(img)[3], 1), function (i)
  {
    if (nd == 3) i2d <- img[,,i] else i2d <- img
    pix <- gv_sci_img2pix(i2d, cal$idx, cal$corner, dim(cal$idx))
    vis <- gv_sci_pix2vis(pix[,,p], dark, p2vm, blseq=blseq)
    return(vis)
  })
  cv  <- simplify2array(cv, higher=T)
}
