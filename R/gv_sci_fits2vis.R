#' Computes the complex visibility of the fringes
#' 
#' Computes the complex visibilty of the SCI fringes in a given FITS file. 
#'
#' @param fits A GRAVITY FITS file
#' @param cal A list of calibration coefficients
#' @param v2pms A list of V2PM coefficients
#' @param dark Dark pixels and read noise from gv_sci_Nfits2pix().
#' @param p Optional. Polarization states to process.
#' @param hdu Optional. The hdu name that contains the SCI inteferograms.
#' @param blseq Optional. Baseline order and OPD polarity. See gv_blorder().
#' 
#' @return A 3D cube of complex visibility vs. spectral channel vs. baseline vs. number of image frames.
#'
gv_sci_fits2vis <- function (fits, dark, v2pms, cal,
                             p=1, hdu='imaging_data_spe',
                             blseq=NULL)
{
  # TODO: check hdu is correct
  image <- gv_readfits(fits, hdu)
  nd    <- length(dim(image))
  pixs  <- lapply(1:ifelse(nd == 3, dim(image)[3], 1), function (i)
  {
    if (nd == 3) { img <- image[,,i] } else { img <- image }
    gv_sci_img2pix(img, cal$idx, cal$corner)-dark$mu
  })
  n_ch <- dim(pixs[[1]])[2]
  n_pl <- dim(pixs[[1]])[3]
  cv  <- lapply(1:length(pixs), function (i)
  {
    vis <- gv_sci_pix2vis(pixs[[i]][,,p], dark$var[,,p], v2pms[(p-1)*n_ch+(1:n_ch)], blseq=blseq)
  })
}
