#' Read multiple fits
#' 
#' @param fits FITS files
#' @param cal List of calibration matrices
#' @param hdu Name of the HDU table to be read
#' 
#' @export
#' 
gv_sci_Nfits2pix <- function (fits, cal, hdu='')
{
  n_fits <- length(fits)
  n_frames <- 0
  pix_mu <- pix_var <- 0
  for (i in 1:n_fits)
  {
    img <- gv_readfits(fits[i], hdu)
    if (length(dim(img)) == 3)
    {
      pixs <- lapply(1:dim(img)[3], function (i) gv_sci_img2pix(img[,,i], cal$idx, cal$corner))
      pixs <- simplify2array(pixs, higher=T)
      pix_mu  <- pix_mu  + apply(pixs, c(1,2,3), sum)
      pix_var <- pix_var + apply(pixs^2, c(1,2,3), sum)
      n_frames<- n_frames + dim(img)[3]
      
    } else {
      pix <- gv_sci_img2pix(img, cal$idx, cal$corner)
      pix_mu  <- pix_mu+pix
      pix_var <- pix_var+pix^2
      n_frames<- n_frames + 1
    }
  }
  pix_mu  <- pix_mu/n_frames
  pix_var <- pix_var/(n_frames-1)-pix_mu^2*n_frames/(n_frames-1)
  return(list(mu=pix_mu, var=pix_var))
}
