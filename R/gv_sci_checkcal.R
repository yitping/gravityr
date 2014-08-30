#' Check SCI detector calibration text file
#' 
#' Cross-check the indices in the calibration file, which are shared between
#' R and C/C++ functions, with an actual FITS file.
#' 
#' 2014-05-29, Yitping - Initial rev.
#' 2014-06-03, Yitping - To support the GRAVITY merged FITS format
#'
#' @param file_cal The calibration file
#' @param file_fits A FITS file containing images from the SCI detector
#' @param file_eps Optional. An eps file to export plot. By default, no.
#' @param rev Optional. For development purpose.
#' @param hdu Optional. The name of the HDU containing the images. For primary HDU, set it to an empty string.
#' 
gv_sci_checkcal <- function (file_cal, file_fits, file_dark=NULL, file_eps=NULL,
		       rev=3, hdu='imaging_data_spe')
{
  imDat <- gv_readfits(file_fits, hdu)
  if (length(dim(imDat)) > 2) imDat <- apply(imDat, c(1,2), mean)
  if (!is.null(file_dark))
  {
    imDrk <- gv_readfits(file_dark, hdu)
    if (length(dim(imDrk)) > 2) imDrk <- apply(imDrk, c(1,2), mean)
    imDat <- imDat-imDrk
  }
  cal <- gv_sci_readcal(file_cal)
  if (rev == 1)
  {
    i <- as.matrix(read.csv(file_cal, header=F))
    corner <- i[1:2]
    nch <- i[3]
    nsl <- i[4]
    i <- i[-(1:4)]
    dim(i) <- c(nch, nsl)
    i <- apply(i, 2, mean)

  } else if (rev == 2) {
    corner <- as.numeric(read.csv(file_cal, header=F, skip=3, nrow=1))
    i <- as.numeric(read.csv(file_cal, header=F, skip=5, nrow=1))
    nch <- i[1]
    nsl <- i[2]
    i <- as.numeric(read.csv(file_cal, header=F, skip=7, nrow=1))
    # adjust corner if the index of the 1st spectral channel is not 1
    corner[1] <- corner[1] + i[1]
    i <- as.matrix(read.csv(file_cal, header=F, skip=9, nrow=nsl))
    i <- apply(i, 1, mean)

  } else if (rev == 3) {
    corner <- cal$corner
    dim_idx <- dim(cal$idx)
    nch <- dim_idx[1]
    nsl <- prod(dim_idx[2:3])
    i <- t(apply(cal$idx, c(2,3), mean))
    dim(i) <- NULL
  }
  # corner[1] - defines the index offset of the 1st spectral channel
  # corner[2] - defines the index offset of the 1st IO slit
  xa <- corner[1]+1-50
  xb <- corner[1]+nch+49
  x  <- ifelse(xa<1,1,xa):ifelse(xb>dim(imDat)[1],dim(imDat)[1],xb)
  y  <- (corner[2]+1):dim(imDat)[2]
  sl <- apply(imDat[x,y], 2, mean);
  ch <- apply(imDat[x,y], 1, mean);

  if (!is.null(file_eps))
  {
    setEPS(width=20, height=8)
    postscript(file_eps, horizontal=F)
  }
  par(mfrow=c(1,2), mar=c(4,2,1,1), oma=c(1,3,0,0))

  plot(x, ch, type='l', xlab='Spectral channels')
  abline(v=corner[1]+1,   lty=3, col='red')
  abline(v=corner[1]+nch, lty=3, col='red')
  if (rev == 3)
  {
    y_axis_max <- max(axTicks(2))
    text(corner[1]+1,   y_axis_max, format(cal$s_um[1], digits=4))
    text(corner[1]+nch, y_axis_max, format(cal$s_um[nch], digits=4))
  }

  plot(y, sl, type='l', xlab='IO slits')
  for (j in 1:nsl)
  {
    abline(v=corner[2]+i[j], lty=3, col='red')
  }

  title(ylab='Flux (ADU)', line=1, outer=T)
  if (length(file_eps) != 0) dev.off()
  return (list(corner=corner, idx=i, nch=nch, nsl=nsl))
}
