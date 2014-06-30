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
#' @param ext Optional. File extension of file_cal.
#' @param rev Optional. For development purpose.
#' @param hdu Optional. The name of the HDU containing the images. For primary HDU, set it to an empty string.
#' @param frnum Optional. The frame of a 3D cube of images to check the calibration with.
#' 
gv_sci_checkcal <- function (file_cal, file_fits, file_eps='', ext='.csv',
		       rev=3, hdu='imaging_data_spe', frnum=1)
{
  imDat <- gv_readfits(file_fits, hdu)
  cal <- gv_sci_readcal(file_cal)
  if (length(dim(imDat)) > 2) imDat <- imDat[,,frnum]
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
    i <- as.numeric(read.csv(file_cal, header=F, skip=5, nrow=1))
    nch <- dim(cal$idx)[2]
    nsl <- length(cal$idx)/nch
    i <- t(apply(cal$idx, c(1,3), mean))
    dim(i) <- NULL
  }
  # corner[1] - defines the index offset of the 1st spectral channel
  # corner[2] - defines the index offset of the 1st IO slit
  x <- (corner[1]-50):(corner[1]+nch-1+49)
  sl <- apply(imDat[corner[1]:(corner[1]+nch-1),], 2, mean);
  ch <- apply(imDat[x,], 1, mean);

  if (length(file_eps) != 0)
  {
#     paths <- unlist(strsplit(file_cal, '/')) 
#     file_eps <- gsub(ext, '.eps', paths[length(paths)])
    setEPS(width=10, height=4)
    postscript(file_eps, horizontal=F)
  }
  par(mfrow=c(1,2), mar=c(4,2,1,1), oma=c(1,3,0,0))

  plot(x, ch, type='l', xlab='Spectral channels')
  abline(v=corner[1],       lty=3, col='red')
  abline(v=corner[1]+nch-1, lty=3, col='red')

  plot(sl, type='l', xlab='IO slits')
  for (j in 1:nsl)
  {
    abline(v=corner[2]+i[j], lty=3, col='red')
  }

  title(ylab='Flux (ADU)', line=1, outer=T)
  if (length(file_eps) != 0) dev.off()
  return (list(corner=corner, idx=i, nch=nch, nsl=nsl))
}
