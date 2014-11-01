#' Parse SCI detector calibration text file
#' 
#' Read the SCI detector calibration file and format the content properly for
#' usage with other GRAVITY R functions, e.g. to convert the SCI raw images to
#' pixel arrays, etc.
#' 
#' @param file The calibration file
#' 
#' @export
#' 
gv_sci_readcal <- function (file, n_io=gv_const()$sci_io_out)
{
  cal <- gv_readcsv(file)
  # reformat idx
  n_pl <- dim(cal$idx)[1]/n_io
  ii <- t(cal$idx)
  dim(ii) <- c(dim(ii)[1], n_io, n_pl)
  cal$idx <- ii
  # readjust corner
  if (cal$s_pix[1] != 0) cal$corner[1] <- cal$corner[1] + cal$s_pix[1]
  return (cal)
}