#' Parse SCI P2VM text file
#' 
#' Read the SCI detector P2VM file and format the content properly for usage
#' with other GRAVITY R functions, e.g. to convert the SCI pixel arrays to OPDs,
#' etc.
#' 
#' @param file_p2vm The P2VM file
#' 
gv_sci_readp2vm <- function (file_p2vm, n_io=gv_const()$sci_io_out)
{
  p2vm <- gv_readcsv(file_p2vm)
  n_pl <- dim(p2vm$k)[1]/n_io
  if (n_pl < 2) return (p2vm)
  # reformat these arrays
  for (a in c('Ii', 'Ij', 'ri', 'rj', 'sx', 'k'))
  {
    i <- which(names(p2vm) == a)
    arr <- t(p2vm[[i]])
    dim(arr) <- c(dim(arr)[1], n_io, n_pl)
    p2vm[[i]] <- arr
  }
  return (p2vm)
}