#' Parse SCI V2PM text file
#' 
#' Read the SCI detector V2PM file and format the content properly for usage
#' with other GRAVITY R functions, e.g. to convert the SCI pixel arrays to OPDs,
#' etc.
#' 
#' @param file_p2vm The P2VM file
#' 
gv_sci_readv2pm <- function (file_v2pm)
{
  v2pms <- gv_readcsv(file_v2pm)
  return (v2pms)
}