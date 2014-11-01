#' Parse SCI V2PM text file
#' 
#' Read the SCI detector V2PM file and format the content properly for usage
#' with other GRAVITY R functions, e.g. to convert the SCI pixel arrays to OPDs,
#' etc.
#' 
#' @param file P2VM file
#' 
#' @export
#' 
gv_sci_readv2pm <- function (file)
{
  v2pms <- gv_readcsv(file)
  return (v2pms)
}