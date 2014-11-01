#' Parse FT P2VM text file
#' 
#' @param file P2VM file
#' 
#' @export
#'
gv_ft_readp2vm <- function (file)
{
	cal  <- gv_readcsv(file)
	p2vm <- c(lapply(cal[2:7], function (m) { dim(m) <- c(24,16); m }),
						cal['Dark'])
	return (p2vm)
}
