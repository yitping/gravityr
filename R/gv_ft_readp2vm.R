#' Parse FT P2VM text file
#' 
#'
gv_ft_readp2vm <- function (file_p2vm)
{
	cal  <- gv_readcsv(file_p2vm)
	p2vm <- c(lapply(cal[2:7], function (m) { dim(m) <- c(24,16); m }),
						cal['Dark'])
	return (p2vm)
}
