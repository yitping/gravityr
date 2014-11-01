#' Compute the differential OPD of the metrology path
#'
#' @param cv A list of OPD per telescope
#' 
#' @export
#' 
gv_met_compdcv <- function (cv)
{
	n <- length(cv)
  m <- do.call('cbind', sapply(1:(n-1), function (i)
                               sapply((i+1):n, function (j)
                                      cv[[i]]$mu*Conj(cv[[j]]$mu))
                               ))
	v <- do.call('cbind', sapply(1:(n-1), function (i)
															 sapply((i+1):n, function (j)
																			cv[[i]]$var+cv[[j]]$var)
															 ))
	return (list(mu=m, var=v))
}
