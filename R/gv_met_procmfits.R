#' Compute complex visibility for a time series of voltages
#'
#' @param volt Vector of diode voltages
#' @param N Number of shifts per phase estimates (must be multiple of 4)
#' @param abcd Initial phase shift
#' @param ka Phase shift calibration
#' 
#' @export
#' 
gv_met_volt2vis <- function (volt, N=4, abcd=1:4, ka=c(0,2,2,0))
{
	vis <- rep(NA, length(volt))
	idx <- gv_circshift(abcd)
	for (i in N:length(volt))
	{
		idx <- gv_circshift(idx, rev=T)
		v <- volt[(i-N+1):i]
		dim(v) <- c(4,N/4)
		u <- apply(v, 1, mean)
		vis[i] <- gv_abcd2vis(u[idx], kx=ka, sx=c(4,4,4,4), fx=sum(u))
	}
	return (vis)
}

#' Compute PD of MET fringes
#'
#' @param fits FITS file to process
#' @param sel Optional. Columns within a HDU/volt table/column to process
#' @param Nps Optional. Number of shifts per phase estimates (must be multiple of 4)
#' @param Ncpu Optional. Use multiple processors
#' @param ka Optional. Phase shift calibration
#' @param hdu Optional. The name of the HDU to process
#' 
#' @export
#' 
gv_met_procmfits <- function (fits, sel=1:4, Nps=4, Ncpu=1,
															ka=c(0,2,2,0), hdu='metrology')
{
	message(sprintf('processing %s [hdu=%s]...', fits, hdu))
	volts <- gv_readfits(fits, hdu=hdu, col='volt')[,sel]
	# find where is 'A'
	phset <- gv_readfits(fits, hdu=hdu, col='volt_ps_set')[1:4]
	ps    <- 1:4
	phA   <- which(phset == 0)
	if (length(phA) == 0) stop('fatal error: metrology phase A not found')
	if (phA > 1)
	{
		for (i in 1:(phA-1)) { ps <- gv_circshift(ps) }
	}

	if (Ncpu > 1)
	{
		cv <- mclapply(1:length(sel), function (j)
									 gv_met_volt2vis(volts[,j], N=Nps, abcd=ps, ka=ka),
									 mc.cores=Ncpu)
	}
	else
	{
		cv <- lapply(1:length(sel), function (j)
								 gv_met_volt2vis(volts[,j], N=Nps, abcd=ps, ka=ka))
	}
	return (t(simplify2array(cv, higher=T)))
}

