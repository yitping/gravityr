#' Compute the power spectral density (PSD) of a signal
#'
#' The DC component is removed. The PSD starts from frequency fs/(N-1)
#'
#' @param sig Signal
#'
#' @export
#'
gv_psd <- function (sig)
{
	n   <- length(sig)
	ftv <- fft(sig)[1:(n/2)]
	psd <- 2*(Mod(ftv[-1])/n)^2
	return (psd)
}
