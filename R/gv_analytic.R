#' Compute the analytic signal
#' 
#' @param sig The signal
#' @param fc Optional. Cut-off frequency of a low-pass filter.
#'
#' @export
#' 
gv_analytic <- function (sig, fc=NULL, debug=F)
{
	N  <- length(sig)
	Ny <- (N+(N%%2))/2
	ma <- c(0, rep(2,Ny-1), rep(0,N-Ny))
	if (!is.null(fc) && (fc > 2) && (fc < Ny))
	{
		ma <- ma*c(rep(1,fc), rep(0,N-(fc+1)+1))
	}
	S  <- ma*fft(sig)
	if (debug) { plot(log10(Mod(S[1:Ny])), type='l') }
	return (fft(S, inv=T)/N)
}

