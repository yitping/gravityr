#' Reorder a vector after an 1D FFT routine
#' 
#' @param m Vector to be reorder
#' 
#' @export
#' 
gv_fftshift <- function(m)
{
  N <- length(m)
  iz <- (N + N%%2)/2
  ms <- m[c((iz+1):N, 1:iz)]
  return(ms)
}

#' Reorder a matrix after an 2D FFT routine
#' 
#' @param m Matrix to be reorder
#' 
#' @export
#' 
gv_fftshift2d <- function(m)
{
  N <- dim(m);
  iz <- (N[2] + N[2]%%2)/2;
  ms <- m[,c((iz+1):N[2], 1:iz)];
  iz <- (N[1] + N[1]%%2)/2;
  ms <- ms[c((iz+1):N[1], 1:iz),];
  return(ms);
}

#' Circular shift
#' 
#' @param a Vector to be shifted
#' @param rev Reverse the direction of shift
#' 
#' @export
#' 
gv_circshift <- function (a, rev=F)
{
	if (!rev)
	{
		ca <- c(a, a[1])
		ca <- ca[-1]
	}
	else
	{
		ca <- c(tail(a,1), a)
		ca <- ca[-length(ca)]
	}
	return (ca)
}

