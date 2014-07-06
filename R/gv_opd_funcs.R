#' Rounding up phases
#' 
#' @param p A variable in radians
#' @export
#' 
gv_round_phase <- function (p)
{
  2*pi*((floor(abs(p)/pi) + floor(abs(p)/pi)%%2)/2)*(p/abs(p))
}

#' @export
#' 
gv_modulo_phase <- function (p) { p-gv_round_phase(p) }

#' Estimate OPD from spectrally dispersed fringes.
#' 
#' Visibility
#' 
#' @param cv Complex visibility vector
#' @param R Optional. Spectral resolution of the spectrometer.
#' @param lam Optional. Wavelength in the middle of the spectrum.
#' @param type Optional. Selects the gallery.
#' @return Estimated OPD.
#' @export
#' 
gv_vis2opd <- function (cv, R=1, lam=1, type='pd', debug=F)
{
  n  <- length(cv)
  n2 <- as.integer(n/2)+1
  gd <- sum(cv[1:(n-1)]*Conj(cv[2:n]))
  pd <- cumprod(rep(gd,n))
  pd <- pd*Conj(pd[n2])
  ### gvspc ###
  #unph <- unwrap(Arg(cv))
  #gd <- median(unlist(sapply(1:(n-1), function (i)
  #													 sapply((i+1):n, function (j)
  #																	(unph[j]-unph[i])/(j-i)))))
  #pd <- exp(complex(im=-gd*(0:(n-1)-as.integer(n/2))))
  if (debug)
  {
    plot(Arg(cv), ylim=c(-3,3)); grid()
    lines(Arg(Conj(pd)*sum(cv*pd)), col=2)
    points(Arg(cv*pd), col=2)
    points(n2, Arg(sum(cv*pd)), pch=0, col=4)
  }
  GD <- Arg(gd)*(n-1)*R # inverted
  PD <- Arg(sum(cv*pd))
  opd <- switch(type,
                gd=GD/2/pi*lam,
                pd=(GD+PD+gv_modulo_phase(GD))/2/pi*lam)
  return(opd)
}