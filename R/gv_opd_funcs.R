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
#' @param type Optional. Selects between the following:
#' \item{pd}{Phase delay (default)}
#' \item{gd}{Group delay}
#' \item{op}{OPD (experimental)}
#' @return Estimated OPD.
#' @export
#' 
gv_vis2opd <- function (cv, R=1, lam=1, type='pd', debug=F)
{
  n  <- length(cv)
  n2 <- as.integer(n/2)+1
  gd <- sum(cv[1:(n-1)]*Conj(cv[2:n]), na.rm=T)
  pd <- cumprod(rep(gd/Mod(gd),n))
  pd <- pd*Conj(pd[n2])
	rs <- cv*pd
	p0 <- sum(rs, na.rm=T)
  if (debug)
  {
    plot(Arg(cv), ylim=c(-3,3)); grid()
    lines(Arg(Conj(pd)*p0), col=2)
    points(Arg(rs), col=2)
    points(n2, Arg(p0), pch=0, col=4)
  }
  PD <- Arg(p0)
  GD <- Arg(gd)*(n-1)*R # inverted
	OP <- (GD+PD+gv_modulo_phase(GD))
	# see derivation in logbook (Apr 2014, pg. 2014-08-20)
	var_PD <- var(Arg(rs*Conj(p0)), na.rm=T)
	var_GD <- var_PD/var(1:n)*R
  opd <- switch(type,
                pd=PD/2/pi*lam,
                gd=GD/2/pi*lam,
								op=OP/2/pi*lam)
	var_opd <- switch(type,
										pd=var_PD/2/pi*lam,
										gd=var_GD/2/pi*lam,
										op=var_GD/2/pi*lam)
  return(list(opd=opd, var_opd=var_opd))
}

