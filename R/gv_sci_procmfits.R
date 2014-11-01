#' Compute GD/PD of SC fringes for a given FITS file
#'
#' @param fits FITS file to be processed
#' @param v2pms List of v2pm matrices
#' @param drk Dark pixels and read noise from gv_sci_Nfits2pix().
#' @param cal A list of calibration coefficients
#' @param p Optional. Polarization states to process.
#' @param hdu Optional. The hdu name that contains the SCI inteferograms.
#' 
#' @export
#' 
gv_sci_procmfits <- function (fits, v2pms, drk, cal, pol=1,
															hdu='imaging_data_sc')
{
	message(sprintf('processing %s [hdu=%s]...', fits, hdu))
	cvf  <- gv_sci_fits2vis(fits, v2pms, drk, cal, p=pol, hdu=hdu)
	flux <- sapply(cvf, function (cvi) cvi$flux)
	N  <- dim(cvf[[1]]$var_pd)[2]
	pd <- list( mu=sapply(cvf, function (cvi) apply(cvi$vis, 1, sum)),
						 var=sapply(cvf, function (cvi) apply(cvi$var_pd, 1, sum)/N^2))
  # TODO: gv_vis2opd is still a hack for now [2014-10-31]
	gd <- list( mu=sapply(cvf, function (cvi)
												apply(cvi$vis, 1, function (x)
															gv_vis2opd(x, R=1, lam=1, type='gd')$opd)),
						 var=sapply(cvf, function (cvi)
												apply(cvi$vis, 1, function (x)
															gv_vis2opd(x, R=1, lam=1, type='gd')$var_opd)))
	v2 <- Mod(pd$mu)^2
	v2 <- v2/sapply(cvf, function (cvi) apply(Mod(cvi$vis), 1, sum))^2
	cp <- apply(pd$mu, 2, gv_vis2cp)
	return (list(flux=flux, v2=v2, pd=pd, gd=gd, cp=cp))
}

