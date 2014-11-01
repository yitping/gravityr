#' Compute weighted mean and variance of column vectors/matrices
#'
#' @param x Values
#' @param v Variance of each value
#' 
#' @export
#' 
gv_weighted.mean <- function (x, v)
{
	sz <- dim(x)
	wt <- apply(1/v, 1, sum, na.rm=T)
	wt_x <- apply(x/v, 1, sum, na.rm=T)/wt
	wt_v <- apply((x-matrix(rep(wt_x,sz[2]), nrow=sz[1], ncol=sz[2]))^2/v, 1, sum, na.rm=T)/wt
	return (list(m=wt_x, v=wt_v))
}

#' Compute the close loop OPD between the FT, SC and MET fringes
#'
#' @param fits Filename of the GRAVITY merged FITS
#' @param v2pms A list of V2PM coefficients to reduce SC fringes
#' @param drk Dark pixels and read noise from gv_sci_Nfits2pix().
#' @param cal A list of calibration coefficients
#' @param p (optional) Polarization states to process.
#' @param lam (optional) A vector of mean wavelengths for the FT, SC and MET fringes in microns.
#' 
#' @export
#' 
gv_calib_intopd <- function (fits, v2pms, drk, cal, pol=1, debug=F, high.res=F,
														 ka=c(0,2,2,0),
														 hdu_sc='imaging_data_sc',
														 hdu_ft='imaging_data_ft',
														 hdu_met='metrology',
														 hdu_ddl='fddl')
{
	# compute the complex visibilities of the SC fringes and all other relevant
	# parameters
	all_sc<- gv_sci_procmfits(fits, v2pms, drk, cal, pol=pol, hdu=hdu_sc)
	pd_sc <- list( mu=t(apply(Arg(all_sc$pd$mu), 1, unwrap)), var=all_sc$pd$var)

	# read the GD of the FT fringes
	all_ft<- gv_ft_readmfits(fits)
	ts_ft <- gv_readfits(fits, hdu=hdu_ft, col='time')

	# compute the complex visibility of the MET fringes
	# CV data starts from Nth index
	cv_met<- gv_met_procmfits(fits, sel=17:20, ka=ka)
	ts_met<- gv_readfits(fits, hdu=hdu_met, col='time')

	# read the FDDL commands
	message(sprintf('reading %s [hdu=%s]...', fits, hdu_ddl))
	cm_ddl<- t(cbind(gv_readfits(fits, hdu=hdu_ddl, col='ft_ddl_cmd'),
									 gv_readfits(fits, hdu=hdu_ddl, col='sc_ddl_cmd')))
	ts_ddl<- gv_readfits(fits, hdu=hdu_ddl, col='time')


	N_bl  <- dim(all_sc$pd$mu)[1]
	N_sc  <- dim(all_sc$flux)[2]

	# align ts=0
	N_ft  <- length(ts_ft[-which(ts_ft<0)])
	N_met <- length(ts_met[-which(ts_met<0)])
	N_ddl <- length(ts_ddl[-which(ts_ddl<0)])

	# how many data per bin
	bin_ft  <- floor(N_ft/N_sc)
	bin_met <- floor(N_met/N_sc)
	bin_ddl <- floor(N_ddl/N_sc)

	# find the indices of data where each time series should take
	ix_ft  <- 1:(bin_ft*N_sc) + tail(which(ts_ft<0),1)
	ix_met <- 1:(bin_met*N_sc)+ tail(which(ts_met<0),1)
	ix_ddl <- 1:(bin_ddl*N_sc)+ tail(which(ts_ddl<0),1)

	ts_ft_sel <- ts_ft[ix_ft]
	pd_ft_sel <- all_ft$pd[,ix_ft]
	gd_ft_sel <- all_ft$gd[,ix_ft]
	dim(pd_ft_sel) <- dim(gd_ft_sel) <- c(N_bl, bin_ft, N_sc)
	dim(ts_ft_sel) <- c(bin_ft, N_sc)

	ts_met_sel<- ts_met[ix_met]
	cv_met_sel<- lapply(1:4, function (i)
											{
												a <- cv_met[i,ix_met]
												dim(a) <- c(bin_met, N_sc)
												return (a)
											})
	cv_met_sel<- simplify2array(cv_met_sel, higher=T)
	dim(ts_met_sel)<- c(bin_met, N_sc)
	
	ts_ddl_sel<- ts_ddl[ix_ddl]
	cm_ddl_sel<- cm_ddl[,ix_ddl]
	dim(cm_ddl_sel)<- c(dim(cm_ddl)[1], bin_ddl, N_sc)
	dim(ts_ddl_sel)<- c(bin_ddl, N_sc)

	ts_ft_binned   <- apply(ts_ft_sel, 2, mean)
	temp <- apply(pd_ft_sel, c(1,3), sum, na.rm=T)
	pd_ft_binned   <- list( mu=t(apply(Arg(temp), 1, unwrap)),
												 var=sapply(1:N_sc, function (i)
																		sapply(1:N_bl, function (j)
																					 var(Arg(pd_ft_sel[j,,i]*Conj(temp[j,i])), na.rm=T)
																					 ))
												 )
	gd_ft_binned   <- list( mu=apply(gd_ft_sel, c(1,3), mean),
												 var=apply(gd_ft_sel, c(1,3), var))

	ts_met_binned  <- apply(ts_met_sel, 2, mean)
	cv_met_binned  <- apply(cv_met_sel, 3, function (cvi)
													{
														m <- apply(cvi, 2, sum, na.rm=T)
														v <- sapply(1:length(m), function (j)
																				{
																					var(Arg(cvi[,j]*Conj(m[j])), na.rm=T)
																				})
														list(mu=m, var=v)
													})
	pd_met_binned  <- gv_met_compdcv(cv_met_binned)
	pd_met_binned  <- list( mu=t(apply(pd_met_binned$mu, 2, function (pdi)
																		 unwrap(Arg(pdi)))),
												 var=t(pd_met_binned$var))
	td_met_binned  <- list( mu=t(Arg(sapply(cv_met_binned, function (cvi) cvi$mu))),
												 var=t(sapply(cv_met_binned, function (cvi) cvi$var)))

	ts_ddl_binned  <- apply(ts_ddl_sel, 2, mean)
	cm_ddl_binned  <- list( mu=apply(cm_ddl_sel, c(1,3), mean),
												 var=apply(cm_ddl_sel, c(1,3), var))

	# maximum error in the time stamp of each bin between FT and MET
	ts_bin <- rbind(ts_ft_binned, ts_met_binned, ts_ddl_binned)
	ts_err <- max(apply(ts_bin, 2, function (tsi) diff(range(tsi))))*1e-3

	if (!high.res)
	{
		return (list(pd_sc=gv_weighted.mean(pd_sc$mu, pd_sc$var),
								 pd_ft=gv_weighted.mean(pd_ft_binned$mu, pd_ft_binned$var),
								 pd_met=gv_weighted.mean(pd_met_binned$mu, pd_met_binned$var),
								 gd_sc=gv_weighted.mean(all_sc$gd$mu*2*pi, all_sc$gd$var*2*pi),
								 gd_ft=gv_weighted.mean(gd_ft_binned$mu, gd_ft_binned$var),
								 cm_ddl=list(mu=apply(cm_ddl_binned$mu, 1, mean), var=rep(0,dim(cm_ddl_binned$var)[1])),
								 td_met=gv_weighted.mean(td_met_binned$mu, td_met_binned$var),
								 ts_err=ts_err))
	}
	else
	{
		return (list(pd_sc=pd_sc,
								 pd_ft=pd_ft_binned,
								 pd_met=pd_met_binned,
								 gd_sc=lapply(all_sc$gd, function (gdi) gdi*2*pi),
								 gd_ft=gd_ft_binned,
								 cm_ddl=cm_ddl_binned,
								 td_met=td_met_binned,
								 ts_err=ts_err))
	}
}


