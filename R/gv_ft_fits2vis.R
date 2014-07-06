#'
#'
gv_ft_fits2vis <- function (fits, p2vm,
														dark=0, p=1, hdu='imaging_data_ft',
														blseq=NULL, ma=5)
{
	tel_idx <- 1:4
	flux <- gv_readfits(fits, hdu, 'pix')
	pol_sel <- (p-1)*120+(1:120)
	if ((length(dark) == 1) && (dark == 0) &&
			(length(which(names(p2vm) == 'Dark')) != 0))
	{
		dark <- p2vm$Dark[pol_sel]
	}
	#vis_pre <- lapply(1:dim(flux)[1], function (i)
	#									gv_ft_pix2vis(flux[i,pol_sel], dark, p2vm, blseq=blseq))
	#vis_pre <- simplify2array(vis_pre, higher=T)
	#vis <- lapply(ma:dim(flux)[1], function (i)
	#							apply(vis_pre[,,(i-ma+1):i], c(1,2), sum))
	#vis <- simplify2array(vis, higher=T)
	vis <- lapply(ma:dim(flux)[1], function (i)
								gv_ft_pix2vis(apply(flux[(i-ma+1):i,pol_sel], 2, mean),
															dark, p2vm, blseq=blseq))
	vis <- simplify2array(vis, higher=T)
}
