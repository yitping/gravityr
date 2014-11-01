#' Read GD/PD of FT fringes for a given FITS file
#'
#' @param fits A FITS file to be read
#' @param hdu Optional. Name of HDU in the FITS file to be read
#' 
#' @export
#' 
gv_ft_readmfits <- function (fits, hdu='imaging_data_ft')
{
  message(sprintf('reading %s [hdu=%s]...', fits, hdu))

	p <- gv_readfits(fits, hdu=hdu, col='pd')
	p <- t(apply(p, 2, function (pdi) exp(complex(imag=pdi))))
	g <- gv_readfits(fits, hdu=hdu, col='gd')
	g <- t(g%*%diag(c(-1,-1,1,1,1,1)))  # maps to 1:2, 1:3, 1:4, 2:3, 2:4, 3:4

	#pd <- list(mu=p, var=NA)
	#gd <- list(mu=g, var=NA)
	#cp <- apply(pd$mu, 2, gv_vis2cp)

  return (list(flux=0, v2=0, pd=p, gd=g))
}


