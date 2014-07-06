#'
#'
gv_ft_pix2vis <- function (flux, dark, p2vm, blseq=NULL)
{
	if (!is.list(blseq))
	{
		blseq <- gv_blorder(gv_telmat(c(2,3,1,4)))
		# verified with CP
		blseq$sign[c(1,2)] <- 1
	}
	px <- flux - dark
	dim(px) <- c(4, 6, 5)
	px <- px[,blseq$order,]
  #fl <- t(apply(px, c(2,3), sum))
	dim(px) <- c(24, 5)
	# column arrangement of px is good - check with flux of all spec. channel
	pix <- cbind(apply(px,1,sum), px)
	phr <- sapply(1:6, function (i) pix[,i]%*%p2vm[[i]])
	bsl_idx <- 5:10
	vis <- t(sapply(2:dim(phr)[2], function (i)
								complex(real=phr[bsl_idx,i], imag=blseq$sign*phr[bsl_idx+6,i])))
	return(vis)
}
