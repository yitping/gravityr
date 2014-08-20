#' Computes the complex visibilty of fringes
#' 
#' Computes the complex visibilty of the SCI fringes in a pixel array.
#' 
#' This is a wrapping function for \code{\link{gv_sci_pix2vis_unsorted}} and performs
#' some post-processing to the computed complex visibility.
#' 
#' The nominal baseline order of the SCI spectrometer is given as: ...
#' 
#' @inheritParams gv_sci_pix2vis_unsorted
#' @param blseq Optional. An vector that contains indices to rearrange the baseline order.
#' 
#' @export
gv_sci_pix2vis <- function (pixels, v2pms, rdnoiz, bad=NULL, blseq=NULL)
{
  if (!is.null(bad)) coh <- gv_sci_pix2vis_unsorted(pixels, v2pms, rdnoiz, bad)
  else  coh <- gv_sci_pix2vis_unsorted(pixels, v2pms, rdnoiz)
  if (!is.list(blseq))
  {
    blseq <- gv_blorder(gv_telmat())
    # verified with CP
    blseq$sign[2] <- -1
  }
  coh$vis <- coh$vis[blseq$order,]
  coh$var_v2 <- coh$var_v2[blseq$order,]
  coh$var_pd <- coh$var_pd[blseq$order,]
  i   <- which(blseq$sign == -1)
  if (length(i) != 0) { coh$vis[i,] <- Conj(coh$vis[i,]) }
  return (coh)
}