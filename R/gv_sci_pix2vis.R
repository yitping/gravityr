#' Computes the complex visibilty of fringes
#' 
#' Computes the complex visibilty of the SCI fringes in a pixel array
#' 
#' @export
gv_sci_pix2vis <- function (pixels, rdnoiz, v2pms, blseq=NULL)
{
  coh <- gv_sci_pix2vis_unsorted(pixels, rdnoiz, v2pms)
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