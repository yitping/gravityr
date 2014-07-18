#' Computes the complex visibilty of fringes
#' 
#' Computes the complex visibilty of the SCI fringes in a pixel array
#' 
#' @export
gv_sci_pix2vis <- function (flux, dark, p2vm, blseq=NULL)
{
  vis <- gv_sci_pix2vis_unsorted(flux, dark, p2vm)
  if (!is.list(blseq))
  {
    blseq <- gv_blorder(gv_telmat())
    # verified with CP
    bl_seq$sign[2] <- -1
    vis <- vis[,blseq$order]
    j   <- which(blseq$sign == -1)
    if (length(j) != 0) { vis[,j] <- Conj(vis[,j]) }
  }
  return (vis)
}