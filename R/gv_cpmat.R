#' Closed triangle matrix
#' 
#' Maps a vector of baselines to a vector of closed triangles.
#' 
#' @param Optional. Flip the matrix.
#' 
#' @export
#' 
gv_cpmat <- function (flip=F)
{
  # assumes baseline vec = 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
  M <- rbind(c( 1,-1, 0, 1, 0, 0),
						 c( 1, 0,-1, 0, 1, 0),
						 c( 0, 0, 0, 1,-1, 1))
  if (flip) M <- t(apply(M, 1, rev))
  return(M)
}
