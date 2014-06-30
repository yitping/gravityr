#' Remap an arbitrary baseline arrangement to the default
#' 
#' A default baseline arrangment in GRAVITY is given as:
#' T1-T2, T1-T3, T1-T4, T2-T3, T2-T4, T3-T4
#' 
#' @param M telescope matrix of size 6x4
#' 
#' @export
gv_blorder <- function (M)
{
  M_default <- rbind(c( 1,-1, 0, 0),
                     c( 1, 0,-1, 0),
                     c( 1, 0, 0,-1),
                     c( 0, 1,-1, 0),
                     c( 0, 1, 0,-1),
                     c( 0, 0, 1,-1))
  idx <- apply(abs(M)%*%abs(t(M_default)), 2, function (m) which(m==2))
  pol <- diag(M[idx,]%*%t(M_default))/2
  return(list(order=idx, sign=pol))
}