#' Telescope matrix
#'
#' Maps a vector of telescopes to a vector of baselines.
#' By default, the order of telescope pair follows the arrangement as seen on
#' the SCI detector.
#'
#' @param map  Telescope vector
#'
#' @export
gv_telmat <- function (map=gv_const()$sci_telarr, flip=T)
{
  N <- length(map)
  M <- do.call('rbind', sapply(1:(N-1), function (i)
    t(sapply((i+1):N, function (j)
      {
      m <- rep(0,N)
      m[c(map[i],map[j])] <- c(1,-1)
      return (m)
      }))))
  # the current image of the IO is "upside down"
  if (flip) M <- apply(M, 2, rev)
  return (M)
}
