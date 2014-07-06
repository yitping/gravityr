#' Compute closure phases or OPDs
#' 
#' By default the vector baseline whose closure phases or OPDs are to be computed
#' are arranged in the default GRAVITY sequence, i.e. T1-T2, T1-T3, and so on.
#' 
#' @param vis Complex visiblity of a vector of baselines.
#' @param opd OPD of a vector of baselines.
#' @return A vector of closure phases or OPDs
#' 
#' @export
gv_vis2cp <- function (vis, cpi=rbind(c(1,4,2), c(1,5,3), c(4,6,5)))
{
  apply(cpi, 1, function (bli) vis[bli[1]]*vis[bli[2]]*Conj(vis[bli[3]]))
}

#' @export
gv_opd2cp <- function (opd, cpi=rbind(c(1,4,2), c(1,5,3), c(4,6,5)))
{
  apply(cpi, 1, function (bli) opd[bli[1]]+opd[bli[2]]-opd[bli[3]])
}

