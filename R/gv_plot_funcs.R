#' Plot closure phase
#' 
#' @param cp Matrix of closure phases. One triangle per row.
#' 
#' @export
#' 
gv_plot_cp <- function (cp, x=NULL, d=1, e_cp=NULL, epsfile=NULL, ...)
{
	mu <- apply(cp, d, sum)
	y  <- t(apply(cp, d, function (cpi) Arg(cpi*Conj(sum(cpi)))))
	n  <- dim(cp)[2]
	ylims <- c(floor(-pi),ceiling((3*2-1)*pi))
	if (is.null(x)) { x <- 1:n }
	ss <- sqrt(apply(y^2, d, sum)/(n-1))
	if (!is.null(epsfile))
	{
		setEPS(width=16, height=8)
		postscript(epsfile, horizontal=F)
	}
	if (is.null(e_cp)) plot(x, y[1,], ylim=ylims, ...)
	else gv_plot_errorbars(x, y[1,], ey=e_cp[1,], ylim=ylims, ...)
	grid()
	message(paste('sd cp[', 1, '] = ',
								format(ss[1], digits=2), ' (',
								format(ss[1]/pi*180, digits=2), ')', sep=''))
	for (i in 2:3)
	{
		offset <- (i-1)*2*pi
		if (is.null(e_cp)) points(x, y[i,]+offset, col=i, ...)
		else {
			par(new=T);
			gv_plot_errorbars(x, y[i,]+offset, ey=e_cp[i,], ylim=ylims, col=i, ...)
		}
		abline(h=offset, col=i, lty=3)
		message(paste('sd cp[', i, '] = ',
									format(ss[i], digits=2), ' (',
									format(ss[i]/pi*180, digits=2), ')', sep=''))
	}
	legend('bottomright', bg='white', seg.len=1,
				 format(Arg(mu), digit=3), col=1:3, pch=1, cex=0.6)
	if (!is.null(epsfile)) { dev.off() }
}

#' Plot GD or GD+PD
#' 
#' @param opd Data to be plot. One set per row.
#' 
#' @export
#' 
gv_plot_opd <- function (opd, x=NULL, d=1, offset=0, e_opd=NULL, epsfile=NULL, noleg=F, ...)
{
	mu <- apply(opd, d, mean)
	ss <- apply(opd, d, sd)
	n  <- dim(opd)[2]
	m <- length(mu)
	ylims <- c(ifelse(offset == 0,min(mu)-3*ss[which.min(mu)],-3*ss[1]),
						 ifelse(offset == 0,max(mu)+3*ss[which.max(mu)],(m-1)*offset+3*ss[m]))
	if (is.null(x)) { x <- 1:n }
	if (offset == 0) y <- opd
	else y <- opd - cbind(mu)%*%rbind(rep(1,dim(opd)[2]))
	if (!is.null(epsfile))
	{
		setEPS(width=16, height=8)
		postscript(epsfile, horizontal=F)
	}
	if (is.null(e_opd)) plot(x, y[1,], ylim=ylims, ...)
	else gv_plot_errorbars(x, y[1,], ey=e_opd[1,], ylim=ylims, ...)
	abline(h=ifelse(offset == 0, mu[1], 0), col=1, lty=3)
	grid()
	for (i in 2:m)
	{
		b <- (i-1)*offset
		if (is.null(e_opd)) points(x, y[i,]+b, col=i, ...)
		else {
			par(new=T)
			gv_plot_errorbars(x, y[i,]+b, ey=e_opd[i,], ylim=ylims, col=i, ...)
		}
		abline(h=ifelse(offset == 0, mu[i], b), col=i, lty=3)
	}
	if (!noleg) legend('bottomright', bg='white', seg.len=1,
										 format(mu, digit=3), col=1:m, pch=1, cex=0.6)
	if (!is.null(epsfile)) { dev.off() }
}

#' Plot errorbars
#' 
#' @param x x
#' @param y y
#' @param ex Error bar for x
#' @param ey Error bar for y
#' 
#' @export
#' 
gv_plot_errorbars <- function (x, y, ex=NULL, ey=NULL, w=0.05, ...)
{
  #cl <- match.call()
  #print(cl)
  #args <- c(list(x=x, y=y), list(...))
  pch_num_id <- which(names(list(...)) == 'pch')
  if (length(pch_num_id) == 0)
  {
    pch_num <- NULL
  } else {
    pch_num <- list(...)[[pch_num_id]]
  }

  #plot(x, y, par_user)
  #do.call('plot', c(args, list(pch=3)))
  plot(x, y, ...)
  if (!is.null(ex))
  {
    x0 <- x-ex
    x1 <- x+ex
    arrows(x0, y, x1, y, code=3, angle=90, length=w, ...)
  }
  if (!is.null(ey))
  {
    y0 <- y-ey
    y1 <- y+ey
    arrows(x, y0, x, y1, code=3, angle=90, length=w, ...)
  }
  if (!is.null(ex) || !is.null(ey))
  {
    #points(x, y, par_user)
    #do.call('points', c(args, list(pch=pch_num)))
    #if (is.null(pch_num)) { points(x, y) }
    #else { points(x, y, pch=pch_num) }
    points(x, y, ...)
  }
}

