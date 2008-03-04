require(lattice)
mciplot <- function(
    formula,
    data,
    allow.multiple = TRUE,
	as.table = TRUE,
	panel = panel.superpose.mciplot,
	h.displ=0,
	c.width=unit(.1,"inches"),
	c.col=NULL,
	c.lty=NULL,
    guidelines=NULL,...
) {
	levelplot(x=formula, data=data, allow.multiple=allow.multiple, as.table = as.table, panel=panel, contour=FALSE,region=FALSE,h.displ=h.displ,c.width=c.width,c.col=c.col,c.lty=c.lty,guidelines=guidelines,...)
}

panel.superpose.mciplot <- function(x, y, z, groups, subscripts, panel.superpose = panel.mciplot, col,pch,lty,guidelines,h.displ,c.width,c.col,c.lty,...) {
    if(!is.null(guidelines)) {
        par <- trellis.par.get("add.line")
        for(i in 1:length(guidelines)) {
            panel.abline(h=guidelines[i],lty=par$lty,col=par$col,alpha=par$alpha,lwd=par$lwd)
        }
    }
	if (missing(groups)) {
        if(is.null(c.lty)) c.lty <- vector(,length=nvals); c.lty[1:nvals] <- trellis.par.get("superpose.line")$lty[1]
        if(is.null(c.col)) c.col <- vector(,length=nvals); c.col[1:nvals] <- trellis.par.get("superpose.symbol")$col[1]
		panel.superpose(x = x, y = y, z = z, subscripts = subscripts, h.displ=h.displ,c.width=c.width,c.col=c.col,c.lty=c.lty,...)
		return()
	}

	x <- x[subscripts]; y <- y[subscripts]; z <- z[subscripts]
	groups <- groups[subscripts]

	if (is.factor(groups)) {
		vals <- levels(groups)
	} else {
		vals <- sort(unique(groups))
	}
	nvals <- length(vals)

    if(missing(col)) {col <- vector(,length=nvals); col[1:nvals] <- trellis.par.get("superpose.symbol")$col[1:(max(nvals,7))]}
    if(missing(pch)) {pch <- vector(,length=nvals); pch[1:nvals] <- trellis.par.get("superpose.symbol")$pch[1:(max(nvals,7))]}
    if(missing(lty)) {lty <- vector(,length=nvals); lty[1:nvals] <- trellis.par.get("superpose.line")$lty[1:(max(nvals,7))]}

    if(length(col)!=nvals) {length(col) <- nvals; col[is.na(col)] <- col[!is.na(col)]}
    if(length(pch)!=nvals) {length(pch) <- nvals; pch[is.na(pch)] <- pch[!is.na(pch)]}
    if(length(lty)!=nvals) {length(lty) <- nvals; lty[is.na(lty)] <- lty[!is.na(lty)]}

    if(is.null(c.lty)) c.lty <- lty
    if(is.null(c.col)) c.col <- col
    
    displ <- h.displ*(seq(1:nvals)-mean(seq(1:nvals)))
	for (i in seq(along = vals)) {
		id <- (groups == vals[i])
		panel.superpose(x = x[id], y = y[id], z = z[id], groups = groups[id], subscripts = TRUE,h.displ=displ[i],c.width=c.width,c.col=c.col[i],c.lty=c.lty[i],col=col[i],pch=pch[i],lty=lty[i],ylim=ylim,...)
	}
}

panel.mciplot <- function(x,y,z,subscripts,h.displ,c.width,c.col,c.lty,col,pch,lty,...) {
    x <- x[subscripts]; y <- y[subscripts]; z <- z[subscripts];
    means <- y
    upper <- y+z
    lower <- y-z
    for(j in 1:length(x)) {
        panel.arrows(x0=x[j]+h.displ,x1=x[j]+h.displ,y0=means[j],y1=upper[j],angle=90,length = c.width,col=rep(c.col,2),lty=c(c.lty,1))
        panel.arrows(x0=x[j]+h.displ,x1=x[j]+h.displ,y0=means[j],y1=lower[j],angle=90,length = c.width,col=rep(c.col,2),lty=c(c.lty,1))
    }
    panel.lines(x+h.displ,means,col=col,lty=lty)
    panel.points(x+h.displ,means,col=col,pch=pch,fill="white")
}