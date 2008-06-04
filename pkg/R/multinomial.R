# 
# MAKE FAMILY OBJECT FOR MULTINOMIAL RESPONSES WITH LOGISTIC LINK
# 

mlogit <-
function(base=1) {
	# 	matrix formulation is possibly very inefficient?!?!?
	# 	moreover it does not admit of bases being different from 1??!?!?
	linkfun <- function(p,base) {
		lfun <- function(p,base) {
			p <- p/sum(p)
			beta <- numeric(length(p))
			if(any(p==1)) beta[which(p==1)]=Inf
			else beta[-base] <- log(p[-base]/p[base])
			return(beta)
		}
		if(is.matrix(p)) {
			beta <- t(apply(p,1,lfun,base=base))
		} else {
			beta <- lfun(p,base)
		}
		return(beta)
	}
	linkinv <- function(eta,base) {
		linv <- function(eta,base) {
			pp <- numeric(length(eta))
			if(any(is.infinite(eta)) || any(eta > log(.Machine$double.xmax)) || any(eta < log(.Machine$double.xmin))) {
				pp[which(is.infinite(eta))] <- 1
				pp[which(eta > log(.Machine$double.xmax))] <- 1 # change this to something better!
			} else {
				expb <- exp(eta)
				sumb <- sum(expb)
				pp[base] <- 1/sumb
				pp[-base] <- expb[-base]/sumb
			}
			return(pp)
		}
		if(is.matrix(eta)) {
			if(ncol(eta)==1) {
				pp <- as.matrix(apply(eta,1,linv,base=base)) # fixes problem with column matrix eta
			} else pp <- t(apply(eta,1,linv,base=base)) 	
		} else {
			pp <- linv(eta,base)
		}
		return(pp)
	}
	mu.eta <- function(eta) {
		if(length(eta)==1) return(eta-eta^2)
		if(is.vector(eta)) return(diag(eta)-outer(eta,eta))
	}
	valideta <- function(eta) {
		TRUE # fix me
	}
	
	name <- "mlogit"
	structure(list(linkfun=linkfun,
			linkinv=linkinv,
			mu.eta=mu.eta,
			valideta=valideta,
			name=name,
			base=base),
		class="link-glm")
}



multinomial <-
function(link="mlogit",base=1) {
	# adapted from gaussian()
	linktemp <- substitute(link)
	if (!is.character(linktemp)) {
		linktemp <- deparse(linktemp)
		if (linktemp == "link") {
			warning("use of multinomial(link=link) is deprecated\n",
				domain = NA)
			linktemp <- eval(link)
			if (!is.character(linktemp) || length(linktemp) !=1)
			stop("'link' is invalid", domain = NA)
		}
	}
	okLinks <- c("mlogit")
	if (linktemp %in% okLinks) {
		if(linktemp == "mlogit") stats <- mlogit() else stats <- make.link(linktemp)
	} else {
		if (is.character(link)) {
			stats <- make.link(link)
			linktemp <- link
		} else {
			if (inherits(link, "link-glm")) {
				stats <- link
				if (!is.null(stats$name))
				linktemp <- stats$name
			} else {
				stop(gettextf("link \"%s\" not available for multinomial family; available links are %s",
						linktemp, paste(sQuote(okLinks), collapse = ", ")),
					domain = NA)
			}
		}
	}
	variance <- function(mu) {
		n <- length(mu)
		v <- diag(n)*outer(mu,1-mu) - (1-diag(n))*outer(mu,-mu)
	}
	validmu <- function(mu) {
		all(mu > 0) && all(mu < 1)
	}
	dev.resids <- function(y,mu,wt) {
		
	}
	initialize <- expression()
	structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun,
			linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
			mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu, valideta = stats$valideta, base=base),
			class = "family")
}

