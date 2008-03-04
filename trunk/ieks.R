#####################################################
# ITERATED EXTENDED KALMAN FILTER
# 
# Families supported
#  - Binomial / identity link
#  - Binomial / logit link
#  - Binomial / probit link
#  - Poisson / identity link
#  - Poisson / log link
#  - Gaussian / identity link
#  - Multinomial / canonical link
#  - Multinomial / pom link
#
# ASSUMING MULTIVARIATE LATENT
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

ieks <- function(ssm, m.start = NA, max.iter = 50, eps = 0.0001)
{
	# require
	# ssm: object of class state space model
	# m.start: optional. Used as starting point in the taylor expansion
	# max.iter: Maximum number of iteration (default 50)
	# eps: Convergence criteria (default 0.0001)
	
	# Returns updated ssm object with
	# mt: Filtered states
	# Ct: Filtered MSE
	# m.tilde: Smoother states
	# C.tilde: Smoothed MSE
	# llh: log-likelihood (unadjusted)
	
	if (is.na(ssm$Yt[1])){
		cat("NOTE: No observations supplied\n")
		cat("      Generating observation using simulate.ssm\n")
		ssm <- simulate.ssm(ssm)
		cat("      Done.\n\n")
		}
	cat("Beginning iteration\n\n")
	msge <- ""
	count <<- 0
	obj <- 10^10
	goon <- TRUE
	if (!is.na(ssm$smoothed$m.tilde[1])) m.start<- ssm$smoothed$m.tilde
	
	objective <- function(ssm, m.start, eps, nosilent=F)
	{
	if(is.na(m.start[1])) obj <- eps + 100
	else obj <- max(abs((ssm$smoothed$m.tilde - m.start)/m.start))
	if (nosilent)
	{
		count <<- count + 1
		cat("Iteration",count,"\n")
		cat("Objective:",obj,"should be smaller than",eps,"\n")
		cat(paste("Log-likelihood:",ssm$filtered$llh,"\n"))
	}	
	obj
	}
	
	#####################################################
	# BINOMIAL / LOGIT LINK CASE
	#####################################################
	if(ssm$fam == "binomial" & ssm$link == "logit") {
		ssm <- kalman.smoother(kal.fil.bin.logit(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.bin.logit(ssm))
			}
	}
	#####################################################
	# BINOMIAL / PROBIT LINK CASE
	#####################################################
	if(ssm$fam == "binomial" & ssm$link == "probit") {
		ssm <- kalman.smoother(kal.fil.bin.probit(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.bin.probit(ssm))
			}
	}
	#####################################################
	# BINOMIAL / IDENTITY LINK CASE
	#####################################################
	if(ssm$fam == "binomial" & ssm$link == "identity") {
		ssm <- kalman.smoother(kal.fil.bin.id(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.bin.id(ssm))
			}
	}
	#####################################################
	# POISSON / LOG LINK CASE
	#####################################################
	if(ssm$fam == "Poisson" & ssm$link == "log") {
		ssm <- kalman.smoother(kal.fil.po.log(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.po.log(ssm))
			}
	}
	#####################################################
	# POISSON / IDENTITY LINK CASE
	#####################################################
	if(ssm$fam == "Poisson" & ssm$link == "identity") {
		ssm <- kalman.smoother(kal.fil.po.id(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.po.id(ssm))
			}
	}
	#####################################################
	# GAUSSIAN / IDENTITY LINK CASE
	#####################################################
	if(ssm$fam == "Gaussian" & ssm$link == "identity") {
		ssm <- kalman.smoother(kalman.filter(ssm))
		msge <- paste("No iterations for Gaussian model\n")
	}
	#####################################################
	# MULTINOMIAL / CANONICAL LINK CASE
	#####################################################
	if(ssm$fam == "multinomial" & ssm$link == "canonical") {
		ssm <- kalman.smoother(kal.fil.mul.can(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.mul.can(ssm))
			}
	}
	#####################################################
	# MULTINOMIAL / POM LINK CASE
	#####################################################
	if(ssm$fam == "multinomial" & ssm$link == "pom") {
		ssm <- kalman.smoother(kal.fil.mul.pom(ssm))
		while(objective(ssm,m.start,eps,T)	>eps & count<max.iter) {
			m.start <- ssm$smoothed$m.tilde
			ssm <- kalman.smoother(kal.fil.mul.pom(ssm))
			}
	}
	
	#####################################################
	# WRAP UP
	#####################################################

	if(objective(ssm,m.start,eps,F) < eps) msge <- paste("Converged in", count, "iterations\n")
	if(count > max.iter) msge <- paste("WARNING! max number of iterations (", max.iter, ")  exceeded\n", sep = "")
	if(msge == "") msge <- paste("WARNING! Distribution:", ssm$fam, "with link:", ssm$link, "not supported.\n")

	cat(msge)
	cat("\nIterated extended Kalman filter and smoother finished.\n")
	ssm
}

#####################################################
# SSM - A STATE SPACE MODEL OBJECT
# 
# FOR USE WITH THE 
# ITERATED EXTENDED KALMAN FILTER AND SMOOTHER
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################


ssm <- function(Ft=function(i,x,phi){},
                Gt=function(i,x,phi){},
                Vt=NA,
                Wt=function(i,x,phi){},
                Xt=NA,
                phi=NA,
                m0,
                C0,
                Yt=NA,
                nt=NA, 
                fam="Gaussian", 
                link="identity",
                m.start=NA)
{
	result <- list(Ft=Ft,
	               Gt=Gt,
	               Vt=Vt,
	               Wt=Wt,
	               Xt=Xt,
	               phi=phi,
	               m0=m0,
	               C0=C0,
	               Yt=Yt,
	               nt=nt,
	               fam=fam,
	               link=link,
	               m.start=m.start,
	               filtered=list(mt=NA,Ct=NA,llh=NA),
	               smoothed=list(m.tilde=NA,C.tilde=NA),
	               simulated=list(theta=NA, lambda=NA)
	               )
	class(result) <- "ssm"
	result
	}

#####################################################
# SIMULATE.SSM
# 
# SIMULATES OBSERVATIONS BASED ON A SSM OBJECT
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

simulate.ssm <- function(ssm, n=50)
{
# requires
# ssm: object of class state space model
# n  : number of observations to generate
# note: No part of latent process allowed to be deterministic

if(!is.na(ssm$Yt[1])) {
	cat("You cannot simulate observations for this \n") 
	cat("state space model since observations are \n")
	cat("already associated!\n\n")
	break}


# Now generate observations

if (ssm$fam!="multinomial")
{
	# generate latent process and signal
	p <- length(ssm$m0)
	theta <- matrix(NA, ncol=p, nrow=n+1)	
	lambda <- rep(NA, n+1)
	theta[1,] <- ssm$m0
	for(i in 2:(n+1))
	{
		theta[i,] <- ssm$Gt(i,Xt[i,],phi)%*%theta[i-1,] + t(rmvnorm(1, mean=rep(0,p), ssm$Wt(i,Xt[i,],phi)))
		lambda[i] <- t(ssm$Ft(i,Xt[i,],phi))%*%theta[i,]
	}
	theta <- theta[-1,]
	lambda <- lambda[-1]
	
	if (ssm$fam=="Gaussian" & ssm$link=="identity")
	{
		ssm$Yt <- rep(NA, n)
		for (i in 1:n)
			ssm$Yt[i] <- rnorm(1,lambda[i],sqrt(ssm$Vt(i,Xt[i],sqrtphi)))
	}

	if (ssm$fam=="Poisson" & ssm$link=="log") ssm$Yt <- rpois(n,exp(lambda))

	if (ssm$fam=="Poisson" & ssm$link=="identity") ssm$Yt <- rpois(n,lambda)

	if (ssm$fam=="binomial" & ssm$link=="identity") ssm$Yt <- rbinom(n,ssm$nt,lambda/ssm$nt)

	if (ssm$fam=="binomial" & ssm$link=="logit") ssm$Yt <- rbinom(n,ssm$nt,exp(lambda)/(1+exp(lambda)))

	if (ssm$fam=="binomial" & ssm$link=="probit") ssm$Yt <- rbinom(n,ssm$nt,pnorm(lambda,mean=0,sd=1))


}
if (ssm$fam=="multinomial")
{
	# number of categories - 1
	k <- ncol(ssm$Ft(1,Xt[1,],phi))
	
	p <- length(ssm$m0)
	theta <- matrix(NA, ncol=p, nrow=n+1)	
	lambda <- matrix(NA, ncol=k,nrow=n+1)
	theta[1,] <- ssm$m0
	for(i in 2:(n+1))
	{
		theta[i,] <- ssm$Gt(i,Xt[i,],phi)%*%theta[i-1,] + t(rmvnorm(1, mean=rep(0,p), ssm$Wt(i,Xt[i,],ssm$phi)))
		lambda[i,] <- t(ssm$Ft(i,Xt[i,],phi))%*%theta[i,]
	}
	theta <- theta[-1,]
	lambda <- lambda[-1,]
	
	if (ssm$link=="canonical")
	{
		# response function
		h <- function(etas)	exp(etas)/(1+sum(exp(etas)))
		# extract cell probabilities
		pis <- matrix(NA, ncol=k,nrow=n)
		for (i in 1:n) 
			pis[i,] <- h(lambda[i,])
	}
	if (ssm$link=="pom")
	{
		# response function
		h <- function(etas)	exp(etas)/(1+exp(etas))

		# EXTRACT PIES
		pies <- function(gam) gam-c(0,gam[-length(gam)])

		# extract cell probabilities
		pis <- matrix(NA, ncol=k,nrow=n)
		for (i in 1:n)
			pis[i,] <- pies(h(lambda[i,]))	
	}

	# generate multinomials
	ssm$Yt <- matrix(NA, ncol=k+1, nrow=n)
	for (i in 1:n)
		ssm$Yt[i,] <- rmultinom(ssm$nt[i],c(pis[i,], 1-sum(pis[i,])))	
	
}


# Wrap up

#if (is.na(ssm$Yt[1])) 
#	{
#	cat(paste("WARNING!\n Simulating",ssm$fam,"observations using",ssm$link,"link not supported!\n"))
#	break
#}

ssm$simulated <- list(theta=theta, lambda=lambda)
ssm
}

#########################################################
# SUMMARY FUNCTION FOR A SSM OBJECT
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

summary.ssm <- function(ssm)
{
	cat("----------------------------------------------------\n")
	cat("       Summary of specified state space model       \n")
	cat("----------------------------------------------------\n")
	cat("\n")
	cat(paste("Observation are assumed",ssm$fam))
	if (ssm$fam=="multinomial") 
	{
		k <- length(t(ssm$Ft(1,Xt[1,],phi))%*%ssm$Gt(1,Xt[1,],phi)%*%ssm$m0)+1
		cat(paste("with",k,"categories\n"))
		if (is.na(ssm$Yt[1])) cat("\bNo observations supplied")
		else cat(paste("\nNumber of observations:",nrow(ssm$Yt),"\n"))
	}
	else
	{	if (is.na(ssm$Yt[1])) cat("\nNo observations supplied")
		else cat(paste("\nNumber of observations:",length(ssm$Yt),"\n"))
	}
	cat(paste("\nLink function:",ssm$link,"\n\n"))
	if (!is.na(ssm$simulated$theta[1])) cat("NOTE:Data has been simulated.\n\n")
	cat("----------------------------------------------------\n")
	cat("Observation equation:\n\n")
	cat("Design matrix Ft:\n")
	print(ssm$Ft)
	cat("\n")
	if (ssm$fam=="Gaussian")
	{
	cat("Observation variance Vt:\n")
	print(ssm$Vt)
	cat("\n")
	} 
	cat("----------------------------------------------------\n")
	cat("System equation:\n\n")
	cat(paste("Dimension of latent process:",length(ssm$m0),"\n\n"))
	cat("Evolution transfer matrix Gt:\n")
	print(ssm$Gt)
	cat("\nEvolution variance Wt:\n")
	print(ssm$Wt)
	cat("\n")
	cat("----------------------------------------------------\n")
	cat("Covariates:\n\n")
	if (is.na(ssm$Xt)) cat("No covariates specified\n")
	else cat("Dimension of covariates:",dim(Xt),"\n")
	cat("----------------------------------------------------\n")
	cat("Parameters:\n\n")
	if (is.na(ssm$phi[1])) cat("No parameters specified.\n")
	else cat("phi:",ssm$phi,"\n")
	cat("----------------------------------------------------\n")
	cat("Priors:\n\n")
	cat("Prior mean, m0:",ssm$m0,"\n")
	cat("Prior variance, C0:\n")
	print(ssm$C0)
	cat("\n")
	cat("----------------------------------------------------\n")
	cat("Status:\n\n")
	if (is.na(ssm$filtered$mt[1])) cat("Filter has not been applied.\n")
		else cat("Filter has been applied.\n")
	if (is.na(ssm$smoothed$m.tilde[1])) cat("Smoother has not been applied.\n\n")
		else cat("Smoother has been applied.\n\n")
	cat("----------------------------------------------------\n")
	cat("")
}


#####################################################
# KALMAN FILTER FOR BINOMIAL DATA 
# AS IN FAHRMEIR AND TUTZ
# ASSUMING IDENTITY LINK AND MULTIVARIATE LATENT
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.bin.id <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	m.start <- ssm$smoothed$m.tilde

# Dimensions
n <- length(Yt)
latent.dim <- dim(C0)[1]

# initialize
mt <- matrix(NA, ncol=latent.dim, nrow=n)
Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
llh <- -n*log(2*pi)/2

# Filter first observation 
G <- Gt(1,Xt[1,],phi)
FF <- Ft(1,Xt[1,],phi)

at <- G%*%m0
Rt[,,1] <- G%*%C0%*%t(G) + Wt(1,Xt[1,],phi)

if (is.na(m.start[1]))
{
 	lambda.ping <- FF%*%at
}
else
{
	lambda.ping <- FF%*%m.start[1,]
}
ft <- t(FF)%*%at
Qt <- t(FF)%*%Rt[,,1]%*%FF + lambda.ping*(1-lambda.ping/nt[1])

At <- as.vector(Rt[,,1]%*%FF)/Qt
et <- Yt[1]-ft
	
mt[1,] <- at + At*et
Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)

llh <- llh -0.5*(log(Qt)+et^2/Qt)

for(i in 2:n)
{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)

	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi) 
	
	if (is.na(m.start[i]))
	{
		lambda.ping <- FF%*%at
	}
	else
	{
		lambda.ping <- FF%*%m.start[i,]
	}

	if (nt[i]==0)
	{
		mt[i,] <- at 
		Ct[,,i] <- Rt[,,i]
	}
	else
	{	
	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,i]%*%FF + lambda.ping*(1-lambda.ping/nt[i])
	
	At <- as.vector(Rt[,,i]%*%FF)/Qt
	et <- Yt[i]-ft
	
	mt[i,] <- at + At*et
	Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)

	llh <- llh -0.5*(log(Qt)+et^2/Qt)
	}
} 


ssm$filtered<- list(mt=mt,Ct=Ct,Rt=Rt, llh=llh)
ssm
}



#####################################################
# KALMAN FILTER FOR BINOMIAL DATA 
# AS IN FAHRMEIR AND TUTZ
# ASSUMING LOGIT LINK AND MULTIVARIATE LATENT
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.bin.logit <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	m.start <- ssm$smoothed$m.tilde

# Dimensions
n <- length(Yt)
latent.dim <- dim(C0)[1]

# initialize
mt <- matrix(NA, ncol=latent.dim, nrow=n)
Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
llh <- -n*log(2*pi)/2
expit <- function(x) {exp(x)/(1+exp(x))}

# Filter first observation 
G <- Gt(1,Xt[1,],phi)
FF <- Ft(1,Xt[1,],phi)

at <- G%*%m0
Rt[,,1] <- G%*%C0%*%t(G) + Wt(1,Xt[1,],phi)

if (is.na(m.start[1]))
{
 	lambda.ping <- FF%*%at
}
else
{
	lambda.ping <- FF%*%m.start[1,]
}

H.ping <- ((1+exp(lambda.ping))^2)/(nt[1]*exp(lambda.ping))
Y.ping <- lambda.ping + H.ping*Yt[1] - (1+exp(lambda.ping))
	
ft <- t(FF)%*%at
Qt <- t(FF)%*%Rt[,,1]%*%FF + H.ping

At <- as.vector(Rt[,,1]%*%FF)/Qt
et <- (Y.ping-ft)

mt[1,] <- at + At*et
Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)

llh <- llh -0.5*(log(Qt)+et^2/Qt)

for(i in 2:n)
{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)

	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi) 
	
	if (is.na(m.start[i]))
	{
		lambda.ping <- FF%*%at
	}
	else
	{
		lambda.ping <- FF%*%m.start[i,]
	}

	if (nt[i]==0)
	{
		mt[i,] <- at 
		Ct[,,i] <- Rt[,,i]
	}
	else
	{	
		H.ping <- ((1+exp(lambda.ping))^2)/(nt[i]*exp(lambda.ping))
		Y.ping <- lambda.ping + H.ping*Yt[i] - (1+exp(lambda.ping))
	
		ft <- t(FF)%*%at
		Qt <- t(FF)%*%Rt[,,i]%*%FF + H.ping

		At <- as.vector(Rt[,,i]%*%FF)/Qt
		et <- (Y.ping-ft)

		mt[i,] <- at + At*et
		Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)
		
		llh <- llh -0.5*(log(Qt)+et^2/Qt)
	}
	} 
ssm$filtered <- list(mt=mt,Ct=Ct,Rt=Rt,llh=llh)
ssm
}


#####################################################
# KALMAN FILTER FOR BINOMIAL DATA 
# AS IN FAHRMEIR AND TUTZ
# ASSUMING PROBIT LINK AND MULTIVARIATE LATENT
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.bin.probit <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	m.start <- ssm$smoothed$m.tilde

# Dimensions
n <- length(Yt) 
latent.dim <- dim(C0)[1]

# initialize
mt <- matrix(NA, ncol=latent.dim, nrow=n)
Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
llh <- -n*log(2*pi)/2
expit <- function(x) {exp(x)/(1+exp(x))}

# Filter first observation 
G <- Gt(1,Xt[1,],phi)
FF <- Ft(1,Xt[1,],phi)

at <- G%*%m0
Rt[,,1] <- G%*%C0%*%t(G) + Wt(1,Xt[1,],phi)

if (is.na(m.start[1]))
{
 	lambda.ping <- FF%*%at
}
else
{
	lambda.ping <- FF%*%m.start[1,]
}

H.ping <- nt[1]*exp(lambda.ping)/((1+exp(lambda.ping))^2)/(dnorm(lambda.ping)^2*nt[1]^2)
Y.ping <- lambda.ping + (Yt[1] - nt[1]*pnorm(lambda.ping))/(nt[1]*dnorm(lambda.ping))#*H.ping

ft <- t(FF)%*%at
Qt <- t(FF)%*%Rt[,,1]%*%FF + H.ping

At <- as.vector(Rt[,,1]%*%FF)/Qt
et <- (Y.ping-ft)

mt[1,] <- at + At*et
Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)

llh <- llh -0.5*(log(Qt)+et^2/Qt)

for(i in 2:n)
{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)

	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi) 
	
	if (is.na(m.start[i]))
	{
		lambda.ping <- FF%*%at
	}
	else
	{
		lambda.ping <- FF%*%m.start[i,]
	}

	if (nt[i]==0)
	{
		mt[i,] <- at 
		Ct[,,i] <- Rt[,,i]
	}
	else
	{	
		H.ping <- nt[i]*exp(lambda.ping)/((1+exp(lambda.ping))^2)/(dnorm(lambda.ping)^2*nt[i]^2)
		Y.ping <- lambda.ping + (Yt[i] - nt[i]*pnorm(lambda.ping))/(nt[i]*dnorm(lambda.ping))#H.ping
	
		ft <- t(FF)%*%at
		Qt <- t(FF)%*%Rt[,,i]%*%FF + H.ping

		At <- as.vector(Rt[,,i]%*%FF)/Qt
		et <- (Y.ping-ft)

		mt[i,] <- at + At*et
		Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)

		llh <- llh -0.5*(log(Qt)+et^2/Qt)
	}
	} 
ssm$filtered <- list(mt=mt,Ct=Ct,Rt=Rt,llh=llh)
ssm
}


#####################################################
# KALMAN FILTER FOR MULTINOMIAL DATA 
#
# ASSUMING CANONICAL LINK AND MULTIVARIATE LATENT
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.mul.can <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	m.start <- ssm$smoothed$m.tilde

	# Dimensions
	n <- nrow(Yt)
	d <- ncol(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
	h <- function(etas) exp(etas)/(1+sum(exp(etas)))
	llh <- -n*d*log(2*pi)/2

	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G) + Wt(1,Xt[1,],phi)

	if (is.na(m.start[1]))
	{
 	lambda.ping <- t(FF)%*%at
	}
	else
	{
	lambda.ping <- t(FF)%*%m.start[1,]
	}
	pis <- h(lambda.ping)
	Sigma <- nt[1]*(diag(c(pis))-pis%*%t(pis))
	H.ping <- solve(Sigma)

	Y.ping <- lambda.ping + H.ping%*%(Yt[1,-d] - nt[1]*pis)
	
	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,1]%*%FF + H.ping

	At <- Rt[,,1]%*%FF%*%solve(Qt)
	et <- (Y.ping-ft)

	mt[1,] <- at + At%*%et
	Ct[,,1] <- Rt[,,1] - At%*%Qt%*%t(At)

	llh <- llh -0.5*(log(det(Qt))+t(et)%*%solve(Qt)%*%(et))

	for(i in 2:n)
	{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)

	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi) 
	
	if (is.na(m.start[i]))
	{
		lambda.ping <- t(FF)%*%at
	}
	else
	{
		lambda.ping <- t(FF)%*%m.start[i,]
	}

	if (nt[i]==0)
	{
		mt[i,] <- at 
		Ct[,,i] <- Rt[,,i]
	}
	else
	{	
	pis <- h(lambda.ping)

	Sigma <- nt[i]*(diag(c(pis))-pis%*%t(pis))
	H.ping <- solve(Sigma)

	Y.ping <- lambda.ping + H.ping%*%(Yt[i,-d] - nt[i]*pis)

	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,i]%*%FF + H.ping

	At <- Rt[,,i]%*%FF%*%solve(Qt)
	et <- (Y.ping-ft)

	mt[i,] <- at + At%*%et
	Ct[,,i] <- Rt[,,i] - At%*%Qt%*%t(At)
		
	llh <- llh -0.5*(log(det(Qt))+t(et)%*%solve(Qt)%*%et)
	}
	} 
	ssm$filtered <- list(mt = mt, Ct = Ct, Rt = Rt, llh = llh)
	ssm
}

#####################################################
# KALMAN FILTER FOR MULTINOMIAL DATA 
#
# ASSUMING POM LINK AND MULTIVARIATE LATENT
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.mul.pom <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	m.start <- ssm$smoothed$m.tilde

	# Dimensions
	n <- nrow(Yt)
	d <- ncol(Yt) - 1
	latent.dim <- length(m0)

	# initialize
	mt <- matrix(NA, ncol = latent.dim, nrow = n)
	Ct <- array(NA, dim = c(latent.dim, latent.dim, n))
	Rt <- array(NA, dim = c(latent.dim, latent.dim, n))
	llh <-  - n * d * log(2 * pi)/2

	# Cumsum observations
	Zt <- t(apply(Yt, 1, cumsum))

	# Misc functions
	dexpit <<- function(x) expit(x) * (1 - expit(x))
	expit <<- function(x) exp(x)/(1 + exp(x))

	# RESPONS FUNCTION
	h <- function(lambda, n) n * expit(lambda)
	
	# DERIVATIVE OF RESPONSE FUNCTION
	dh <- function(lambda, n)	n * diag(dexpit(lambda))

	# EXTRACT GAMMAS FROM LAMBDA
	gammas <- function(lambda, n) exp(lambda)/(1 + exp(lambda))

	# EXTRACT PIS
	pis <- function(gam)	gam - c(0, gam[ - length(gam)])

	# VARIANCE - without zero row and column
	Gam <- function(gam, n)
	{
		# Create lower triangle
		d <- length(gam)
		result <- matrix(0, d, d)
		for(i in 1:(d - 1))
			result[((i + 1):d), i] <- gam[i] * (1 - gam[ - c(1:i)])
		result <- result + t(result)
		diag(result) <- gam * (1 - gam)
		result <- n * result
		result
	}

	# Filter first observation 
	G <- Gt(1, Xt[1,  ], phi)
	FF <- Ft(1, Xt[1,  ], phi)

	at <- G %*% m0
	Rt[,  , 1] <- G %*% C0 %*% t(G) + Wt(1, Xt[1,  ], phi)

	if(is.na(m.start[1])) {
		lambda.ping <- t(FF) %*% at
	}
	else {
		lambda.ping <- t(FF) %*% m.start[1,  ]
	}
	gam <- gammas(lambda.ping, nt[1])
	h.ping <- dh(c(lambda.ping), nt[1])
	h.ping.inv <- diag(1/diag(h.ping))

	Vt <- h.ping.inv %*% Gam(gam, nt[1]) %*% h.ping.inv 
	Y.ping <- lambda.ping + h.ping.inv %*% (Zt[1,  - (d + 1)] - h(lambda.ping, nt[1]))

	ft <- t(FF) %*% at
	Qt <- t(FF) %*% Rt[,  , 1] %*% FF + Vt

	At <- Rt[,  , 1] %*% FF %*% solve(Qt)
	et <- (Y.ping - ft)

	mt[1,  ] <- at + At %*% et
	Ct[,  , 1] <- Rt[,  , 1] - At %*% Qt %*% t(At)
	llh <- llh -0.5*(log(det(Qt))+t(et)%*%solve(Qt)%*%(et))

	for(i in 2:n) {
		G <- Gt(i, Xt[i,  ], phi)
		FF <- Ft(i, Xt[i,  ], phi)

		at <- G %*% mt[i - 1,  ]
		Rt[,  , i] <- G %*% Ct[,  , i - 1] %*% t(G) + Wt(i, Xt[i,  ], phi)

		if(is.na(m.start[i])) {
			lambda.ping <- t(FF) %*% at
		}
		else {
			lambda.ping <- t(FF) %*% m.start[i,  ]
		}

		if(nt[i] == 0) {
			mt[i,  ] <- at
			Ct[,  , i] <- Rt[,  , i]
		}
		else {
			gam <- gammas(lambda.ping, nt[i])
			h.ping <- dh(c(lambda.ping), nt[i])
			h.ping.inv <- diag(1/diag(h.ping))

			Vt <- h.ping.inv %*% Gam(gam, nt[i]) %*% h.ping.inv 
			Y.ping <- lambda.ping + h.ping.inv %*% (Zt[i,  - (d + 1)] - h(lambda.ping,nt[i]))

			ft <- t(FF) %*% at
			Qt <- t(FF) %*% Rt[,  , i] %*% FF + Vt

			At <- Rt[,  , i] %*% FF %*% solve(Qt)
			et <- (Y.ping - ft)

			mt[i,  ] <- at + At %*% et
			Ct[,  , i] <- Rt[,  , i] - At %*% Qt %*% t(At)
			llh <- llh -0.5*(log(det(Qt))+t(et)%*%solve(Qt)%*%(et))
		}
	}
	ssm$filtered <- list(mt = mt, Ct = Ct, Rt = Rt, llh = llh)
	ssm
}

#####################################################
# KALMAN FILTER FOR POISSON DATA 
# AS IN FAHRMEIR AND TUTZ
# ASSUMING LOG LINK AND MULTIVARIATE LATENT
# Used by iterated extended kalman filter/smoother
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.po.log <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	m.start <- ssm$smoothed$m.tilde

# Dimensions
n <- length(Yt)
latent.dim <- dim(C0)[1]

# initialize
mt <- matrix(NA, ncol=latent.dim, nrow=n)
Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
llh <- -n*log(2*pi)/2

# Filter first observation 
G <- Gt(1,Xt[1,],phi)
FF <- Ft(1,Xt[1,],phi)

at <- G%*%m0
Rt[,,1] <- G%*%C0%*%t(G) + Wt(1, Xt[1,], phi)

if (is.na(m.start[1]))
	{
	lambda.ping <- FF%*%at
	}
else
	{
	lambda.ping <- FF%*%m.start[1,]
	}
H.ping <- exp(-lambda.ping)	
Y.ping <- lambda.ping + H.ping*Yt[1] - 1

ft <- t(FF)%*%at
Qt <- t(FF)%*%Rt[,,1]%*%FF + H.ping

At <- as.vector(Rt[,,1]%*%FF)/Qt
et <- Y.ping-ft

mt[1,] <- at + At*et
Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)

llh <- llh-0.5*(log(Qt)-et^2/Qt)

for(i in 2:n)
{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)
	
	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i, Xt[i,], phi) 

	if (is.na(m.start[1]))
		{
		lambda.ping <- FF%*%at
		}
	else
		{
		lambda.ping <- FF%*%m.start[i,]
		}
	
	H.ping <- exp(-lambda.ping)
	Y.ping <- lambda.ping + H.ping*Yt[i] - 1
	
	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,i]%*%FF + H.ping
	
	At <- as.vector(Rt[,,i]%*%FF)/Qt
	et <- Y.ping - ft
	
	mt[i,] <- at + At*et
	Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)

	llh <- llh-0.5*(log(Qt)+et^2/Qt)
} 

ssm$filtered <- list(mt=mt,Ct=Ct,Rt=Rt, llh=llh)
ssm
}
#####################################################
# KALMAN FILTER FOR POISSON DATA 
# AS IN FAHRMEIR AND TUTZ
# ASSUMING IDENTITY LINK AND MULTIVARIATE LATENT
# Used by iterated extended kalman filter/smoother
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kal.fil.po.id <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	m.start <- ssm$smoothed$m.tilde

	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
	llh <- -n*log(2*pi)/2
	
	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	if (is.na(m.start[1]))
	{
 	lambda.ping <- FF%*%at
	}
	else
	{
	lambda.ping <- FF%*%m.start[1,]
	}
	Rt[,,1] <- G%*%C0%*%t(G) + Wt(1, Xt[1,], phi)

	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,1]%*%FF + lambda.ping

	At <- as.vector(Rt[,,1]%*%FF)/Qt
	
	mt[1,] <- at + At*(Yt[1]-ft)
	Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)
	
	for(i in 2:n)
	{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)
	
	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i, Xt[i,], phi)  

	if (is.na(m.start[1]))
	{
		lambda.ping <- FF%*%at
	}
	else
	{
		lambda.ping <- FF%*%m.start[i,]
	}
	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,i]%*%FF + lambda.ping
	
	At <- as.vector(Rt[,,i]%*%FF)/Qt

	mt[i,] <- at + At*(Yt[i]-ft)
	Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)
	} 

ssm$filtered <- list(mt=mt,Ct=Ct,Rt=Rt, llh=llh)
ssm
}

generate <-F
if (generate)
{
Ft <- function(i,x,phi)
{
	result <- c(1,0)
	result
}

Gt <- function(i,x,phi)
{
	result <- matrix(c(0,1,0,1), ncol=2, byrow=T)
	result
}
Gt(1)

Wt <- function(i,x,phi)
{
	result <- matrix(c(phi,0,0,0),ncol=2)
	result
}

#simulate data
m0 <- c(0,5)
n <- 100
theta <- matrix(NA, ncol=2,nrow=n)
Yt <- rep(NA,n)
theta[1,] <- m0
for( i in 2:n)
{ 
	theta[i,] <- Gt(i)%*%theta[i-1,] + 	t(rmvnorm(1,mean=c(0,0), cov=Wt(1,1,0.01)))
	Yt[i] <- rpois(1,t(Ft(i))%*%theta[i,])
}


tsplot((theta[-1,1]), ylim=range(Yt[-1]))
tspoints(Yt[-1], pch=16)

bk.fil <- kal.fil.po.id(ss)
ss <- ssm(Ft=Ft,
					Gt=Gt,
					Wt=Wt,
					Xt=NA, 
					phi=0.01,
					m0=c(0,5),
					C0=diag(rep(3,2)),
					Yt=Yt[-1], m.start=NA,fam="Poisson",link="identity")


bk.smo <- kalman.smoother(bk.fil)
bk.smo <- bk.smo$smoothed	
bk.fil <- bk.smo$filtered
tsplot(bk.fil$mt[,1])
tsplot(bk.fil$mt[,2])
tsplot(bk.smo$m.tilde[,1], ylim=c(range(theta[-1,1])))
tsplot(bk.fil$mt[,2])

tslines((theta[-1,1]), ylim=range(Yt[-1]), lty=2)

tsplot(bk.smo$m.tilde[,1])
tsplot(theta[-1,1], lty=2)
tspoints(Yt[-1], pch=16)
}

#####################################################
# KALMAN FILTER FOR POISSON DATA 
# AS IN WEST, HARRISON AND MIGON, JASA 1985
# ASSUMING LOG LINK, GAMMA PRIOR AND MULTIVARIATE LATENT
# For use with vandriver example
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

ckal.fil.po.gam.log <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt

	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
	f1 <- function(alpha,qt) {(qt-trigamma(alpha))^2}

	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G)+Wt(1,Xt[1,],phi)

	ft <- FF%*%at
	qt <- t(FF)%*%Rt[,,1]%*%FF

	alpha <- mean(optimize(f1, interval=c(0,100000), qt=qt)$interval)
	beta <- exp(digamma(alpha)-ft)

	ft.star <- digamma(alpha+Yt[1]) - log(beta+1)
	qt.star <- trigamma(alpha+Yt[1])

	mt[1,] <- at + (Rt[,,1]%*%FF)*as.vector((ft.star-ft)/qt)
	Ct[,,1] <- Rt[,,1] - Rt[,,1]%*%FF%*%t(FF)%*%Rt[,,1]*as.vector((1-qt.star/qt)/qt)
	
	for(i in 2:n)
	{ 
		G <- Gt(i, Xt[i,], phi)
		FF <- Ft(i, Xt[i,], phi)

		at <- G%*%mt[i-1,]
		Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi)
		ft <- FF%*%at
		qt <- t(FF)%*%Rt[,,i]%*%FF

		alpha <- mean(optimize(f1, interval=c(0,100000), qt=qt)$interval)
		beta <- exp(digamma(alpha)-ft)

		if (is.na(Yt[i])) 
			{
			mt[i,] <- at
			Ct[,,i] <- Rt[,,i-1]
			}
		else
			{
			ft.star <- digamma(alpha+Yt[i]) - log(beta+1)
			qt.star <- trigamma(alpha+Yt[i])

			mt[i,] <- at+Rt[,,i]%*%FF*as.vector((ft.star-ft)/qt)
			Ct[,,i] <- Rt[,,i]-Rt[,,i]%*%FF%*%t(FF)%*%Rt[,,i]*as.vector((1-qt.star/qt)/qt)
			}
	}

	ssm$filtered <- list(mt=mt,Ct=Ct, Rt=Rt)
	ssm
}

#####################################################
# KALMAN FILTER FOR POISSON DATA 
# AS IN WEST, HARRISON AND MIGON, JASA 1985
# ASSUMING IDENTITY LINK, GAMMA PRIOR AND MULTIVARIATE LATENT
# For use with vandriver example
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

ckal.fil.po.gam.id <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt

	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
	f1 <- function(alpha,qt) {(qt-trigamma(alpha))^2}

	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G)+Wt(1,Xt[1,],phi)

	ft <- FF%*%at
	qt <- t(FF)%*%Rt[,,1]%*%FF

	alpha <- ft^2/qt
	beta <- ft/qt
	ft.star <- (alpha+Yt[1])/(beta+1)
	qt.star <- (alpha+Yt[1])/((beta+1)^2)

	mt[1,] <- at + (Rt[,,1]%*%FF)*as.vector((ft.star-ft)/qt)
	Ct[,,1] <- Rt[,,1] - Rt[,,1]%*%FF%*%t(FF)%*%Rt[,,1]*as.vector((1-qt.star/qt)/qt)
	
	for(i in 2:n)
	{ 
		G <- Gt(i, Xt[i,], phi)
		FF <- Ft(i, Xt[i,], phi)

		at <- G%*%mt[i-1,]
		Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi)
		ft <- FF%*%at
		qt <- t(FF)%*%Rt[,,i]%*%FF

		alpha <- ft^2/qt
		beta <- ft/qt

		if (is.na(Yt[i])) 
			{
			mt[i,] <- at
			Ct[,,i] <- Rt[,,i-1]
			}
		else
			{
			ft.star <- (alpha+Yt[i])/(beta+1)
			qt.star <- (alpha+Yt[i])/((beta+1)^2)

			mt[i,] <- at+Rt[,,i]%*%FF*as.vector((ft.star-ft)/qt)
			Ct[,,i] <- Rt[,,i]-Rt[,,i]%*%FF%*%t(FF)%*%Rt[,,i]*as.vector((1-qt.star/qt)/qt)
			}
	}

	ssm$filtered <- list(mt=mt,Ct=Ct, Rt=Rt)
	ssm
}

#####################################################
# KALMAN FILTER FOR BINOMIAL DATA 
# AS IN WEST, HARRISON AND MIGON, JASA 1985
# ASSUMING ID LINK, BETA PRIOR AND MULTIVARIATE LATENT
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

ckal.fil.bin.beta.id <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt
	
	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))

	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G)+Wt(1,Xt[1,],phi)

	ft <- FF%*%at
	qt <- t(FF)%*%Rt[,,1]%*%FF

	alpha <- ft*(ft*(nt[1]-ft)/qt-1)/nt[1]
	beta <- (1-ft/nt[1])*(ft*(nt[1]-ft)/qt-1)
	
	ft.star <- nt[1]*(alpha+Yt[1])/(alpha+beta+nt[1])
	qt.star <- ft.star*(nt[1]-ft.star)/(alpha+beta+nt[1]+1)
	
	mt[1,] <- at + (Rt[,,1]%*%FF)*as.vector((ft.star-ft)/qt)
	Ct[,,1] <- Rt[,,1] - Rt[,,1]%*%FF%*%t(FF)%*%Rt[,,1]*as.vector((1-qt.star/qt)/qt)
	
	for(i in 2:n)
	{ 
		G <- Gt(i, Xt[i,], phi)
		FF <- Ft(i, Xt[i,], phi)

		at <- G%*%mt[i-1,]
		Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi)
		ft <- FF%*%at
		qt <- t(FF)%*%Rt[,,i]%*%FF

		alpha <- ft*(ft*(nt[i]-ft)/qt-1)/nt[i]
		beta <- (1-ft/nt[i])*(ft*(nt[i]-ft)/qt-1)

		if (is.na(Yt[i])) 
			{
			mt[i,] <- at
			Ct[,,i] <- Rt[,,i-1]
			}			
		else
			{
			ft.star <- nt[i]*(alpha+Yt[i])/(alpha+beta+nt[i])
			qt.star <- ft.star*(nt[i]-ft.star)/(alpha+beta+nt[i]+1)

			mt[i,] <- at+Rt[,,i]%*%FF*as.vector((ft.star-ft)/qt)
			Ct[,,i] <- Rt[,,i]-Rt[,,i]%*%FF%*%t(FF)%*%Rt[,,i]*as.vector((1-qt.star/qt)/qt)
			}
	}

	ssm$filtered <- list(mt=mt,Ct=Ct,	 Rt=Rt)
	ssm
}

#####################################################
# KALMAN FILTER FOR BINOMIAL DATA 
# AS IN WEST, HARRISON AND MIGON, JASA 1985
# ASSUMING LOGIT LINK, BETA PRIOR AND MULTIVARIATE LATENT
#
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

ckal.fil.bin.beta.logit <- function(ssm)
{
	# require
	# ssm: object of class ssm

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt
	nt <- ssm$nt

	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))

	# Function used to determine prior
	f <- function(vec) 
		{
		rt <- vec[1]	
		st <- vec[2]
		ff <- digamma(rt)-digamma(st)
		qq <- trigamma(rt)+trigamma(st)
		result <- sqrt((ft-ff)^2+(qt-qq)^2)
		result	
	}

	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G)+Wt(1,Xt[1,],phi)

	ft <<- FF%*%at
	qt <<- t(FF)%*%Rt[,,1]%*%FF

	par <- nlm(p=c((1+exp(ft))/qt,(1+exp(-ft))/qt), 
		             f=f)$estimate

		alpha <- par[1]
		beta <- par[2]	
	

	ft.star <- digamma(alpha+Yt[1])-digamma(beta+nt[1]-Yt[1])
	qt.star <- trigamma(alpha+Yt[1])+trigamma(beta+nt[1]-Yt[1])

	mt[1,] <- at + (Rt[,,1]%*%FF)*as.vector((ft.star-ft)/qt)
	Ct[,,1] <- Rt[,,1] - Rt[,,1]%*%FF%*%t(FF)%*%Rt[,,1]*as.vector((1-qt.star/qt)/qt)

	for(i in 2:n)
	{ 
		G <- Gt(i, Xt[i,], phi)
		FF <- Ft(i, Xt[i,], phi)

		at <- G%*%mt[i-1,]
		Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi)
		

		ft <<- FF%*%at
		qt <<- t(FF)%*%Rt[,,i]%*%FF

		par <- nlm(p=c((1+exp(ft))/qt,(1+exp(-ft))/qt), 
		             f=f)$estimate

		alpha <- par[1]
		beta <- par[2]	
	
		if (is.na(Yt[i])) 
			{
			mt[i,] <- at
			Ct[,,i] <- Rt[,,i-1]
			}			
		else
			{
			ft.star <- digamma(alpha+Yt[i])-digamma(beta+nt[i]-Yt[i])
			qt.star <- trigamma(alpha+Yt[i])+trigamma(beta+nt[i]-Yt[i])

			mt[i,] <- at+Rt[,,i]%*%FF*as.vector((ft.star-ft)/qt)
			Ct[,,i] <- Rt[,,i]-Rt[,,i]%*%FF%*%t(FF)%*%Rt[,,i]*as.vector((1-qt.star/qt)/qt)
			}
	}

	ssm$filtered <- list(mt=mt,Ct=Ct,	 Rt=Rt)
	ssm
}

#####################################################
# KALMAN FILTER FOR GAUSSIAN DATA 
# AS IN WEST AND HARRISON
# For use with iterated extended kalman filter, ieks
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kalman.filter <- function(ssm)
{
	# require
	# ssm: object of class state space model

	# Returns an updated ssm object
	# mt: Filtered states
	# Ct: Filtered MSE
	# Rt: Prior variances (needed in the Kalman smoother)
	# llh: Log-likelihood

	# Extract information from ssm object
	Ft <- ssm$Ft
	Gt <- ssm$Gt
	Wt <- ssm$Wt
	Vt <- ssm$Vt
	Xt <- ssm$Xt
	phi <- ssm$phi
	m0 <- ssm$m0
	C0 <- ssm$C0
	Yt <- ssm$Yt

	# Dimensions
	n <- length(Yt)
	latent.dim <- dim(C0)[1]

	# initialize
	mt <- matrix(NA, ncol=latent.dim, nrow=n)
	Ct <- array(NA, dim=c(latent.dim,latent.dim,n))
	Rt <- array(NA, dim=c(latent.dim,latent.dim,n))
	llh <- -n*log(2*pi)/2
	
	# Filter first observation 
	G <- Gt(1,Xt[1,],phi)
	FF <- Ft(1,Xt[1,],phi)

	at <- G%*%m0
	Rt[,,1] <- G%*%C0%*%t(G) + Wt(1,Xt[1,],phi)		

	ft <- t(FF)%*%at
	Qt <- (t(FF)%*%Rt[,,1]%*%FF) + Vt(1,Xt[1,],phi)

	At <- as.vector(Rt[,,1]%*%FF)/Qt
	et <- Yt[1]-ft
	
	mt[1,] <- at + At*et
	Ct[,,1] <- Rt[,,1] - At%*%t(At)*as.vector(Qt)
	llh <- llh -0.5*(log(Qt)+et^2/Qt)

	for(i in 2:n)
	{
	G <- Gt(i,Xt[i,],phi)
	FF <- Ft(i,Xt[i,],phi)

	at <- G%*%mt[i-1,]
	Rt[,,i] <- G%*%Ct[,,i-1]%*%t(G) + Wt(i,Xt[i,],phi)
	
	if (is.na(Yt[i]))
	{
		mt[i,] <- at
		Ct[,,i] <- Rt[,,i]
	}
	else
	{
	ft <- t(FF)%*%at
	Qt <- t(FF)%*%Rt[,,i]%*%FF + Vt(i,Xt[i],phi)

	At <- as.vector(Rt[,,i]%*%FF)/Qt
	et <- Yt[i]-ft
	
	mt[i,] <- at + At*et
	Ct[,,i] <- Rt[,,i] - At%*%t(At)*as.vector(Qt)
	llh <- llh -0.5*(log(Qt)+et^2/Qt)
	} 
}
ssm$filtered <- list(mt=mt,Ct=Ct,Rt=Rt,llh=llh)
ssm
}


#####################################################
# KALMAN SMOOTHER  
# Used by iterated extended Kalman filter/smoother
# bjarke klein 2003
# bjarke@statdem.sdu.dk
#####################################################

kalman.smoother <- function(ssm) {
	# Requires
	# ssm: object of class state space model
	# requires that the filter has been run!
	
	# Returns an updated ssm object with
	# m.tilde
	# C.tilde
	
	if(is.na(ssm$filtered$mt[1])) {
		cat("You have not run the filter yet - please do that\n")
		break
	}
	
	# Extract information from ssm object 
	Gt <- ssm$Gt
	m.tilde <- ssm$filtered$mt
	C.tilde <- ssm$filtered$Ct
	Rt <- ssm$filtered$Rt
	
	# number of observations
	n <- dim(m.tilde)[1]
	
	# Dimension of latent process 
	latent.dim <- dim(m.tilde)[2]
	
	# dimension of observation
	obs.dim <- 1
	
	# Initialize
	# First smoothed is last filtered
	G <- Gt(n, Xt[n,  ], phi)
	
	for(i in (n - 1):1) {
		Gtp1 <- G
		G <- Gt(i, Xt[i,  ], phi)
		Bt <- C.tilde[,  , i] %*% t(Gtp1) %*% solve(Rt[,  , i + 1])
		m.tilde[i,  ] <- m.tilde[i,  ] + Bt %*% (m.tilde[i + 1,  ] - Gtp1 %*% m.tilde[i,  ])
		C.tilde[,  , i] <- C.tilde[,  , i] + Bt %*% (C.tilde[,  , i + 1] - Rt[,  , i + 1]) %*% t(Bt)
	}
	ssm$smoothed <- list(m.tilde = m.tilde, C.tilde = C.tilde)
	ssm
}


###############################
# Generate single random Multinomial(n,pr)
# source http://www.stat.cmu.edu/~hseltman/files/multinom.q
#################################
rmultinom <- function(n=5, pr=c(0.5,0.5), long=F)
{
  k <- length(pr)
  if (abs(1-sum(pr))>0.000001) stop("rmultinom: parameter pr must be the k probabilities (summing to 1)")

  if(long) {
    y <- runif(n, 0, 1)
    p <- cumsum(pr)
    Seq <- 1:n
    x <- sapply(y, function(y, Seq, p) {Seq[y <= p][1]}, Seq=Seq, p=p)
  } else {
    x <- rep(NA,k)
    p <- pr/c(1,(1-cumsum(pr[1:(k-1)])))
    for (i in 1:(k-1)) {
      if (n==0) {
        x[i] <- 0
        if (i==k-1) x[k] <- 0
        next
      }
      y <- rbinom(1,n,p[i])
      x[i] <- y
      if (i==k-1) x[k] <- n-y
      n <- n-y
    }
  }
  return(x)
}

rmvnorm <- function(n, mean=rep(0, nrow(sigma)),
                      sigma=diag(length(mean))){

  if(nrow(sigma) != ncol(sigma)){
    stop("sigma meanst be a square matrix")
  }

  if(length(mean) != nrow(sigma)){
    stop("mean and sigma have non-conforming size")
  }
  
  sigsvd <- svd(sigma)
  retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  retval
}


dmvnorm <- function(x, mean, sigma){

  if(is.vector(x)){
    x <- matrix(x, ncol=length(x))
  }

  if(missing(mean)){
    mean <- rep(0, length=ncol(x))
  }
  
  if(missing(sigma)){
    sigma <- diag(ncol(x))
  }

  if(ncol(x) != ncol(sigma)){
    stop("x and sigma have non-conforming size")
  }
  
  if(nrow(sigma) != ncol(sigma)){
    stop("sigma meanst be a square matrix")
  }
  if(length(mean) != nrow(sigma)){
    stop("mean and sigma have non-conforming size")
  }

  retval <- exp(-mahalanobis(x, center=mean, cov=sigma)/2)
  det <- prod(eigen(sigma, sym=TRUE)$values)
  retval<- retval / (sqrt(det) * sqrt(2*pi)^ncol(x))

  retval
}
  

vandrivers <- matrix(
c(12,6,12,8,10,13,11,6,10,16,13,14,14,6,8,11,7,13,13,11,11,14,16,14,17,16,15,13,13,15,12,6,9,13,14,15,14,3,12,13,12,8,8,15,8,5,17,14,13,5,8,
5,12,11,13,15,11,11,10,13,8,6,8,14,12,14,13,9,4,13,6,15,12,16,7,12,10,9,9,6,7,13,14,13,14,11,11,10,4,8,9,10,10,5,13,12,10,9,7,5,10,5,
6,8,6,12,15,7,14,4,10,8,7,11,3,5,11,10,10,7,10,11,9,7,8,13,8,5,8,7,12,10,7,4,10,4,8,8,7,10,8,14,8,9,8,6,7,6,5,4,5,10,7,
10,12,7,4,5,6,4,4,8,8,3,7,12,2,7,8,3,2,6,3,7,
6,8,8,4,3,5,5,3,4,3,6,6,7,5,7,7,4,
7,rep(0,169),rep(1,23)),ncol=2,byrow=F)