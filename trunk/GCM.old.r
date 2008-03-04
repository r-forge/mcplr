
gcm.cont.fit <- function(x,y,w,lambda,r,q) {
  gcm.distance.matrix <- function(x,w,r) {
    nc <- ncol(x)
    nr <- nrow(x)
    dis <- matrix(0,ncol=nr,nrow=nr)
    for(i in 1:nc) {
      dis <- dis + w[i]*abs(outer(x[,i],x[,i],"-"))^r
    }
    dis^(1/r)
  }
  gamm <- .2
  lamb <- 1e-2
  sdy <- 1

  dis <- gcm.distance.matrix(x,w=w,r=r)
  sim <- exp(-lambda*dis^q)
  sim[lower.tri(sim,diag=T)] <- 0

  sw <- exp(


tmp <- gcm.distance.matrix(x,w=c(.25,.25,.25,.25),r=2)
tmp <- exp(-lamb*tmp^2)
tmp[lower.tri(tmp,diag=T)] <- 0
sam <- exp(gamm*1:nt)
tmp2 <- sam*tmp
tmp2 <- t(t(tmp2)/colSums(tmp2))
tmp2[is.nan(tmp2)] <- 0

}

gcm.discr.fit <- function(x,y,w,lambda) {


}

gcm <- function(


setClass("gcm",
  representation(
    x="matrix",
    y="vector",
    distance="function",
    similarity="function",
    sampling="function",
    parameters="list"),
  prototype(
    x=matrix(),
    y=vector(),
    distance=function(x,y,w,r) sum(w*abs(x-y)^r)^(1/r),
    similarity=function(d,lambda,q) exp(-lambda*(d^q)),
    sampling=
    parameters=list(),
)

gcm.distance <- function(x,y,w,r) {
    sum(w*abs(x-y)^r)^(1/r)
}


###

gcm.distance.matrix <- function(x,w,r) {
  nc <- ncol(x)
  nr <- nrow(x)
  dis <- matrix(0,ncol=nr,nrow=nr)
  for(i in 1:nc) {
    dis <- dis + w[i]*abs(outer(x[,i],x[,i],"-"))^r
  }
  dis^(1/r)
}

gcm.sampling.matrix <- function(nt,gamma) {
  k <- seq(1,nt,by=1)
  apply(k
}

nt <- 200
nx <- 4
x <- matrix(rnorm(nt*nx,sd=5),ncol=nx)
y <- x%*%c(.1,.2,.3,.4) + rnorm(nt,sd=.1)
r <- seq(-10,10,length=200) # grid to evaluate Lik(r)



gamm <- .2
lamb <- 1e-2
sdy <- 1



tmp <- gcm.distance.matrix(x,w=c(.25,.25,.25,.25),r=2)
tmp <- exp(-lamb*tmp^2)
tmp[lower.tri(tmp,diag=T)] <- 0
sam <- exp(gamm*1:nt)
tmp2 <- sam*tmp
tmp2 <- t(t(tmp2)/colSums(tmp2))
tmp2[is.nan(tmp2)] <- 0

tmp3 <- outer(r,as.vector(y),dnorm,sd=sdy)

lik.grid <- apply(tmp2,2,"%*%",t(tmp3))
#norm <- colSums(lik.grid)
#lik.grid <- t(t(lik.grid)/colSums(lik.grid))
#lik.grid[is.nan(lik.grid)] <- 0

ids <- c(30,50,70,90)
plot(c(min(r),max(r)),c(0,max(lik.grid[,ids])),type="n")
for(i in 1:length(ids)) {
  lines(r,lik.grid[,ids[i]],lty=i)
}




gcm.similarity <- function(d,lamda,q) {
  exp(-lambda*(d^q))
}

gcm.lik <- function(y,x,w,q=1,r=1,lambda=1,gamma=1,sigma=1,nt) {
  for(i in 2:nt) {
    dis <- apply(x[1:(i-1),],1,gcm.distance,y=x[i,],w=w,r=r)
    sim <- gcm.similarity(dis)
    p <- exp(-gamma*(1:(i-1)))


x <- cbind(c(3,4,1,8,19),c(5,77,3,1,12))

x <- matrix(rnorm(800),ncol=4)
tmp <- gcm.distance.matrix(x,w=c(.3,.3,.2,.2),r=2)

y <- 1:12




