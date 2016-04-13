\name{GCM}
\alias{GCM}
\alias{GCM-class}
\alias{GCMlearning-class}
\alias{GCMresponse-class}
\alias{fit,GCM-method}
\alias{getPars,GCM-method}
\alias{predict,GCM-method}
\alias{runm,GCM-method}
\alias{setPars,GCM-method}
\alias{logLik,GCM-method}

\title{Generalized Context Model}
\description{
Constructs a Generalized Context Model.
}
\usage{
GCM(learning,response,parameters=list(w=NULL,lambda=1,r=1,q=1,gamma=NULL),
fixed,data,subset,ntimes=NULL,replicate=TRUE,remove.intercept=FALSE)
}
\arguments{
\item{learning}{an object of class \code{formula} (or one that can be coerced to 
that class): a symbolic description of the learning model to be fitted. For more 
details of model speecification, see \code{lm} or \code{glm}.}
\item{response}{an object of class \code{formula} (or one that can be coerced to 
that class): a symbolic description of the response model to be fitted. Only the
left hand side of the formula will be used. For more details of model 
specification, see \code{lm} or \code{glm}.}
\item{parameters}{an (optional) named list with (starting) values of the 
parameters. If no values are supplied, defaults are used. }
\item{fixed}{a logical vector indicating whether model parameters are fixed or 
free}
\item{data}{data frame containing the variables listed in the formula argument.}
\item{subset}{subset of data to fit model to.}
\item{ntimes}{a vector with number of observations for each replication.}
\item{replicate}{logical to indicate whether model parameters are identical for 
each replication in ntimes.}
\item{remove.intercept}{(logical) should the intercept term be removed from 
the x matrix of the model?}
}
\details{The Generalized Context Model (Nosofsky, 1986) is an exemplar model.
It predicts the value of a criterion at t based on the similarity of a probe cue 
\eqn{x_t}{x[t]} to stored cues \eqn{x_{t-k}}{x[t-k]}, \eqn{k = 1,\ldots,t-1}{k=1,...,t-1}.

The similarity between two cues is computed as
\deqn{s(x_t,x_k) = e^(-\lambda d(x_t,x_k)^{q}}{s(x[t],x[k]) = exp(-lambda * d(x[t],x[k])^(q) )}
where \eqn{d(x_t,x_k)}{d(x[t],x[k])} is the generalized Minkowski distance
\deqn{d(x_t,x_k) = (\sum_{j} w_j | x_{jt} - x_{jk} |^{r})^{\frac{1}{r}}}{d(x[t],x[k]) = (sum[j] w[j] | x[jt] - x[jk] |^(r))^(1/r)}
When making a response to cue \eqn{x_t}{x[t]}, first the overall similarity to each of the previously encountered exemplars of each of the levels of outcome y is determined as
\deqn{S(j) = \sum s}{} 

The GCM can be seen as a mixture model. Each encountered exemplar adds a
new component to the mixture. The mixture proportions are defined by the
similarities. See the package manual for more information.

A more general version of the GCM is implemented by \code{\link{gGCM}}.
}

\value{An object of class \code{GCM}.}

\references{
Nosofsky, R.M. (1986). Attention, Similarity, and the 
Identification-Categorization Relationship. \emph{Journal of Experimental 
Psychology: General}, \emph{115}, 39-57.}
\examples{
## open weather prediction data
data(WPT)
controls <- subset(WPT,id \%in\% paste("C",1:16,sep=""))
## initialize model, use remove.intercept=TRUE so that the x matrix will contain
##   four columns.
mod <- GCM(y~x1+x2+x3+x4,response=r~1,data=controls,
   fix=list(r=TRUE,q=TRUE),ntimes=rep(200,16),
   remove.intercept=TRUE)
\dontrun{
## estimate free parameters
## discount (ignore) first 5 responses in each series
## as these can give deterministic predictions
mod <- fit(mod,discount=5)
}
summary(mod)
}
\keyword{models}

