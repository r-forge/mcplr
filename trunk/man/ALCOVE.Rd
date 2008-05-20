\name{alcove}
\alias{alcove}
\title{Attention Learning COVEring map model}
\description{
Constructs an Attention Learning COVEring map model.
}
\usage{
alcove(formula,level=c("nominal","interval"),distance=gcm.distance("cityblock"),
 similarity=gcm.similarity("exponential"),sampling=gcm.sampling("uniform"),
 parameters=list(),fixed,parStruct,data,subset,ntimes=NULL,replicate=TRUE,
 base=NULL)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to
 that class): a symbolic description of the model to be fitted. For more details
 of model speecification, see \code{lm} or \code{glm}.}
\item{parameters}{an (optional) list with (starting) values of the parameters.
 If no values are supplied, defaults are used. }
\item{humble}{logical. If TRUE, humble teaching signal is used.}
\item{exemplar.locations}{}
\item{random.locations}{}
\item{n.locations}{}
\item{fixed}{}
\item{parStruct}{}
\item{data}{}
\item{subset}{}
\item{ntimes}{}
\item{replicate}{}
\item{base}{}
}

\details{ALCOVE (Kruschke, 1992) is based on the \code{\link{gcm}} model, but
has a mechanism to learn the attention weights. It is formulated as an ANN.
}
\value{A (fitted) object of class \code{ALCOVE}}
\references{Krusche (1992).}
\examples{
## open weather prediction data
data(WP)
## initialize model
mod <- alcove(y~x1+x2+x3+x4-1,data=WP,ntimes=c(200,200))
## estimate free parameters
mod <- estimate(mod)
summary(mod)
}
\keyword{models}