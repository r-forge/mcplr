\name{McplModel}
\alias{McplModel}
\alias{getPars,McplModel,ANY-method}

\title{McplModel}
\description{
Constructs an McplModel.
}
\usage{
McplModel(learningModel,responseModel)
}
\arguments{
\item{learningModel}{an object of class \code{LearningModel}}
\item{responseModel}{an object of class \code{ResponseModel}}
}

\details{}
\value{An object of class \code{McplModel}.}

\examples{
## open weather prediction data
data(WPT)
## first specify a learning model, then a response model, and then combine
##   them into an McplModel.
lMod <- RescorlaWagner(y~x1+x2+x3+x4,data=WPT,ntimes=rep(200,25),remove.intercept=TRUE)
rMod <- RatioRuleResponse(r~predict(lMod)-1,data=WPT,ntimes=rep(200,25))
mMod <- McplModel(lMod,rMod)
## estimate free parameters
\dontrun{mMod <- fit(mMod)}
summary(mMod)
}
