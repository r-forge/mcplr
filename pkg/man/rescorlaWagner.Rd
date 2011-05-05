\name{RescorlaWagner}
\alias{RescorlaWagner}
\alias{RescorlaWagner-class}
\alias{ContinuousRescorlaWagner-class}
\alias{fit,RescorlaWagner-method}
\alias{runm,RescorlaWagner-method}
\alias{runm,ContinuousRescorlaWagner-method}
\alias{getPars,RescorlaWagner-method}
\alias{predict,RescorlaWagner-method}
\alias{setPars,RescorlaWagner-method}

\title{Rescorla-Wagner Model}
\description{
Constructs a Rescorla-Wagner learning model.
}
\usage{
RescorlaWagner(formula,parameters=list(alpha=.1,beta=1,lambda=1,ws=0),data,subset,
fixed=list(alpha=FALSE,beta=TRUE,lambda=TRUE,ws=TRUE),parStruct,remove.intercept=FALSE,
base=NULL,ntimes=NULL,replicate=TRUE)
}
\arguments{
\item{formula}{an object of class \code{formula} (or one that can be coerced to 
 that class): a symbolic description of the model to be fitted. For more details
 of model speecification, see \code{lm} or \code{glm}.}
\item{parameters}{an optional named list with (starting) values of the 
parameters. Missing parameters in the list will be set to default values.}
\item{data}{(optional) data frame for evaluation of the formula.}
\item{subset}{(optional) subset of the data.}
\item{fixed}{(optional) logical vector indicating whether parameters are fixed 
(TRUE) or freely estimable (FALSE).}
\item{parStruct}{(optional) ParStruct object. Note that if parStruct is given, 
the `fixed' argument above will be ignored.}
\item{remove.intercept}{should an intercept be included in the predictors? If 
set to TRUE, an intercept term (if present) will be removed from x. This can 
be useful when the predictors are categorical variables, where removing the 
intercept in the formula will result in an extra dummy variable.}
\item{base}{which level of the criterion variable is considered the base 
category? If set to a positive integer, that level will not be included in the
criterion matrix y.}
\item{ntimes}{an optional vector with, for each repetition in the data, the 
total number of trials.}
\item{replicate}{are the repeated series true replications, i.e., are the model
parameters considered identical for each series?}
}
\details{The Rescorla-Wagner model (Rescorla & Wagner, 1972) is an associative
learning model.

The associative strength of the connection between an CS j and US k is updated as:

\deqn{w_{jk}(t) = w_{jk}(t-1) + [\alpha_j \beta_k](\lambda_k y_k - \sum_i w_{ik}(t-1) x_i(t))}{w[j,k,t] = w[j,k,t-1] + (alpha[j]*beta[k,m])*(lambda[k] y[t,k] - sum[i] w[i,k,t-1] x[t,i])*x[t,j]}

In this expression \eqn{\lambda_k}{lambda[k]} is a scaling parameter; when \eqn{y_k}{y[k]} takes values 0 and 1 when US k is absent/present respectively, then \eqn{\lambda_k}{lambda[k]} reflects the maximum conditioning US k can produce, i.e.,  it represents the limit of learning. The \eqn{\alpha_j}{alpha[j]} and \eqn{\beta_k}{beta[k]} are rate parameters dependent, respectively, on the CS \eqn{x_j}{x[j]} and US \eqn{y_k}{y[k]}. These parameters are viewed as having fixed values based on the physical properties of the particular CS and US. On any given trial the current associative strength, \eqn{\sum_i w_{ik}x_i(t)}{sum[i] w[ik]x[i](t)}, is compared with \eqn{\lambda_k y_k}{lambda[k]*y[k]} and the difference is treated like an error to be corrected; this happens by producing a change in associative strength accordingly. Consequently, the Rescorla-Wagner model is an error-correction model.

}
\value{An object of class \code{LearningModel}}
\references{
Rescorla, R. A., and Wagner, A. R. (1972). A theory of Pavlovian conditioning: 
Variations in the effectiveness of reinforcement and nonreinforcement. In A. H. 
Black and W. F. Prokasy (Eds.), \emph{Classical Conditioning II}, 
Appleton-Century-Crofts. pp. 64-99.
}
\author{Maarten Speekenbrink}
\examples{
## Open weather prediction data
data(WPT)
## Initialize model
lMod <- RescorlaWagner(y~x1+x2+x3+x4,data=WPT,ntimes=c(200,200),
remove.intercept=TRUE)
## Compute associative strengths
lMod <- runm(lMod)
## Estimate free parameters
\dontrun{lMod <-fit(lMod)}
summary(lMod)
\dontrun{
## Now include a response model
rMod <- RatioRuleResponse(r~1,data=WPT)
mMod <- McplModel(lMod,rMod)
fMod <- fit(mMod)
summary(fMod)
}
## The following model only learns to expect the presence/absence of the second 
## outcome, which completely determines the presence/absence of the first 
## outcome.
lMod2 <- RescorlaWagner(as.numeric(y==1)~x1+x2+x3+x4,data=WPT,ntimes=c(200,200),
remove.intercept=TRUE)
# compute associative strengths
lMod2 <- runm(lMod2)
\dontrun{
## Now include a response model. Note that the response in the model should 
## have the same format as the outcome in the learning model. 
rMod2 <- GlmResponse(as.numeric(r==1)~predict(lMod2)-1,family=binomial(),data=WPT)
mMod2 <- McplModel(lMod2,rMod2)
fMod2 <- fit(mMod2)
summary(fMod2)
}


}
