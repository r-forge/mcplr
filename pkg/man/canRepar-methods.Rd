\name{canRepar-methods}
\docType{methods}
\alias{canRepar-methods}
\alias{canRepar}

\title{Can a model be reparametrized?}
\description{
The \code{canRepar} method is used to check whether a model can be 
reparametrized. For certain models, reparametrization can be useful to speed up 
convergence. By default, the method returns FALSE. By specifying the method 
for particular models, it can be made to return TRUE. In that case, if \code{setTransPars} and 
\code{getTransPars} methods exists for the model, they will be used when 
\code{fit} method is called on the model.
}
\section{Methods}{
\describe{

\item{object = "ANY"}{Generic function.}

}}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
