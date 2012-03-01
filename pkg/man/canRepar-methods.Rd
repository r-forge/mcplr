\name{canRepar-methods}
\docType{methods}
\alias{canRepar}
\alias{canRepar-methods}
\alias{canRepar,ANY-method}
\alias{getTransPars}
\alias{getTransPars-methods}
\alias{setTransPars}
\alias{setTransPars-methods}
\alias{getTransPars,ANY-method}
\alias{canRepar,McplBaseModel-method}
\alias{canRepar,McplModel-method}
\alias{canRepar,LearningModel-method}
\alias{canRepar,ResponseModel-method}
\alias{setTransPars,ANY,ANY-method}

\alias{canRepar,GCMlearning-method}
\alias{canRepar,RescorlaWagner-method}
\alias{getTransPars,GCMlearning-method}
\alias{getTransPars,RescorlaWagner-method}
\alias{setTransPars,GCMlearning-method}
\alias{setTransPars,McplModel-method}
\alias{setTransPars,RescorlaWagner-method}



\title{Reparametrization methods}
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
