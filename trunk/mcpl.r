mcpl <- function(learningModel,responseModel,fit=TRUE) {
  mod <- new("mcplModel",
    learningModel=learrningModel,
    responseModel=responseModel)
  if(fit) mod <- fit(mod)
  mod
}