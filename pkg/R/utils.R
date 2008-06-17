fixedListToVec <- function(fixed,parameters) {
  tfix <- relist(rep(FALSE,length(unlist(parameters))),skeleton=parameters)
  fix.names <- names(fixed)
  for(i in fix.names) {
    if(is.logical(fixed[[i]])) {
      if(length(tfix[[i]]) > 0) {
        if(is.matrix(tfix[[i]])) {
          tfix[[i]] <- matrix(rep(fixed[[i]],length=length(tfix[[i]])),ncol=ncol(tfix[[i]]),nrow=nrow(tfix[[i]]))
        } else
          tfix[[i]] <- rep(fixed[[i]],length=length(tfix[[i]]))
      } else {
        warning("element",i,"in fixed not in parameter list")
      }
    } else {
      stop("Recursive lists not implemented (yet) in reconstruction of fixed vector from fixed list. Please provide fixed in vector form")
    }
  }
}