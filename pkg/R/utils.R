fixedListToVec <- function(fixed,parameters) {
  fix <- relist(rep(FALSE,length(unlist(parameters))),skeleton=parameters)
  fix.names <- names(fixed)
  for(i in fix.names) {
    if(is.logical(fixed[[i]])) {
      if(length(fix[[i]]) > 0) {
        if(is.matrix(fix[[i]])) {
          fix[[i]] <- matrix(rep(fixed[[i]],length=length(fix[[i]])),ncol=ncol(fix[[i]]),nrow=nrow(fix[[i]]))
        } else
          fix[[i]] <- rep(fixed[[i]],length=length(fix[[i]]))
      } else {
        warning("element ",i," in fixed not in parameter list")
      }
    } else {
      stop("Recursive lists not implemented (yet) in reconstruction of fixed vector from fixed list. Please provide fixed in vector form")
    }
  }
  fix
}

