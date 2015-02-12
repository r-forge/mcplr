genListToVec <- function(given,skeleton,def.value,def.class,name) {
  gv <- relist(rep(def.value,length(unlist(skeleton))),skeleton=skeleton)
  gv.names <- names(given)
  for(i in gv.names) {
    if(is(given[[i]],def.class)) {
      if(length(gv[[i]]) > 0) {
        if(is.matrix(gv[[i]])) {
          gv[[i]] <- matrix(rep(given[[i]],length=length(gv[[i]])),ncol=ncol(gv[[i]]),nrow=nrow(gv[[i]]))
        } else
          gv[[i]] <- rep(given[[i]],length=length(gv[[i]]))
      } else {
        warning("element ",i," in ",name," not in parameter list")
      }
    } else {
      stop("Recursive lists not implemented (yet) in reconstruction of ",name," vector from ",name," list. Please provide ",name," in vector form")
    }
  }
  gv
}

genRepListToVec <- function(given,skeleton,cases,def.value,def.class,name) {
  gv <- vector("list",length=cases)
  if(is.list(given)) {
    if(is.null(names(given)) && length(given) == cases) {
      for(i in 1:cases) {
        if(is.list(given[[i]])) {
          gv[[i]] <- genListToVec(given[[i]],skeleton,def.value,def.class,name)
        } else {
          gv[[i]] <- given[[i]]
        }
      }
    } else {
      for(i in 1:cases) {
        gv[[i]] <- genListToVec(given,skeleton[[i]],def.value,def.class,name)
      }
    }
  } #else {
    #if(length(given) == length(unlist(skeleton))) {
    #  gv <- given
    #} else {
    #  gv <- rep(given,length=length(unlist(skeleton)))
    #}
  #}
  gv
}

fixedListToVec <- function(fixed,parameters) {
  out <- genListToVec(fixed,parameters,FALSE,"logical","fixed")
  out
}

fixedRepListToVec <- function(fixed,parameters,cases) {
  out <- genRepListToVec(fixed,parameters,cases,FALSE,"logical","fixed")
  out
}

idListToVec <- function(id,parameters) {
  out <- genListToVec(id,parameters,-999,"numeric","id")
  rem.ids <- seq(from=1,to=(sum(out == -999) + length(unique(out))))
  rem.ids <- rem.ids[!(rem.ids %in% unique(out))]
  rem.ids[1:sum(out == -999)]
  out[out == -999] <- rem.ids
  out
}

idRepListToVec <- function(id,parameters,cases) {
  out <- genRepListToVec(id,parameters,cases,-999,"numeric","id")
  rem.ids <- seq(from=1,to=(sum(out == -999) + length(unique(out))))
  rem.ids <- rem.ids[!(rem.ids %in% unique(out))]
  rem.ids[1:sum(out == -999)]
  out[out == -999] <- rem.ids
  out
}

minListToVec <- function(minimum,parameters) {
  out <- genListToVec(minimum,parameters,-Inf,"numeric","minimum")
  out
}

minRepListToVec <- function(minimum,parameters,cases) {
  out <- genRepListToVec(minimum,parameters,cases,-Inf,"numeric","minimum")
  out
}
  
maxListToVec <- function(maximum,parameters) {
  out <- genListToVec(maximum,parameters,Inf,"numeric","maximum")
  out
}

maxRepListToVec <- function(maximum,parameters,cases) {
  out <- genRepListToVec(maximum,parameters,cases,Inf,"numeric","maximum")
  out
}