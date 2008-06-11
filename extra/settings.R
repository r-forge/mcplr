# These functions are modelled after those in Lattice

mcpl.defaultOptions <- function() {
  list(
    mcpl.options=list(
      estimate = list(
        default.method="Nelder-Mead"
      ),
      logLik = list(
        na.action=c("keep","remove")
      )
    )
  )
}
mcpl.getOption <- function(name) {
    get("mcpl.options", envir = .McplEnv)[[name]]
}
mcpl.options <- function(...) {
    updateList <- function(x, val) {
        if (is.null(x)) x <- list()
        modifyList(x, val)
    }
    new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
    old <- .McplEnv$mcpl.options
    if (length(new) == 0) return(old)
    nm <- names(new)
    if (is.null(nm)) return(old[unlist(new)])
    isNamed <- nm != ""
    if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])
    retVal <- old[nm]
    names(retVal) <- nm
    nm <- nm[isNamed]
    .McplEnv$mcpl.options <- updateList(old, new[nm])
    invisible(retVal)
}