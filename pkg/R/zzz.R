.onLoad <- function(lib,pkg)
    library.dynam("mcplR", pkg,lib)

.onUnload <- function(libpath)
    library.dynam.unload("mcplR", libpath)
