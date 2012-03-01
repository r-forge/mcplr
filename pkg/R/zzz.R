.onLoad <- function (libname, pkgname) {
    library.dynam(pkgname, pkgname, lib.loc=libname)
}

.onUnload <- function (libpath) {
    library.dynam.unload("mcplR", libpath)
}

