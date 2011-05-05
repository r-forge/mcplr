.onLoad <- function (libname, pkgname) {
    library.dynam(pkgname, pkgname, lib.loc=libname)
}

.onUnload <- function (libname, pkgname) {
    library.dynam.unload(pkgname)
}

