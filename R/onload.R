.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Conley package loaded...")
}

.onUnload <- function (libpath) {
  library.dynam.unload("conleyerrors", libpath)
}
