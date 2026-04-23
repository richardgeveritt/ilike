# Package-level environment.
.ilike_env <- new.env(parent = emptyenv())
.ilike_env$libname <- ""

.onLoad <- function(libname, pkgname) {
  # Store libname so compile() can search for sibling packages there later.
  # On R CMD check, libname is the check temp library where LinkingTo packages
  # are installed alongside ilike.
  .ilike_env$libname <- libname
}

.onUnload <- function (libpath) {
  library.dynam.unload("ilike", libpath)
}
