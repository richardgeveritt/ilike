# Package-level environment.
.ilike_env <- new.env(parent = emptyenv())
.ilike_env$libname <- ""

.onLoad <- function(libname, pkgname) {
  # Store libname so compile() can include it in the filesystem search for
  # package include directories (belt-and-suspenders fallback).
  .ilike_env$libname <- libname
}

.onUnload <- function (libpath) {
  library.dynam.unload("ilike", libpath)
}
