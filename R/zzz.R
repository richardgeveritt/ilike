# Package-level environment to cache include/lib flags captured at load time,
# when all LinkingTo dependencies are guaranteed to be on .libPaths().
.ilike_env <- new.env(parent = emptyenv())
.ilike_env$sourcecpp_cppflags <- ""
.ilike_env$sourcecpp_libs <- ""

.onLoad <- function(libname, pkgname) {
  depends_pkgs = c(pkgname, "Rcpp", "RcppArmadillo", "BH", "dqrng", "sitmo")
  include_flags = character(0)
  lib_flags = character(0)
  for (pkg in depends_pkgs)
  {
    inc = system.file("include", package = pkg)
    if (nchar(inc) > 0)
      include_flags = c(include_flags, paste0("-I", inc))
    if (pkg == "RcppArmadillo")
    {
      ld = tryCatch(RcppArmadillo::RcppArmadilloLdFlags(), error = function(e) "")
      if (nchar(ld) > 0)
        lib_flags = c(lib_flags, ld)
    }
  }
  .ilike_env$sourcecpp_cppflags = paste(include_flags, collapse = " ")
  .ilike_env$sourcecpp_libs = paste(lib_flags, collapse = " ")
}

.onUnload <- function (libpath) {
  library.dynam.unload("ilike", libpath)
}
