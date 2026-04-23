# Package-level environment to cache include/lib flags captured at load time,
# when all LinkingTo dependencies are guaranteed to be available.
.ilike_env <- new.env(parent = emptyenv())
.ilike_env$sourcecpp_cppflags <- ""
.ilike_env$sourcecpp_libs <- ""

.onLoad <- function(libname, pkgname) {
  find_pkg_include = function(pkg) {
    # Primary: system.file works when the package is on .libPaths().
    inc = system.file("include", package = pkg)
    if (nchar(inc) > 0 && dir.exists(inc)) return(inc)
    # Fallback: on R CMD check, LinkingTo-only packages are installed to the
    # same library as ilike (libname) but may not be on .libPaths() during
    # tests.  Search libname and all known library paths explicitly.
    for (lib in unique(c(libname, .libPaths()))) {
      candidate = file.path(lib, pkg, "include")
      if (dir.exists(candidate)) return(candidate)
    }
    return("")
  }

  depends_pkgs = c(pkgname, "Rcpp", "RcppArmadillo", "BH", "dqrng", "sitmo")
  include_flags = character(0)
  lib_flags = character(0)
  for (pkg in depends_pkgs)
  {
    inc = find_pkg_include(pkg)
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
