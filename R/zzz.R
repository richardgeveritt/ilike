# Package-level environment.
.ilike_env <- new.env(parent = emptyenv())
.ilike_env$libname <- ""

.onLoad <- function(libname, pkgname) {
  # Store libname so compile() can include it in the filesystem search for
  # package include directories (belt-and-suspenders fallback).
  .ilike_env$libname <- libname

  # On Linux/Windows, R loads package shared libraries with RTLD_LOCAL by
  # default, meaning ilike's symbols (including the vtable for `Parameters`)
  # are NOT in the global symbol table.  When user code compiled by
  # sourceCpp() creates a `Parameters` object its lazy binding resolves to
  # a *different* vtable copy, causing a segfault when the object crosses the
  # DLL boundary into ilike.so.  Re-loading ilike.so with RTLD_GLOBAL puts
  # every ilike symbol into the global table so all DLLs share one vtable.
  # On macOS the dynamic linker uses a flat namespace by default, so this
  # is a harmless no-op there.
  if (.Platform$OS.type != "windows") {
    ilike_so = tryCatch({
      so_path = system.file("libs", paste0("ilike.so"), package = "ilike")
      if (!nzchar(so_path)) {
        # Some platforms put the .so in a sub-arch directory
        r_arch = .Platform$r_arch
        if (nzchar(r_arch))
          so_path = system.file("libs", r_arch, paste0("ilike.so"), package = "ilike")
      }
      so_path
    }, error = function(e) "")
    if (nzchar(ilike_so) && file.exists(ilike_so)) {
      tryCatch(
        dyn.load(ilike_so, local = FALSE, now = FALSE),
        error = function(e) invisible(NULL)   # already loaded is fine
      )
    }
  }
}

.onUnload <- function (libpath) {
  library.dynam.unload("ilike", libpath)
}
