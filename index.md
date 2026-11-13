# ilike

Implementations of a number of ABC algorithms and pseudo-marginal
methods for dynamic models. Documentation at
<https://richardgeveritt.github.io/ilike_documentation/>

## Installation

`ilike` is installed from source, so it needs a working compiler
toolchain for R, including the Fortran libraries that R links against
(ilike uses LAPACK via RcppArmadillo), and the HDF5 C library.

### 1. Install the compiler toolchain

**macOS**

1.  Install the Xcode command line tools:

    ``` bash
    xcode-select --install
    ```

2.  Install the gfortran build that CRAN uses, from
    <https://mac.r-project.org/tools/>. It installs into
    `/opt/gfortran`, which is where R looks for the Fortran libraries.
    The gfortran **version must match the one your R was built with**:
    run

    ``` bash
    R CMD config FLIBS
    ```

    and look at the version number in the first path, e.g.
    `-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0` means
    gfortran 14.2 (used by R 4.5 and later). An older gfortran in
    `/opt/gfortran` gives the same linker errors as a missing one; the
    CRAN page explains how to remove it before installing the right
    version.

    Homebrew’s `gfortran` (`brew install gcc`) is **not** a substitute:
    it lives in a different place, so linking fails with errors such as

        ld: warning: search path '/opt/gfortran/lib' not found
        ld: library 'emutls_w' not found

    (You can see where R expects the libraries with
    `R CMD config FLIBS`.)

**Linux**

Install R’s development toolchain, which includes gfortran, e.g.
`sudo apt-get install r-base-dev` (Debian / Ubuntu) or
`sudo dnf install R-devel` (Fedora / RHEL).

**Windows**

Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
matching your R version; it includes the compilers and Fortran
libraries.

### 2. Install the HDF5 system library

`ilike` requires the HDF5 C library. Install it for your platform before
installing the R package.

**macOS**

``` bash
brew install hdf5
```

**Linux (Debian / Ubuntu)**

``` bash
sudo apt-get install libhdf5-dev
```

**Linux (Fedora / RHEL / CentOS)**

``` bash
sudo dnf install hdf5-devel
```

**Windows**

HDF5 is distributed via the
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) toolchain. Open
an Rtools MSYS2 shell (e.g. `C:\rtools45\ucrt64.exe`) and run:

``` bash
pacman -Sy mingw-w64-x86_64-hdf5
```

> **Note:** Rtools must already be installed. `pacman` is the package
> manager bundled with Rtools 4.x.

If you upgrade HDF5 to a new major version later (e.g. with
`brew upgrade hdf5`), reinstall ilike so that it links against the new
library.

### 3. Install ilike

``` r

# install.packages("remotes")
remotes::install_github("richardgeveritt/ilike")
```
