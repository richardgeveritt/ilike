# ilike
Implementations of a number of ABC algorithms and pseudo-marginal methods for dynamic models. Documentation at https://richardgeveritt.github.io/ilike_documentation/

<!-- badges: start -->
[![R-CMD-check](https://github.com/richardgeveritt/ilike/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/richardgeveritt/ilike/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

### 1. Install the HDF5 system library

`ilike` requires the HDF5 C library. Install it for your platform before installing the R package.

**macOS**

```bash
brew install hdf5
```

**Linux (Debian / Ubuntu)**

```bash
sudo apt-get install libhdf5-dev
```

**Linux (Fedora / RHEL / CentOS)**

```bash
sudo dnf install hdf5-devel
```

**Windows**

HDF5 is distributed via the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) toolchain. Open an Rtools MSYS2 shell (e.g. `C:\rtools45\ucrt64.exe`) and run:

```bash
pacman -Sy mingw-w64-x86_64-hdf5
```

> **Note:** Rtools must already be installed. `pacman` is the package manager bundled with Rtools 4.x.

### 2. Install ilike

```r
# install.packages("remotes")
remotes::install_github("richardgeveritt/ilike")
```
