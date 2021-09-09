# Introduction

This is a completely revamped version of a package implementing two algorithms for two dimensional phase unwrapping that can be set up and solved as [linear programs](https://en.wikipedia.org/wiki/Linear_programming). The revisions were motivated both for reasons of speed and maintainability. The LP solver I previously used appears no longer to be in active development, while the one ([CBC](https://github.com/coin-or/Cbc))I chose for this version is well funded and supported by the [COIN-OR](https://www.coin-or.org/) project. In limited testing I found the new solver can be orders of magnitude faster than the previous one. The speed difference is especially noticeable in the problem solved by `netflowpuw`, which could previously take some minutes or even hours to execute.

# Installation

Two non-CRAN packages are required and should be installed before this one. First is my package [zernike](https://github.com/mlpeck/zernike), which can be installed from source or from the Windows binary provided in the Releases section of my Github project. Also required is [rcbc](https://dirkschumacher.github.io/rcbc/index.html), which provides an R interface to the CBC LP solver. Please contact me at mlpeck54 -at- gmail.com if you need a Windows binary as none is provided at present.

Both of the functions in this package use sparse matrices defined in the package Matrix, which is part of the standard R installation. Finally you need either `remotes` or `devtools` to install this package. Both are available on CRAN.

This package is written entirely in R and can easily be installed from source without any additional tools. From the R command prompt in a console window or terminal just type `remotes::install_github("mlpeck/lppuw")`. This will download and install the source package. The installed binary can be loaded into your workspace with `load(lppuw)`.

# Basic usage

Usage is simple: given a modulo 2&pi; phase map `phi` and optionally a modulation map `mod` in matrices of the same dimension the two available functions are called with:

```
wf.bc <- brcutpuw(phi)
wf.nf <- netflowpuw(phi, wts=mod)
```

The `wts` argument in `netflowpuw` is optional but in my experience using the modulation map in the cost function gives superior results to the default.

There are sample phase and modulation maps in the data directory (data courtesy Vladimir Galogaza) which can be loaded with `data("phasemaps", package="lppuw")`. There is also a usage example that can be run with `example(brcutpuw, package="lppuw", ask=FALSE)`.

