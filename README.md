# SWTest

An R package which implements the methods in [Testing for the Network Small-World Property](https://arxiv.org/pdf/2103.08035.pdf). The package provides the function `sw_test()` which takes as input a network and a baseline model and performs a hypothesis test as described in the paper.

## Installation

``` r
library(devtools)
install_github(KartikSL/SWTest)
```

## Example Usage

``` r
library(SWTest)
data(dolphins)
sw_test(dolphins, c("ER", "DCSBM"), c(500, 50))
```
