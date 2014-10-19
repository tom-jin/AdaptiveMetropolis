AdaptiveMetropolis
==================

An R package implementing a series of adaptive metropolis hastings Markov chain algorithms.

Written originally as a project to satisfy the assesment requirements of the Statistical Computing module from the Oxford-Warwick Statistics Programme (OxWaSP) EPSRC Center for Doctoral Traning.

## Install

The devtools package is required to install development versions of this package with the following command:

```R
devtools::install_github("tom-jin/AdaptiveMetropolis")
```

## Demo

Start with a Normal(0,1) distribution explore a Unif(0,100) distribution for 10 steps. Then adapt the variance of the proposal distribution for enhanced mixing whilst exploring the target distribution.

```R
data <- AM(function(x) {dunif(x, 0,100)}, 1000, 0, 1, 10)
plot(data)
```
