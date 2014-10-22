AdaptiveMetropolis
==================

[![Build Status](https://travis-ci.org/tom-jin/AdaptiveMetropolis.png?branch=master)](https://travis-ci.org/tom-jin/AdaptiveMetropolis)

An R package implementing a series of adaptive metropolis hastings Markov chain algorithms.

Written originally as a project to satisfy the assesment requirements of the Statistical Computing module from the [Oxford-Warwick Statistics Programme](http://www2.warwick.ac.uk/fac/sci/statistics/oxwasp/) [EPSRC Center for Doctoral Traning](http://www.epsrc.ac.uk/skills/students/centres/).

## Install

The devtools package is required to install development versions of this package with the following command:

```R
devtools::install_github("tom-jin/AdaptiveMetropolis")
```

## Demo

Start with a Normal(0,1) distribution explore a Unif(0,100) distribution for 10 steps. Then adapt the variance of the proposal distribution for enhanced mixing whilst exploring the target distribution.

```R
data <- AM(function(x) {dnorm(x, 0,100)}, 5000, 0, 1)
plot(data[1001:5000])
```
