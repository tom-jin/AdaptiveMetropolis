% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ASWAM}
\alias{ASWAM}
\title{Adaptive Scaling Within Adaptive Metropolis Sampler}
\usage{
ASWAM(target, n, init.state, init.cov, burn.in = 1000)
}
\arguments{
\item{target}{The target distribution where the unscaled densities can be
evaluated.}

\item{n}{The number of samples to produce.}

\item{init.state}{The initial state of the sampler.}

\item{init.cov}{The covariance of the sampler during burn in.}

\item{burn.in}{(Optional) Number of iterations to burn in for. Default: 1000}
}
\value{
\code{n} samples from the target distribution.
}
\description{
An implementation of the Adaptive Scaling within Adaptive Metropolis algorithm
as published by Andrieu and Thoms (2008).
The adapted covariance matrix which encodes knowledge of the second moments
of the target distribution is scaled to achieve a target acceptance ratio.
}
\details{
The proposal distribution is a multivariate normal distribution of the
appropriate dimension. The covariance of this distribution is computed from
previous samples and then scaled to aim for an acceptance rate of 0.234. This
is done to stike a balence between a high acceptance rate with poor mixing
and a low acceptance rate with periods of stationarity. Using diffusion
limits an acceptance rate of 0.234 has been found to be optimal.

The dimension is detected using the length of the \code{init.state}
parameter.

Adaptive samplers do not attempt to solve the problem of multi modal
distributions. If the target distribution is known to have multiple modes,
particually if they are far apart, more advanced techniques are recomended
such as nesting adaptive Metropolis within parallel tempering.
}
\examples{
# 4D Multivariate normal
library(mvtnorm)
data <- ASWAM(function(x) {dmvnorm(x, rep(10, 4), 10*diag(4))}, 8000, rep(0, 4), diag(4))
pairs(data)
}
\author{
Tom Jin <sjin@stats.ox.ac.uk>
}
\references{
Andrieu, C., Thoms, J. (2008) \emph{A tutorial on adaptive MCMC.}
Statistics and Computing, Vol. 18, Issue 4
}

