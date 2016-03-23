\name{csm-package}
\alias{csm-package}
\docType{package}
\title{
Causal State Models - inferred by CSSR and structually equivalent to Hidden Markov Models.
}
\description{
This package provides tools for creation of Causal State Models. They are structually equivalent to Hidden Markov Models and inferred from data by CSSR algorithm using a distribution equality statistical test (at the moment only support for chi-square test is implemented). There is "csm" class and methods for it and CSSR algorithm implementation for inference of CSM from data.
}
\details{
\tabular{ll}{
Package: \tab csm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-07\cr
License: \tab GPL (>= 2)\cr
}
Use csm(states, symbols, transArr) to create a predefined Causal State Model, use simulate(csm, nsim) to create a simulation of a sequence of the model, use CSSR(data, L, alpha) to infer Causal State Model from data.
}
\author{
Rudolfs Petrovs

Maintainer: Rudolfs Petrovs <rudolfs.petrovs@gmail.com>
}
\references{
Shalizi, Cosma Rohilla, Kristina Lisa Shalizi, and James P. Crutchfield. "An algorithm for pattern discovery in time series." arXiv preprint cs/0210025 (2002).

Shalizi, Cosma Rohilla, and Kristina Lisa Shalizi. "Blind construction of optimal nonlinear recursive predictors for discrete sequences." In Proceedings of the 20th conference on Uncertainty in artificial intelligence, pp. 504-511. AUAI Press, 2004.
}

\keyword{ package }
\seealso{
Main algorithm: \code{\link{CSSR}}

"csm" class is provided by this package: \code{\link{csm}}

A useful and important method implemented: \code{\link{simulate.csm}}

}
\examples{
even <- csm(data.frame(c("A", "A", "B"),c("A", "B", "A"), c(0, 1, 1), c(0.5, 0.5, 1)))
data <- simulate(even, 1000)$symbols
acsm <- CSSR(data, 3, 0.001)

}
