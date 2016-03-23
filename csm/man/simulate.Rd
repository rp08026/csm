\name{simulate}
\alias{simulate}
\alias{simulate.csm}
\alias{simulate.csmlike}
\title{Simulate a CSM and CSMlike objects}
\usage{
\method{simulate}{csm}(object, nsim = 1, seed = NULL, ...)
\method{simulate}{csmlike}(object, nsim = 1, seed = NULL, ...)
}
\arguments{
\item{object}{an object of class "csm" or "csmlike".}

\item{nsim}{a number of observations to simulate. Defaults to 1.}

\item{seed}{an object specifying if and how the random number generator should be initialized ('seeded').}

\item{\dots}{additional optional arguments.}
}
\description{
Simulate an observable sequence from the Causal State Model
}
\details{
The result of simulation is a sequence of hidden states under $states
and a sequence of observable emitted symbols under $symbols. Additional information
about used indexing is available under $statesInd, $symbolsInd and $transInd.
}
\examples{
even <- csm(data.frame(c("A", "A", "B"),c("A", "B", "A"), c(0, 1, 1), c(0.5, 0.5, 1)))
simulate(even, 1000)
}

