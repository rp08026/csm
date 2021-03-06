% Generated by roxygen2 (4.0.2): do not edit by hand
\name{csm}
\alias{csm}
\title{Create a CSM object}
\usage{
csm(trans.df)
}
\arguments{
\item{trans.df}{a dataframe with transitions.}
}
\description{
Creates a Causal State Model object of class "csm" or "csmlike".
}
\details{
'trans.df' row specifies one transition. It is c("From","To","Symbol", "Probability") - first two are state names,
then emitted symbol and then this transition probability.
}
\examples{
coin <- csm(data.frame(rep("COIN", 2), rep("COIN", 2), c("H", "T"), rep(0.5, 2)))
}

