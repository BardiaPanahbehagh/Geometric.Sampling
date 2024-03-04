\name{ optimal GFS Sampling}
\alias{GFSoptimal}
\title{Optimal GFS sampling}
\usage{
GFSoptimal(X,v,Lambda,chaos,M.chaos,secondor,plot,epsilon)
}
\description{
Using this function it is possible to implement unequal probability sampling based on GFS method.
}
\arguments{
 \item{X}{vector of auxiliary variable}
 \item{n}{expected sample size}
 \item{v}{size of segments}
 \item{Lambda}{number of tested design to find the optimal}
 \item{loss.f}{loss function, 1:variance, 2:minimax}
 \item{M.chaos}{number of times that we exchange the segmentss for increasing entropy and creating new designs}
 \item{secondor}{1:second order of the sample, 2:second order of the population}
 \item{plot}{1:plot}
 \item{counter}{1:active}
 \item{epsilon}{Shifting the auxiliary variable to be used as a new auxiliary variable}
}
\examples{
GFSoptimal(
  X          = runif(100, 10, 50)
  ,n         = 30
  ,v         = 0.01
  ,Lambda     = 5000
  ,loss.f    = 1
  ,M.chaos  = 100
  ,secondor  = 0
  ,plot      = 0
  ,epsilon     = 5
)
}

