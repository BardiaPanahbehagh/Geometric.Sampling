\name{GFS Sampling}
\alias{GFSfixedsize}
\title{Fixed size GFS sampling}
\usage{
GFSfixedsize(pik,v,chaos,M.chaos,secondor,plot,counter)
}
\description{
Using this function it is possible to implement unequal probability sampling based on GFS method.
}
\arguments{
 \item{pik}{vector of inclusion probabilities}
 \item{v}{size of the segments}
 \item{chaos}{increasing entropy, leads to decreasing zero second orders}
 \item{M.chaos}{number of times that we exchange the segmentss for increasing entropy}
 \item{secondor}{1:second order of the sample, 2:second order of the population}
 \item{plot}{1:plot}
 \item{counter}{1:active}
}
\examples{
GFSfixedsize(
  pik=rep(.2,100)
  ,v=.01
  ,chaos=1
  ,M.chaos=100000
  ,secondor=0
  ,plot=1
)
}

