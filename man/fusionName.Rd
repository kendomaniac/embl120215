\name{fusionName}
\alias{fusionName}
\title{A function to extract fusion names for a list of fSet object}
\description{A function allowing extract fusion names from a list of fSet objects.}
\usage{fusionName(list, parallel=FALSE)}
\arguments{
  \item{list}{a list of fSet object}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
fusion.names <- fusionName(tmp)
fusion.names
}
\keyword{utilities}


