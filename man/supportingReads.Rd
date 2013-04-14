\name{supportingReads}
\alias{supportingReads}
\title{A function to extract supporting reads values from a list of fSet object}
\description{A function extracting supporting reads values from a list of fSet objects.}
\usage{supportingReads(list, parallel=F)}
\arguments{
  \item{list}{a list of fSet objects}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
supporting.reads <- supportingReads(tmp)
supporting.reads
}
\keyword{utilities}

