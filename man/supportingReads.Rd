\name{supportingReads}
\alias{supportingReads}
\title{A function to extract supporting reads values from a list of fSet object}
\description{A function extracting supporting reads values from a list of fSet objects.}
\usage{supportingReads(list)}
\arguments{
  \item{list}{a list of fSet objects}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
supporting.reads <- supportingReads(tmp)
supporting.reads
}
\keyword{utilities}

