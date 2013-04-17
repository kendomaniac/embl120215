\name{supportingReads}
\alias{supportingReads}
\title{A function to extract supporting reads values from a list of fSet object}
\description{A function extracting supporting reads values from a list of fSet objects.}
\usage{supportingReads(list, fusion.reads=c("all","spanning"), parallel=F)}
\arguments{
  \item{list}{a list of fSet objects}
  \item{fusion.reads}{it allows to extract spanning reads associated to the SeedCount slot or all the detected reads associate to the RescuedCount}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
supporting.reads <- supportingReads(tmp, fusion.reads="spanning")
supporting.reads
}
\keyword{utilities}

