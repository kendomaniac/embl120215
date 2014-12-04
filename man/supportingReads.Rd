\name{supportingReads}
\alias{supportingReads}
\title{A function to extract supporting reads values from a list of fSet object}
\description{A function extracting supporting reads values from a list of fSet objects. Please note that not all outputs of supported tools provides both spanning reads, i.e. pair-end reads having one of the two mates spanning over the break point, and encompassing reads, i.e. pair-end reads having the two mates mapping on the exons of the two transcripts involved in the fusion. The presence of both type of reads is mandatory to provide the full number of reads covering the junction region. To know which information are provided by the supported tool please check \code{\link{fSet}}}
\usage{supportingReads(list, fusion.reads=c("encompassing","spanning"), parallel=FALSE)}
\arguments{
  \item{list}{a list of fSet objects}
  \item{fusion.reads}{it allows to extract spanning reads associated to the SeedCount slot or encompassing reads associated to the RescuedCount}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hg19")
supporting.reads <- supportingReads(tmp, fusion.reads="spanning")
supporting.reads
}
\keyword{utilities}


