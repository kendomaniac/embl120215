\name{validateSamFile}
\alias{validateSamFile}
\title{A function to validate SAM or BAM files}
\description{A function to validate SAM or BAM files using picard-tools           
}
\usage{validateSamFile(input, output, mode=c("VERBOSE", "SUMMARY"), max.output="100")}
\arguments{
  \item{input}{SAM/BAM file to be validated}
  \item{output}{file name in which to save the validation information}
  \item{mode}{Mode of output. Default value: VERBOSE}
  \item{max.output}{max number of reported error lines}
}
\value{
    Validation information referring to a SAM/BM file. 
}
\author{Raffaele A Calogero}
\examples{
#validateSamFile(input=paste(find.package(package="chimera"),"/examples/mcf7_trs_accepted_hits.bam",sep=""), output="validation_results.txt", mode="SUMMARY", max.output="100")
}

\seealso{ \code{\link{picardInstallation}}}
\keyword{utilities}


