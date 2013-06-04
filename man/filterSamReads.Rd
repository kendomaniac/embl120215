\name{filterSamReads}
\alias{filterSamReads}
\title{A function to filter SAM or BAM files}
\description{A function to filter SAM or BAM files using picard-tools           
}
\usage{filterSamReads(input, output, filter=c("includeAligned","excludeAligned"))}
\arguments{
  \item{input}{SAM/BAM file to be validated}
  \item{output}{file name in which to save the filtered reselts}
  \item{filter}{type of filter}
}
\value{
    A filtered SAM/BAM. 
}
\author{Raffaele A Calogero}
\examples{
#	filterSamReads(input="kd2_accepted_hits2.sam", output="kd2_accepted_hits2_mapped.sam", filter="includeAligned")
}

\seealso{ \code{\link{picardInstallation}}}
\keyword{utilities}


