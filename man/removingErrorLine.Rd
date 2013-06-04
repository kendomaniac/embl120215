\name{removingErrorLine}
\alias{removingErrorLine}
\title{A function to remove a line stopping SAM to BAM conversion}
\description{A function to remove a line stopping SAM to BAM conversion}
\usage{removingErrorLine(line.number, filenameIn, filenameOut)}
\arguments{
  \item{line.number}{line number to be removed}
  \item{filenameIn}{input file name}
  \item{filenameOut}{output file name}
}
\value{
       SAM file without the error line
}
\author{Raffaele A Calogero}
\examples{
#removingErrorLine(14680066,"kd2_accepted_hits.sam","kd2_accepted_hits1.sam")
}
\keyword{utilities}


