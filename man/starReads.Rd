\name{starReads}
\alias{starReads}
\title{A function to extract reads info from STAR fusion output}
\description{A function producing a GRangeList for the reads information, involved in a fusion event.}
\usage{starReads(fusion.report, parallel=FALSE)}
\arguments{
  \item{fusion.report}{STAR fusion output file}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
#tmp <- starReads("Chimeric.out.junction", parallel=FALSE)

}
\keyword{utilities}

