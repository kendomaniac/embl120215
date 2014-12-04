\name{is.fSet}
\alias{is.fSet}
\title{A function to evaulate if an object belongs to fSet class or not}
\description{A function to evaulate if an object belongs to fSet class or not.}
\usage{is.fSet(x)}
\arguments{
  \item{x}{an object}
}
\value{
       If the object belons to fSet class it returns TRUE, else it returned FALSE
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hg19")
is.fSet(tmp[[1]])
}
\keyword{utilities}


