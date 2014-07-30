\name{breakpointOverlaps}
\alias{breakpointOverlaps}
\title{A function to extract the reads overlapping to fusion breakpoint.}
\description{A function to extract the reads overlapping to fusion breakpoint.}
\usage{breakpointOverlaps(fset, plot=FALSE, ylim=NULL)}

\arguments{
         \item{fset}{An fSet object. The slots fusionRNA and fusionGA needs to be loaded}
		 \item{plot}{If FALSE plot is not printed}
		 \item{ylim}{If NULL it uses the full fusion transcript coverage to define the ylim of the plot. If setted it can be used to zoom in the plot to better see the structure of the coverage at the brek point}
}

\value{
       An object of GAlignment class. A plot of the fusion trascript coverage in blue and of the reads spanning over the break point in yellow.
}
\author{Raffaele A Calogero}
\examples{
load(paste(find.package(package="chimera"), "/examples/fset_ARFGEF2-SULF2.rda", sep=""))
my.seq <- chimeraSeqs(my.fset)
my.fset <- addRNA(my.fset, my.seq)
tmp <- breakpointOverlaps(my.fset)
}
\keyword{utilities}


