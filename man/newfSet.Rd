\name{newfSet}
\alias{newfSet}
\title{A constructor for fSet class objects}
\description{A function to create a new fSet object}
\usage{newfSet(fusionInfo = list(fusionTool = NA,
                                      UniqueCuttingPositionCount = 0,
                                      SeedCount = 0, 
                                      RescuedCount = 0, 
                                      SplicePattern = NA, 
                                      FusionGene = NA,
                                      frameShift = NA), 
	           fusionLoc = GRangesList(),
               fusionRNA = DNAStringSet(),
               fusionGA = GAlignments())
}

\arguments{
         \item{fusionInfo}{A list of the fusion characteristics see fSet class slot fusionInfo}
         \item{fusionLoc}{A GRangesList encompassing the fusion break points}
         \item{fusionRNA}{A DNAStringSet encompassing the fusion transcript}
         \item{fusionGA}{A GAlignments encompassing the location of reads mapping on the fusion transcript}
}

\value{
       An object of fSet class
}
\author{Raffaele A Calogero}
\examples{
tmp <- newfSet()
tmp
}
\keyword{utilities}


