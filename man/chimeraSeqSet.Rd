\name{chimeraSeqSet}
\alias{chimeraSeqSet}
\title{A function to generate s DNAStringSet encompassing fusion sequences}
\description{A function generating the nucleotide sequences of chimeras described in a list of fSet, i.e. 
the list generated using importFusionData function.             
}
\usage{chimeraSeqSet(list, parallel=FALSE)}
\arguments{
  \item{list}{A list of fSet objects.}
  \item{parallel}{If TRUE uses the BioParallel package}
}
\value{
    A DNAStringSet encompassing the fusions described in a list of fSet objects. This object represents the ideal reference to
	remap reads over detected fusions. Remapping is required to validate fusions using GapFiller de novo reconstruction. 
}
\author{Raffaele A Calogero}
\examples{
	tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
	fusion.names <- fusionName(tmp)
	fusion.names
	myset <- tmp[1:3]
	tmp.seq <- chimeraSeqSet(myset, parallel=FALSE)
	#writeXStringSet(tmp.seq, "detected.fusions.fa", format="fasta")
}

\seealso{\code{\link{fusionName}, \link{importFusionData}, \link{gapfillerInstallation}, \link{gapfillerRun}}}
\keyword{utilities}


