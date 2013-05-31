\name{chimeraSeqs}
\alias{chimeraSeqs}
\title{A function to generate the nucleotide sequences of a fusion event}
\description{A function generating the nucleotide sequences of a chimera.             
}
\usage{chimeraSeqs(fset, extend=1000, type="transcripts")}
\arguments{
  \item{fset}{A fSet object.}
  \item{extend}{number of nucleotides used to extend a genomic region that is not an annotated gene. Default is 1000 nts}
  \item{type}{Chimera can be build at transcript level, i.e. transcript level will generate multiple fusion transcripts depending of the number of transcripts associated to each of the two genes involved in the fusion.}
}
\value{
    A DNAStringSet encompassing the fusions generated using all the isoforms for each gene involved in the fusion.The name of each element of the DNAStringSet has the following format: gene1-lengthOfGeneFragment:gene2-lengthOfGeneFragment.
In case the fusion junction is located in an intronic sequence, a warning is provided. The presence of a partial intron in the fusion is an indication that the fusion does not generate an active chimeric peptide.  
}
\author{Raffaele A Calogero}
\examples{
	tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
	fusion.names <- fusionName(tmp)
	fusion.names
	myset <- tmp[[13]]
	tmp.seq <- chimeraSeqs(myset, type="transcripts")
	#writeXStringSet(tmp.seq, paste(sub(":","_",fusion.names[[13]]),".fa",sep=""), format="fasta")
}

\seealso{\code{\link{fusionName}}}
\keyword{utilities}


