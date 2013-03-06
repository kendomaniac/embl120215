\name{plotCoverage}
\alias{plotCoverage}
\title{A function to plot the coverage of a fusion gene}
\description{A function to plot the coverage of a fusion gene.             
}
\usage{plotCoverage(fset, plot.type=c("exons","junctions"), junction.spanning=20, fusion.only=FALSE, xlab="nts", ylab="Coverage", main="", col.box1="red", col.box2="green", ybox.lim=c(-4,-1))
}
\arguments{
  \item{fset}{A fSet object}
  \item{plot.type}{exons plot exons coverage as junctions plot coverage of junction between exons}
  \item{junction.spanning}{number of nucleotides located upstream and downstream the junction/fusion location}
  \item{fusion.only}{if TRUE only fusion coverage is plotted}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{main}{Plot title}
  \item{col.box1}{color of the box describing the first gene}
  \item{col.box2}{color of the box describing the second gene}
  \item{ybox.lim}{y range defining the height of the box representing the exons}
}
\author{Raffaele A Calogero}
\examples{
	tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
	fusion.names <- fusionName(tmp)
	fusion.names
	myset <- tmp[[13]]
	trs <- chimeraSeqs(myset, type="transcripts")
	myset <- addRNA(myset , trs)
	tmp.seq <- fusionRNA(myset)
	myset <- addGA(myset, paste(path.package(package="chimera"),"/examples/mcf7_trs_accepted_hits.bam",sep=""))
	ga <- fusionGA(myset)
	plotCoverage(myset, plot.type="exons", col.box1="red", col.box2="green", ybox.lim=c(-4,-1))
	plotCoverage(myset, plot.type="junctions", col.box1="red", col.box2="yellow", ybox.lim=c(-4,-1))
	plotCoverage(myset, fusion.only=TRUE, col.box1="red", col.box2="yellow", ybox.lim=c(-4,-1))

}

\seealso{\code{\link{fusionName}}, \code{\link{chimeraSeqs}}}
\keyword{utilities}

