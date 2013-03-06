\name{tophatRun}
\alias{tophatRun}
\title{A function to generate a bam file for fusions coverage evaluation}
\description{A function mapping reads to a chimera sequence set. The bam produced by this remapping on a putative fusion will be used to plot the coverage data for all the fused constructs. The function assumes that tophat is installed and located in the path.            
}
\usage{tophatRun(input1, input2, output,cores=1, bowtie= "bowtie", tophat= "tophat",ebwt=paste(getwd(),"mychimera.fa",sep="/"),alignment=c("se","pe"))}
\arguments{
  \item{input1}{The R1 fastq of a pair-end}
  \item{input2}{The R2 fastq of a pair-end}
  \item{output}{Folder in which tophat will generate the output}
  \item{cores}{number of cores to be used by tophat with program name, e.g. /somewhere/tophat}
  \item{bowtie}{full path of bowtie with program name, e.g. /somewhere/bowtie}
  \item{tophat}{full path of tophat}
  \item{ebwt}{full path and name of the fasta file of one of the set of fusions of interest, to be used to build the bowtie 
database. The fusion nucleotide sequences was generated with the function chimeraSeqs}
\item{alignment}{se means that both fastq files from the pair-end run are concatenated, pe means that tophat will use fastq files in pair-end mode}

}
\value{
    TopHat standard output. The bam file of interest is accepted_hits.bam. The bam file will be then loaded in the slot fusionsLoc of the fSetSummary object from which fusions were retrieved. 
}
\author{Raffaele A Calogero}
\examples{
#tophatRun(input1=paste(find.package(package="chimera"),"/examples/mcf7_sample_1.fq",sep=""), input2=paste(find.package(package="chimera"),"/examples/mcf7_sample_2.fq",sep=""), ebwt=paste(find.package(package="chimera"),"/examples/SULF2_ARFGEF2.fa",sep=""), alignment="se", output=getwd(), cores=1)
}

\seealso{ \code{\link{chimeraSeqs}}}
\keyword{utilities}

