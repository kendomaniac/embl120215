\name{subreadRun}
\alias{subreadRun}
\title{A function to generate a bam file for fusions coverage evaluation}
\description{A function mapping reads to a chimera sequence set. The bam produced by this remapping on a putative fusion will be used to plot the coverage data for all the fused constructs. The function uses Rsubread aligner.           
}
\usage{subreadRun(ebwt,input1, input2, outfile.prefix="accepted_hits", alignment=c("se","pe"),cores=1)}
\arguments{
  \item{ebwt}{Full path and name of the fasta file of one of the set of fusions of interest, to be used to build the index 
database. The fusion nucleotide sequences can be generated with the function chimeraSeqs}
  \item{input1}{The R1 fastq of a pair-end}
  \item{input2}{The R2 fastq of a pair-end}
  \item{outfile.prefix}{Prefix of the output bam file. Default is accepted_hits}
\item{alignment}{se means that both fastq files from the pair-end run are concatenated, pe means that tophat will use fastq files in pair-end mode}
\item{cores}{Number of cores to be used by the aligner}
}
\value{
    Standard bam file output. The bam file name by default is accepted_hits.bam. 
}
\author{Raffaele A Calogero}
\examples{
    library(Rsubread)
	subreadRun(ebwt=paste(find.package(package="chimera"),"/examples/SULF2_ARFGEF2.fa",sep=""), 
	input1=paste(find.package(package="chimera"),"/examples/mcf7_sample_1.fq",sep=""), 
	input2=paste(find.package(package="chimera"),"/examples/mcf7_sample_2.fq",sep=""), 
	outfile.prefix="accepted_hits", alignment="se", cores=1)
}

\seealso{ \code{\link{chimeraSeqs}}}
\keyword{utilities}


