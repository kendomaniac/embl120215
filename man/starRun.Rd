\name{starRun}
\alias{starRun}
\title{A function to generate a bam file for fusions coverage evaluation}
\description{A function mapping reads to a chimera sequence set. The bam produced by this remapping on a putative fusion will be used to plot the coverage data for all the fused constructs. The function assumes that STAR is installed and located in the path.            
}
\usage{starRun(input1, input2, cores=1, star= "STAR", samtools="samtools", fa, alignment=c("se","pe"), chimSegmentMin=0, chimJunctionOverhangMin=20)}
\arguments{
  \item{input1}{The R1 fastq of a pair-end}
  \item{input2}{The R2 fastq of a pair-end}
  \item{cores}{number of cores to be used by tophat with program name, e.g. /somewhere/tophat}
  \item{star}{full path of STAR with program name, e.g. /somewhere/STAR}
  \item{samtools}{full path of samtools}
  \item{fa}{full path and name of the fasta file of one of the set of fusions of interest, to be used to build the STAR 
database. The fusion nucleotide sequences was generated with the function chimeraSeqs}
  \item{alignment}{se means that both fastq files from the pair-end run are concatenated, pe means that tophat will use fastq files in pair-end mode}
  \item{chimSegmentMin}{STAR fusion parameter, see STAR manual}
  \item{chimJunctionOverhangMin}{STAR fusion parameter, see STAR manual}

}
\value{
    The function create a folder called chimeraDB_time, where time is the time when the folder was created. STAR output will be located in the folder output_time, where time is the time when the folder was created. The bam file of interest is accepted_hits.bam. 
}
\author{Raffaele A Calogero}
\examples{
#starRun(input1=paste(find.package(package="chimera"),"/examples/mcf7_sample_1.fq",sep=""), input2=paste(find.package(package="chimera"),"/examples/mcf7_sample_2.fq",sep=""), fa=paste(find.package(package="chimera"),"/examples/SULF2_ARFGEF2.fa",sep=""), alignment="se", output=getwd(), cores=1)
}

\seealso{ \code{\link{chimeraSeqs}}}
\keyword{utilities}


