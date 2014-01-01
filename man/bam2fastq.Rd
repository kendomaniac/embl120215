\name{bam2fastq}
\alias{bam2fastq}
\title{A function to extract pair end reads from the bam file generated with subread function}
\description{A function to extract pair end reads from the bam file generated with subread function. The output files are ready to
be used for fusion validaiton with gapfiller}
\usage{bam2fastq(bam, filename="ready4gapfiller",parallel=FALSE)}
\arguments{
  \item{bam}{name of the bam file to be used for PE reads extraction}
  \item{filename}{base name for the PE fastq output data}
  \item{parallel}{option that allow the use of BioParallel package}
}
\value{
       PE fastq files
}
\author{Raffaele A Calogero}
\examples{
#if(require(Rsubread)){
# 	subreadRun(ebwt=paste(find.package(package="chimera"),"/examples/SULF2_ARFGEF2.fa",sep=""), 
#    input1=paste(find.package(package="chimera"),"/examples/mcf7_sample_1.fq",sep=""), 
#    input2=paste(find.package(package="chimera"),"/examples/mcf7_sample_2.fq",sep=""), 
#    outfile.prefix="accepted_hits", alignment="se", cores=1)
#	bam2fastq(bam="accepted_hits.bam", filename="ready4gapfiller", parallel=F)
#}

}
\keyword{utilities}


