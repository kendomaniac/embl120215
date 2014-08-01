\name{bam2fastq}
\alias{bam2fastq}
\title{A function to extract pair end reads from the bam file generated with subreadRun function}
\description{A function to extract pair end reads from the bam file generated with subread function. The output files are ready to
be used for fusion validation with gapfiller}
\usage{bam2fastq(bam, filename="ready4gapfiller",ref,parallel=FALSE)}
\arguments{
  \item{bam}{name of the bam file to be used for PE reads extraction}
  \item{filename}{base name for the PE fastq output data}
  \item{ref}{name of the fusion sequence that was used as reference}
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
#   ref.name <- names(readDNAStringSet(paste(find.package(package="chimera"),"/examples/SULF2_ARFGEF2.fa",sep=""), format="fasta"))
#	bam2fastq(bam="accepted_hits.bam", filename="ready4gapfiller", ref=ref.name, parallel=F)
#}

}
\keyword{utilities}


