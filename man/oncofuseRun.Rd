\name{oncofuseRun}
\alias{oncofuseRun}
\title{A function to annotate fusions with Oncofuse. Oncofuse is a naive bayesian classifier. Its goal is to identify those fusion sequences with higher probability of being driver than passenger events}
\description{A list of fSet objects can be annotated using the Oncofuse java application.           
}
\usage{oncofuseRun(listfSet, tissue=c("EPI","HEM","MES","AVG"), org=c("hg19","hg38"), threads=1, plot=FALSE)}
\arguments{
  \item{listfSet}{A list of fSets}
  \item{tissue}{Type of tissue in which the fusion was detected. EPI epithelial, HEM ematological, MES mesenchimal, AVG average expression}
  \item{org}{Supported genome assembly version. IMPORTANT: the genome version used for the fusion detection and for the oncofuse analysis need to be same}
  \item{threads}{number of threads that Oncofuse will use}
  \item{plot}{plotting the expression gain score versus the passenger mutation probability. For mor einfo please see  Shugay et al. Bioinformatics 2013:29,2539}
}
\value{
    The output is a dataframe containing the output of Oncofuse. For more info on Oncofuse please see Shugay et al. Bioinformatics 2013, 29, 2539. 
}
\author{Raffaele A Calogero}
\examples{
#tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hg19")
#installOncofuse()
#of.out <- oncofuseRun(tmp, tissue="EPI")
}
\keyword{utilities}


