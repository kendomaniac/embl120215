\name{annotatingGenes}
\alias{annotatingGenes}
\title{A function to associate EG and name symbols to rsem-generated-gene-tables, using UCSC as annotation}
\description{A function to associate EG and name symbols to rsem-generated-gene-tables, using UCSC as annotation}
\usage{annotatingGenes(filename, org=c("hg19","mm9"), truncating.expected.counts=FALSE)}
\arguments{
  \item{filename}{Name of the file generated with rsem with the suffix .genes.results}
  \item{org}{A tag indicating the organism and its genomic version}
  \item{truncating.expected.counts}{the expected counts calculated by rsem are nont integer. 
  If this parameter is set to TRUE the expected counts are truncated removing the decimals and converted in integers}
}
\author{Raffaele A Calogero}
\examples{
#	dir <- dir()
#	dir <- dir[grep(".genes.results", dir)]
#	tmp <- annotatingGenes(filename=dir[1], org="mm9")
}
\keyword{utilities}


