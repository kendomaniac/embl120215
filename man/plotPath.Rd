\name{plotPath}
\alias{plotPath}
\title{A function plot a path}
\description{A function to associate EG and name symbols to rsem-generated-gene-tables, using UCSC as annotation}
\usage{plotPath(my.path, path.db=c("kegg","biocarta","nci","spike","humancyc","panther"), type=c("entrez","symbol"))}
\arguments{
  \item{my.path}{Name of the path of interest}
  \item{path.db}{Name of the pathway to be used for plotting}
  \item{type}{Entrez geneid or Symbol}
}
\author{Raffaele A Calogero}
\examples{
#  samples.dir <- paste(find.package(package="thu12feb25"),"/examples/", sep="")
#  samples <- dir(samples.dir)[grep(".genes.results",dir(samples.dir))]
#  counts <- lapply(paste(samples.dir,samples, sep=""), function(x){
#      tmp <- annotatingGenes(filename=x, org="mm9")
#	  counts.tmp <- tmp$expected_count
#  })
# counts.df <- t(data.frame(matrix(unlist(counts), nrow=length(samples), byrow=T)))
# tmp <- annotatingGenes(filename=paste(samples.dir,samples[1], sep=""), org="mm9")
# dimnames(counts.df) <- list(as.character(tmp$Symbol),gsub(".genes.results","",samples))

#   counts.df1 <- counts.df[rowSums(counts.df) > 10,]
#   counts.df1 <- as.matrix(counts.df1)
#converting in uppercase gene symbols
#   dimnames(counts.df1)[[1]] <- toupper(dimnames(counts.df1)[[1]])
#   deg <- DEGraph(counts.df1, conditions, panther, type="RNASeq", IDs = "symbol", norm.method="DESeq2")
#   names(which(deg$res[,1] <= 0.05))
#3
#   names(panther)[which(names(panther)%in%names(which(deg$res[,1] <= 0.05)))]

#	plotPath(my.path="mRNA splicing", path.db="panther", type="symbol")

}
\keyword{utilities}


