\name{plotChr}
\alias{plotChr}
\title{A function plot a genes to chromosome association}
\description{A function to plot genes to chromosome association}
\usage{plotChr(de.df, org, genome, plot=c("bars","chart"))}
\arguments{
  \item{de.df}{A dataframe of differentially expressed genes having columns log2FC, SYMBOLS, ENTREZID, GENE DESCRIPTION}
  \item{org}{Organism database, e.g. Homo.sapiens}
  \item{genome}{Genome version, e.g. mm10}
    \item{plot}{type of output plot. bars provides a plot of the number of genes in any chromosomes. Chart provides a plot of the position of the genes on the chromosomes}
}
\author{Raffaele A Calogero}
\examples{
#de.lst <- plotChrs(de.df=de.selected, org="Mus.musculus", genome="mm10", plot="bar")
#de.lst <- plotChrs(de.df=de.selected, org="Mus.musculus", genome="mm10", plot="chart")
}
\keyword{utilities}


