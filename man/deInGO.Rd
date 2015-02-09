\name{deInGO}
\alias{deInGO}
\title{A function plot a genes to chromosome association}
\description{A function to plot genes to chromosome association}
\usage{deInGO(go.term, de.universe, org=c("Mus.musculus","Homo.sapiens"))}
\arguments{
  \item{go.term}{A dataframe of differentially expressed genes having columns log2FC, SYMBOLS, ENTREZID, GENE DESCRIPTION}
  \item{de.universe}{A vector containing 0 for all genes not DE and 1 for the DE. the vector is named with the gene symbols}
  \item{org}{Organism database, e.g. Homo.sapiens}
}
\author{Raffaele A Calogero}
\examples{
#deInGO(go.term="GO:0030154", de.universe=de.universe, org="Mus.musculus")
}
\keyword{utilities}


