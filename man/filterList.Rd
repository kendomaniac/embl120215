\name{filterList}
\alias{filterList}
\title{A function to filter a list of fSet objects}
\description{A function filtering a list of fSet objects on the basis of supporting reads or fusion names.}
\usage{filterList(x,type=c("supporting.reads","fusion.names", "intronic", "annotated.genes", "read.through"),query, parallel=F)}
\arguments{
  \item{x}{a list of fSet object}
  \item{type}{filtering can be performed on the basis of
              supporting.read: a minimal number of supporting reads,
              fusion.names: a vector list of user defined fusion names,
              intronic: only fusion encompassing exon data are retained
              annotated.genes: only fusion encompassing two annotated genes are retained
              read.through: only fusion different gene names are retained}
  \item{query}{query is number in case supporting.reads is selected or a vector of fusion names if the case fusion.names is selected. In case type is intronic query is NULL. In the latter case fusion having one of the elements located in an intronic region are discarded}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
fusion.names <- fusionName(tmp)
tmp1 <- filterList(tmp, type="fusion.names", fusion.names[c(1,3,7)], parallel=F)
tmp2 <- filterList(tmp, type="supporting.reads", 2, parallel=F)
#tmp3 <- filterList(tmp, type="intronic")
#tmp4 <- filterList(tmp, type="annotated.genes", parallel=F)
#tmp5 <- filterList(tmp, type="read.through", parallel=F)
}
\keyword{utilities}

