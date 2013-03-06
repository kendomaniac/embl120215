\name{filterList}
\alias{filterList}
\title{A function to filter a list of fSet objects}
\description{A function filtering a list of fSet objects on the basis of supporting reads or fusion names.}
\usage{filterList(x,type=c("supporting.reads","fusion.names", "intronic"),query)}
\arguments{
  \item{x}{a list of fSet object}
  \item{type}{filtering can be performed on the basis of a minimal number of supporting reads or on the basis of fusion names}
  \item{query}{query is number in case supporting.reads is selected or a vector of fusion names if the case fusion.names is selected. In case type is intronic query is NULL. In the latter case fusion having one of the elements located in an intronic region are discarded}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
fusion.names <- fusionName(tmp)
tmp1 <- filterList(tmp, type="fusion.names", fusion.names[c(1,3,7)])
tmp2 <- filterList(tmp, type="supporting.reads", 2)
}
\keyword{utilities}

