\name{filterList}
\alias{filterList}
\title{A function to filter a list of fSet objects}
\description{A function filtering a list of fSet objects on the basis of various parameters.}
\usage{filterList(x,type=c("supporting.reads","fusion.names", "intronic", "annotated.genes", "read.through","oncofuse"),oncofuse.output=NULL, query=NULL, oncofuse.type=c("g5CDS", "g3CDS", "g5g3CDS", "g5exon", "g3exon", "g5g3exon","passenger.prob", "expression.gain"), parallel=FALSE)}
\arguments{
  \item{x}{a list of fSet object}
  \item{type}{filtering can be performed on the basis of
              supporting.read: a minimal number of supporting reads,
              fusion.names: a vector list of user defined fusion names,
              intronic: only fusion encompassing exon data are retained
              annotated.genes: only fusion encompassing two annotated genes are retained
              read.through: only fusion with different gene names are retained
			  oncofuse: using the output of oncofuseRun various filtering option can be applied using oncofuse.type and oncofuse.threshold}
  \item{oncofuse.output}{The output generated with oncofuseRun}
  \item{query}{query is a number. 
  In case supporting.reads is selected or a vector of fusion names if the case fusion.names is selected. 
  In case type is intronic query is NULL. In the latter case fusion having one of the elements located in an intronic region are discarded. 
  In the case of oncofuse.type equal to passenger.prob or expression.gain query is the threshold number to be used for the filtering}
  \item{oncofuse.type}{This refers to the filtering based on the output of oncofuse that can be generated using the oncofuseRun function. g5CDS selects only fusions having the breakpoint in 5' end gene CDS, g3CDS selects only fusions having the breakpoint in 3' end gene CDS, g5g3CDS selects only fusions having the breakpoint in both gene CDSs, g5exon selects only fusions having the breakpoint in an exon of 5' end gene, g3exon selects only fusions having the breakpoint in an exon of 3' end gene, g5g3exon selects only fusions having the breakpoint in an exon of both genes. In the case of oncofuse.type equal to driver.prob the filter will use the probability of the fusion to be a tumor driver. In the case of oncofuse.type equal to expression.gain the filter will use the score that suggests the presence of a gain in expression}
  \item{parallel}{option to run a parallel version of the function, default FALSE}
}
\author{Raffaele A Calogero}
\examples{
tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
fusion.names <- fusionName(tmp)
tmp1 <- filterList(tmp, type="fusion.names", fusion.names[c(1,3,7)], parallel=FALSE)
tmp2 <- filterList(tmp, type="supporting.reads", query=2, parallel=FALSE)
#tmp3 <- filterList(tmp, type="intronic")
#tmp4 <- filterList(tmp, type="annotated.genes", parallel=FALSE)
#tmp5 <- filterList(tmp, type="read.through", parallel=FALSE)
#csdf.of <- oncofuseRun(csdf.e, tissue="EPI")
#tmp6 <- filterList(csdf.e[1:100],oncofuse.output=csdf.of, type="oncofuse", oncofuse.type="g5g3CDS", parallel=TRUE)
#tmp7 <- filterList(csdf.e[1:100],oncofuse.output=csdf.of, type="oncofuse", oncofuse.type="passenger.prob", query=0.1, parallel=TRUE)
}
\keyword{utilities}


