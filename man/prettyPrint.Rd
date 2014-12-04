\name{prettyPrint}
\alias{prettyPrint}
\title{A function to represent a list of fSet as a dataframe}
\description{A function reorganizing a list of fSet in a dataframe structure. The dataframe is then saved in a tab delimited file}
\usage{prettyPrint(list, filename, fusion.reads=c("all","spanning"))}
\arguments{
  \item{list}{a list of fSet object}
  \item{filename}{the name of the file in which the dataframe is printed}
   \item{fusion.reads}{it allows to extract spanning reads associated to the SeedCount slot or all the detected reads associate to the RescuedCount}
}
\author{Raffaele A Calogero}
\examples{
#tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hg19")
#fusion.names <- fusionName(tmp)
#tmp1 <- filterList(tmp, type="fusion.names", fusion.names[c(1,3,7)], parallel=FALSE)
#prettyPrint(tmp1, "tmp1.df.txt", fusion.reads="spanning")
}
\keyword{utilities}


