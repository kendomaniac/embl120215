\name{tophatInstallation}
\alias{tophatInstallation}
\title{A function to download tophat, bowtie and samtools}
\description{A function allowing the download and installation of tophat, bowtie and samtools in chimera package folder. The function also creates soft links in the user bin folder to allow the call of the above mentioned programs.            
}
\usage{tophatInstallation(binDir, os=c("unix","mac"))}
\arguments{
  \item{binDir}{The user bin folder}
  \item{os}{The supported operating systems}
}
\author{Raffaele A Calogero}
\examples{
#tophatInstallation(binDir="/somewhere/inyourpc/bin", os="mac")
}
\keyword{utilities}


