\name{starInstallation}
\alias{starInstallation}
\title{A function to download STAR}
\description{A function allowing the download and installation of STAR (Dobin et al. Bioinformatics 2012) in chimera package folder. The function also creates soft links in the user bin folder to allow the call of the above mentioned program.            
}
\usage{starInstallation(binDir, os=c("unix","mac"))}
\arguments{
  \item{binDir}{The user bin folder}
  \item{os}{The supported operating systems}
}
\author{Raffaele A Calogero}
\examples{
#starInstallation(binDir="/somewhere/inyourpc/bin", os="mac")
}
\keyword{utilities}


