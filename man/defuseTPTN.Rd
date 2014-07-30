\name{defuseTPTN}
\alias{defuseTPTN}
\title{A function that genererate a list of fSet objects encopassing 60 experimentally validated fusions and 61 false fusions}
\description{The function uses the file defuse_TP_FP.txt located in the examples folder of chimera to 
generate a list of fSet object encompassing 60 TP fusions and 61 TN fusions extracted from supplementary table 8
of the deFuse paper McPherson et al. (2011) PLoS Comput Biol 7(5): e1001138.
The fset object can be used to evaluate the efficacy of filtering approaches available in chimera package. 
}
\usage{defuseTPTN()
}
\value{
       An list of fSet objects: Experimentally validated fusions  1-60, fish validated fusions 1-14, false fusions 61-121.
}
\author{Raffaele A Calogero}
\examples{
tptn <-defuseTPTN()
}
\keyword{utilities}


