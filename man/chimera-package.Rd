\name{chimera-package}
\alias{chimera-package}
\alias{chimera}
\docType{package}
\title{A package for secondary analysis of fusion products}

\description{
The package imports fusion results from tophat-fusion, mapSplice, deFuse, fusionmap, bellerophontes, fusionfinder, fusionhunter, STAR, Rsubread, fusionCatcher.
The package was design to facilitate the characterisation of fusion products events.
Data upload:
outputs for the above indicated fusion detection tools can be imported using \code{\link{importFusionData}} in a list of \code{\link{fSet}}  objects.
\code{\link{fSet-class}} offers various methods to extract information from the fSet objects.
The fusion names can be extracted with \code{\link{fusionName}} function.
The number of reads supporting a fusion event cam be extracted with the \code{\link{supportingReads}} function.

Filtering:
The imported fusion list can be filtered using \code{\link{filterList}}

Annotation:
Oncofuse can be installed in chimera with the function \code{\link{oncofuseInstallation}}
Various information on the fusions location, on structural and functional domains affected by the fusion event as well as a prediction of the putative functional effect of the fusion on the cell can be obtained by using \code{\link{oncofuseRun}}.

\code{\link{chimeraSeqs}} generates the nucleotide sequence of a fusion transcript described in an fSet object. \code{\link{chimeraSeqSet}} does the same but on a list of \code{\link{fSet}} objects.

\code{\link{fusionPeptides}} allows to investigate if the fusion events generate also a fusion at protein level.

\code{\link{subreadRun}} allows to remap reads on the fused transcripts reconstructed with \code{\link{chimeraSeqs}}

Validation
GapFiller is a seed-and-extend local assembler capable to produce (in-silico) longer and highly accurate sequences from a collection of Next Generation Sequencing reads. It can be installed in chimera with the function \code{\link{gapfillerInstallation}}.
The function \code{\link{gapfillerRun}} allows to check if reads mapped by \code{\link{subreadRun}} over a fused transcript generated with \code{\link{chimeraSeqs}} are able to reconstruct by de novo assembly the fusion break-point.

Export
The function \code{\link{prettyPrint}} converts the information stored in a list of \code{\link{fSet}} objects in a dataframe structure that is saved as tab delimited file.
}

\details{
\tabular{ll}{
Package: \tab chimera\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-16\cr
License: \tab Artistic-2.0\cr
}
}

\author{
Raffaele A Calogero
Maintainer: Raffaele A Calogero <raffaele.calogero@unito.it>
}

\references{
Beccuti M, et al. The structure of state-of-art gene fusion-finder algorithms. OA Bioinformatics 2013;1(1):2.
Carrara, M., et al. State of art fusion-finder algorithms are suitable to detect transcription-induced chimeras in normal tissues? BMC bioinformatics 2013;14 Suppl 7:S2.
Carrara, M., et al. State-of-the-Art Fusion-Finder Algorithms Sensitivity and Specificity. BioMed research international 2013;2013:340620.
Shugay, M., et al. Oncofuse: a computational framework for the prediction of the oncogenic potential of gene fusions. Bioinformatics 2013;29(20):2539-2546.
Nadalin, F., et al. GapFiller: a de novo assembly approach to fill the gap within paired reads. BMC bioinformatics 2012;13 Suppl 14:S8.
}

\keyword{ package }

