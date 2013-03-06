\name{fusionPeptides}
\alias{fusionPeptides}
\title{A function associating to a fusion the peptides involved.}
\description{A function associate the donor and the acceptor peptides involved in the fusion.             
}
\usage{fusionPeptides(fset, which.isoform=1, donor.up=200, acceptor.down=200)}
\arguments{
  \item{fset}{An fSet object, encompassing a fusion including two annotated genes}
  \item{which.isoform}{Since each gene can be represented by multiple transcripts, one of the available transcripts need to be selected}
\item{donor.up}{The number of nucleotides upstream to the donor end, to be used to extract the DNA fragment for PCR validation}
\item{acceptor.down}{The number of nucleotides downstream to the acceptor start, to be used to extract the DNA fragment for PCR validation}
}

\value{
A list encompassing:
\item{selected.alignment}{the alignment used to define the position of the donor end on the transcript} 
\item{transcrip1}{the nucleotide sequence of the donor fusion transcript} 
\item{pep1}{the peptide sequence encoded by the donor transcript}
\item{frame.pep1}{The frame of the encoded peptide}
\item{transcrip2}{the nucleotide sequence of the acceptor fusion transcript} 
\item{selected.alignment2}{the alignment used to define the position of the acceptor start on the transcript} 
\item{pep2}{the peptide sequence encoded by the acceptor transcript}
\item{frame.pep2}{The frame of the encoded peptide}
\item{fusion.status}{if peptides are in frame on the basis of the fused mRNAs. IMPORTANT: this information os only based on fusion junction information, if indels are present results might be different. sequencign of the validation.seq, see below are mandatory to confirm the fusion event.}
\item{validation.seq}{The DNA fragment encompassing the fusion to be used to design PCR primers for sanger sequencing}
\item{junction.ga}{The GappedAlignment object describing the reads encompassing the fusion junction. If no genomic regions are indicated the reads spanning over the junction are missed.}
}
\author{Raffaele A Calogero}
\examples{
#	#creating a fusion report from output of fusionMap
	tmp <- importFusionData("fusionmap", paste(find.package(package="chimera"),"/examples/mcf7.FMFusionReport", sep=""), org="hs")
	#extracting the fusions names
	fusion.names <- fusionName(tmp)
	fusion.names
	#extracting the fSet object ofr one of the fusions
	myset <- tmp[[13]]
	#constructing the fused sequence(s)
	trs <- chimeraSeqs(myset, type="transcripts")
	#adding the sequences to the fSet object
	myset <- addRNA(myset , trs)
	#extracting sequences from an fSet object
	tmp.seq <- fusionRNA(myset)
	#adding reads mapped on the fusion generated using tophatRun function
	myset <- addGA(myset, paste(path.package(package="chimera"),"/examples/mcf7_trs_accepted_hits.bam",sep=""))

   mypeps <- fusionPeptides(fset=myset, which.isoform=1, donor.up=200, acceptor.down=200)
}
\keyword{utilities}

