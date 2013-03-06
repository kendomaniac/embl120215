
###
fusionPeptides <- function(fset, which.isoform=1, donor.up=200, acceptor.down=200){
    rna <- fusionRNA(fset)
    rna <- rna[which.isoform] 
    rna.name <- names(rna)
	gr <- fusionGRL(fset)	
	g1.name <- elementMetadata(gr[[1]])$KnownGene
	g2.name <- elementMetadata(gr[[2]])$KnownGene
	chimera <- paste(g1.name, g2.name, sep=":")
	chr.sym <- as.list(org.Hs.egSYMBOL)
    chimera.tmp <- strsplit(chimera,":")
    if(length(chimera.tmp[[1]])>2){
		     cat("\nAt list one of the fusion element is not an annotated gene.\nExtraction of fusion regions for this specific situation is not implemented, yet.\n")
	         return()
	}
	rna.namelst <- strsplit(rna.name, ":")
	length.donor <-strsplit(rna.namelst[[1]], "-")
	p1.name <- length.donor[[1]][1]
	p2.name <- length.donor[[2]][1]
	chimera.name <- paste(p1.name, p2.name, sep=":")
	length.donor <- as.numeric(length.donor[[1]][2])
	p1.cdna <- rna[[1]][1:length.donor]
    cdna.end <- trunc(length(p1.cdna)/3) * 3
    left1.nts <- length(p1.cdna) - cdna.end
    if(left1.nts == 0) {
	   cat("\nCDS of gene1 is on frame 1\n")
	   frame1 <- 0
	    printed.frame1 <- 1
	}
    if(left1.nts == 1) {
	   cat("\nCDS of gene1 is on frame 2\n")
	   frame1 <- 1
	   printed.frame1 <- 2
    }
    if(left1.nts == 2) {
	   cat("\nCDS of gene1 is on frame 3\n")
	   frame1 <- 2
	   printed.frame1 <- 3
	}
	p2.cdna <- rna[[1]][(length.donor + 1):length(rna[[1]])]
    cdna.start <- trunc(length(p2.cdna)/3) * 3
    left2.nts <- length(p2.cdna) - cdna.start
    if(left2.nts == 0) {
	   cat("\nCDS of gene2 is on frame 1\n")
	   frame2 <- 0
	    printed.frame2 <- 1
	}
    if(left2.nts == 1) {
	   cat("\nCDS of gene2 is on frame 2\n")
	   frame2 <- 1
	    printed.frame2 <- 2
    }
    if(left2.nts == 2) {
	   cat("\nCDS of gene2 is on frame 3\n")
	   frame2 <- 2
	   printed.frame2 <- 3
	}
	tmp1 <- as.character(p1.cdna)
	tmp2 <- as.character(p2.cdna)
	sub1<- substring(tmp1,(nchar(tmp1) - donor.up), nchar(tmp1))
	sub2<- substring(tmp2,1, acceptor.down)
	tmp.ga <- fusionGA(fset)
	tmp.ga <- tmp.ga[which(as.character(seqnames(tmp.ga))== chimera.name)]
	tmp.ga <- tmp.ga[intersect(which(start(tmp.ga) < length(p1.cdna)), which(end(tmp.ga) > length(p1.cdna)))]
	if(sum(left1.nts, left2.nts)== 3) {
		fusion.status <- "In frame"
	}else{
		fusion.status <- "Not in frame"
	}
	
    return(list(selected.alignment=which.isoform, 
	            transcript1=p1.cdna, 
	            pep1=translate(p1.cdna[1:(length(p1.cdna)- frame1)]), 
	            frame.pep1=printed.frame1,
                transcript2=p2.cdna, 
                pep2=translate(p2.cdna[1:length(p2.cdna)]),
                frame.pep2=printed.frame2,
                fusion.status=fusion.status,
                validation.seq=paste(sub1,sub2,collapse=""), 
                junction.ga=tmp.ga)
     )
}
