
###
filterList <- function(x,type=c("supporting.reads","fusion.names", "intronic"),query){
       if(type=="fusion.names"){
	       if(!is.character(query)){
		        cat("\nfiltering by fusion names needs to pass to the method vector of character names\n")
		        return()
		   }
	       tmp <- fusionName(x)
	       loc <- which(tmp %in% query)
       }
       if(type=="supporting.reads"){
	       if(!is.numeric(query)){
		        cat("\nfiltering by supporting reads needs to pass to the method a numerical reads threshold\n")
		        return()
		   }
	       tmp <- supportingReads(x, fusion.reads="spanning")
	       loc <- which(tmp >= query)
       }
       if(type=="intronic"){
	      cat("\n")
          tmp <- sapply(x,.detectIntronic)
          loc <- which(tmp == 0)
       }
       filtered <- x[loc]
       return(filtered)
}
###


.detectIntronic <- function(fset){
  cat(".")
  grl <- fusionGRL(fset)
  #defining the junction as a point object 
  #donor.end
 # start(grl[[1]]) <- end(grl[[1]])
  #acceptor.start
#  end(grl[[2]]) <- start(grl[[2]])
  g1.name <- elementMetadata(grl[[1]])$KnownGene
  g2.name <- elementMetadata(grl[[2]])$KnownGene
  chimera <- paste(g1.name, g2.name, sep=":")
  chr.sym <- as.list(org.Hs.egSYMBOL)
  chimera.tmp <- strsplit(chimera,":")
  if(length(chimera.tmp[[1]]) > 2){return(1)}
  g1 <- chimera.tmp[[1]][1]
  eg1 <- names(chr.sym[which(chr.sym == g1)])
  g2 <- chimera.tmp[[1]][2]
  eg2 <- names(chr.sym[which(chr.sym == g2)])	 
  eg.lst <- list(gene_id=eg1)
  eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
  if(length(eg.trs.n)==0){return(1)}
  #getting only the trs encompassing fusion position
  fusion.trs <- findOverlaps(grl[[1]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
  if(is.na(fusion.trs)){return(1)}
  eg.trs.n <- eg.trs.n[fusion.trs]
  tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
  tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
  tmp.gene1 <- NULL
  for(i in 1:length(tmp.tx)){
		tmp.seq <- .buildFusion(type="donor.end", grl, tmp.tx[i])
        tmp.gene1 <- c(tmp.gene1, tmp.seq$seq)
        if(!is.na(tmp.seq$intron.location)){
	           return(1)
        }
  }
  names(tmp.gene1) <- tmp.name
  #
  eg.lst <- list(gene_id=eg2)
  eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
  if(length(eg.trs.n)==0){return(1)}
  fusion.trs <- findOverlaps(grl[[2]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
  if(is.na(fusion.trs)){return(1)}
  eg.trs.n <- eg.trs.n[fusion.trs]
  tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
  tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
  tmp.gene2 <- NULL
  for(i in 1:length(tmp.tx)){
		  tmp.seq <- .buildFusion(type="acceptor.start", grl, tmp.tx[i])
          tmp.gene2 <- c(tmp.gene2, tmp.seq$seq)
          if(!is.na(tmp.seq$intron.location)){
	           return(1)
        }
  }
  return(0)
}





supportingReads <- function(list, fusion.reads=c("all","spanning"), parallel=F){
	tmp <- list
	if(parallel){
	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
         p <- MulticoreParam()
         if(fusion.reads=="all"){
                nsr <- bplapply(tmp, function(x) x@fusionInfo$RescuedCount, BPPARAM=p)
                nsr <- as.numeric(unlist(nsr))
         }else if(fusion.reads=="spanning"){
	            nsr <- bplapply(tmp, function(x) x@fusionInfo$SeedCount, BPPARAM=p)
	            nsr <- as.numeric(unlist(nsr))
	     }
    }else{
	    if(fusion.reads=="all"){ 
                nsr <- sapply(tmp, function(x) x@fusionInfo$RescuedCount)
		        nsr <- as.numeric(nsr)
	    }else if(fusion.reads=="spanning"){ 
			    nsr <- sapply(tmp, function(x) x@fusionInfo$SeedCount)
			    nsr <- as.numeric(nsr)
	    }
    }
	return(nsr)	
}
###
fusionName <- function(list, parallel=F){ 
	tmp <- list
	if(parallel){
	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
         p <- MulticoreParam()

         g1 <- bplapply(tmp, function(x) x@fusionLoc[[1]]@elementMetadata$KnownGene, BPPARAM=p)
         g1 <- as.character(unlist(g1))
		 g2 <- bplapply(tmp, function(x) x@fusionLoc[[2]]@elementMetadata$KnownGene, BPPARAM=p)
		 g2 <- as.character(unlist(g2))
		 g1g2 <- paste(g1,g2, sep=":")
    }else{ 
        g1 <- sapply(tmp, function(x) x@fusionLoc[[1]]@elementMetadata$KnownGene)
		g2 <- sapply(tmp, function(x) x@fusionLoc[[2]]@elementMetadata$KnownGene)
		g1g2 <- paste(g1,g2, sep=":")
    }
	return(g1g2)
}
####
starReads <- function(fusion.report, parallel=F){
	        if(parallel){ 
		         require(BiocParallel) || stop("\nMission BiocParallel library\n")
	             p <- MulticoreParam()
	        }
		    report <- read.table(fusion.report, sep="\t", header=F)
		    names(report) <- c("gene1.chr","gene1.start", "gene1.strand", "gene2.chr","gene2.start", "gene2.strand","junction.type","repeat.left","repeat.right","reads.name","first.aligned.read1","read1.cigar","first.aligned.read2","read2.cigar")
            #reads GA
            r1.gr <- GRanges(seqnames = as.character(report$gene1.chr), 
                              ranges = IRanges(start = as.numeric(report$first.aligned.read1), 
                              end= as.numeric(report$first.aligned.read1)),  cigar=as.character(report$read1.cigar))
            r2.gr <- GRanges(seqnames = as.character(report$gene2.chr), 
                              ranges = IRanges(start = as.numeric(report$first.aligned.read2), 
                              end= as.numeric(report$first.aligned.read2)),  cigar=as.character(report$read2.cigar))
		    reads.gr <- GRangesList(r1.gr, r2.gr)
		    return(reads.gr)
}            


