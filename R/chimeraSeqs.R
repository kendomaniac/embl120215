#functions:
#importFusionData
#.fmImport import fusions from FusionMap
#.ffImport import fusions from FusionFinder
#.fhImport import fusions from FusionHunter
#.msImport import fusions from MapSplice
#.thfImport import fusions from TopHat-fusion
#.dfImport import fusions from deFuse
#.chimeraSeqs generates the fusion nucleotide sequence
#tophatInstallation
#tophatRun
#plotCoverage
#filterList
#supportingReads
#fusionName
chimeraSeqs <- function(fset, extend=1000, type="transcripts"){
	grl <- fusionGRL(fset)
	#defining the junction as a point object 
    #donor.end
 #   start(grl[[1]]) <- end(grl[[1]])
    #acceptor.start
#    end(grl[[2]]) <- start(grl[[2]])
	g1.name <- elementMetadata(grl[[1]])$KnownGene
	g2.name <- elementMetadata(grl[[2]])$KnownGene
	chimera <- paste(g1.name, g2.name, sep=":")
	chr.sym <- as.list(org.Hs.egSYMBOL)
    chimera.tmp <- strsplit(chimera,":")
    if(length(chimera.tmp[[1]])==2){
	 g1 <- chimera.tmp[[1]][1]
	 eg1 <- names(chr.sym[which(chr.sym == g1)])
	 g2 <- chimera.tmp[[1]][2]
	 eg2 <- names(chr.sym[which(chr.sym == g2)])	 
	 if(type=="transcripts"){
		eg.lst <- list(gene_id=eg1)
		eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
		if(length(eg.trs.n)==0){
			cat("\nERROR: The Entrez gene id returns an empty GenomicRange object\n")
			return(eg.trs.n) 
		}
		#getting only the trs encompassing fusion position
		fusion.trs <- findOverlaps(grl[[1]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
		if(is.na(fusion.trs)){
			cat("\nThe location of the fusion does not seems to fit the location of the donor end gene\n")
			return(GRangesList(donor.fusion=grl[[1]], donor.end.gene=eg.trs.n))
		}
		eg.trs.n <- eg.trs.n[fusion.trs]
		tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
		tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
		tmp.gene1 <- NULL
		for(i in 1:length(tmp.tx)){
		  tmp.seq <- .buildFusion(type="donor.end", grl, tmp.tx[i])
          tmp.gene1 <- c(tmp.gene1, tmp.seq$seq)
          if(!is.na(tmp.seq$intron.location)){
	           cat(paste("\nIntron ",tmp.seq$intron.location," was used as donor end\n",sep=""))
          }
          
        }
        names(tmp.gene1) <- tmp.name
        #
		eg.lst <- list(gene_id=eg2)
		eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
		if(length(eg.trs.n)==0){
			cat("\nERROR: The Entrez gene id returns an empty GenomicRange object\n")
			return(eg.trs.n) 
		}		
		fusion.trs <- findOverlaps(grl[[2]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
		if(is.na(fusion.trs)){
			cat("\nThe location of the fusion does not seems to fit the location of the acceptor start gene\n")
			return(GRangesList(acceptor.fusion=grl[[2]], donor.end.gene=eg.trs.n))
		}		
		eg.trs.n <- eg.trs.n[fusion.trs]
		tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
		tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
		tmp.gene2 <- NULL
		for(i in 1:length(tmp.tx)){
		    tmp.seq <- .buildFusion(type="acceptor.start", grl, tmp.tx[i])
            tmp.gene2 <- c(tmp.gene2, tmp.seq$seq)
            if(!is.na(tmp.seq$intron.location)){
	           cat(paste("\nIntron ",tmp.seq$intron.location," was used as acceptor start\n",sep=""))
            }
        }
        names(tmp.gene2) <- tmp.name
        fusions <- NULL
        fusions.names <- NULL
        for(i in 1: length(tmp.gene1)){
	      	for(j in 1: length(tmp.gene2)){
                  fusions <- c(fusions, paste(as.character(tmp.gene1[[i]]), as.character(tmp.gene2[[j]]),sep="",collapse=""))
                  fusions.names <- c(fusions.names, paste(paste(names(tmp.gene1[i]),length(tmp.gene1[[i]]),sep="-"), paste(names(tmp.gene2[j]), length(tmp.gene2[[j]]), sep="-"), sep=":"))
	        }
        }
        fusions <- DNAStringSet(unlist(fusions))
        names(fusions) <- fusions.names
        
        return(fusions)
      }else if(type=="gene"){
	     cat("\nNot implemented\n")
	     return()
	  }
   }else if(length(chimera.tmp[[1]])==3){
      tmp <- grep("chr",chimera.tmp[[1]])
      if(type=="transcripts"){
        if(tmp == 2){
	         g1 <- chimera.tmp[[1]][1]
	         eg1 <- names(chr.sym[which(chr.sym == g1)])
	         eg.lst <- list(gene_id=eg1)
		     eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
		     #getting only the trs encompassing fusion position
		     fusion.trs <- findOverlaps(grl[[1]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
		     eg.trs.n <- eg.trs.n[fusion.trs]
		     tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
		     tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
		     tmp.gene1 <- NULL
		     for(i in 1:length(tmp.tx)){
               tmp.gene1 <- c(tmp.gene1, .buildFusion(type="donor.end", grl, tmp.tx[i]))
             }
             names(tmp.gene1) <- tmp.name #list of seqs

             grl2 <- new("GRangesList")
	         cat("\nAcceptor gene start does not fit the expected start on the detected fusion\nFusion is located outside gene exons. Fusion sequence is assembled using:")
             end(grl[[2]]) <- end(grl[[2]]) + extend
             cat(paste(as.character(seqnames(grl[[2]])), sep=":",paste(start(grl[[2]]), end(grl[[2]]),sep="-")),"\n")
	         grl2[[1]] <- grl[[2]]
	         names(grl2) <- paste(as.character(seqnames(grl[[2]])), sep=":",paste(start(grl[[2]]), end(grl[[2]]),sep="-"))
	         acceptorGeneList <- list()
	         z <- 1
	         for(i in 1:length(grl2)){
		        gr2.tmp <- grl2[i]
		        acceptorGene <- NULL
		        for(j in 1:length(gr2.tmp[[1]])){
		           acceptorGene[j] <- getSeq(Hsapiens, seqnames(gr2.tmp[[1]][j]), start=start(ranges(gr2.tmp[[1]][j])), end=end(ranges(gr2.tmp[[1]][j])), as.character=TRUE)
		        }
		        acceptorGeneList[[z]] <- paste(acceptorGene, collapse="")
		        z <- z+1
	          }
	          names(acceptorGeneList) <- names(grl2)
	          fusions <- NULL
	          fusions.names <- NULL
              for(i in 1: length(tmp.gene1)){
	            	for(j in 1: length(acceptorGeneList)){
                       fusions <- c(fusions, paste(as.character(tmp.gene1[[i]]), as.character(acceptorGeneList[[j]]),sep="",collapse=""))
                       fusions.names <- c(fusions.names, paste(paste(names(tmp.gene1[i]),length(tmp.gene1[[i]]),sep="-"), paste(names(acceptorGeneList[j]), extend,sep="-"), sep=":"))
	                }
              }
              fusions <- DNAStringSet(unlist(fusions))
              names(fusions) <- fusions.names
	          return(fusions)
	     }else if(tmp == 1) {
	       grl1 <- new("GRangesList")
	       cat("\nDonor gene end does not fit the expected end on the detected fusion.\nFusion is located outside gene exons. Fusion sequence is assembled using:")
           start(grl[[1]]) <- start(grl[[1]]) - extend
           cat(paste(as.character(seqnames(grl[[1]])), sep=":",paste(start(grl[[1]]), end(grl[[1]]),sep="-")),"\n")
	       grl1[[1]] <- grl[[1]]
	       names(grl1) <- paste(as.character(seqnames(grl[[1]])), sep=":",paste(start(grl[[1]]), end(grl[[1]]),sep="-"))
	       #extracting seq
           donorGeneList <- list()
           z <- 1
           for(i in 1:length(grl1)){
	          gr1.tmp <- grl1[i]
	          donorGene <- NULL
	          for(j in 1:length(gr1.tmp[[1]])){
	              donorGene[j] <- getSeq(Hsapiens, seqnames(gr1.tmp[[1]][j]), start=start(ranges(gr1.tmp[[1]][j])), end=end(ranges(gr1.tmp[[1]][j])), as.character=TRUE)
	          }
	          donorGeneList[[z]] <- paste(donorGene, collapse="")
	          z <- z+1
           }
           names(donorGeneList) <- names(grl1)
           ###
	       g2 <- chimera.tmp[[1]][3]
	       eg2 <- names(chr.sym[which(chr.sym == g2)])
	       eg.lst <- list(gene_id=eg2)
		   eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id", "tx_name"))
		   #getting only the trs encompassing fusion position
		   fusion.trs <- findOverlaps(grl[[2]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
		   eg.trs.n <- eg.trs.n[fusion.trs]
		   tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
		   tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
		   tmp.gene2 <- NULL
		   for(i in 1:length(tmp.tx)){
              tmp.gene2 <- c(tmp.gene2, .buildFusion(type="acceptor.start", grl, tmp.tx[i]))
           }
           names(tmp.gene2) <- tmp.name #list of seqs
	       fusions <- NULL
	       fusions.names <- NULL
           for(i in 1: length(tmp.gene2)){
	        	for(j in 1: length(donorGeneList)){
                    fusions <- c(fusions, paste(as.character(donorGeneList[[j]]), as.character(tmp.gene2[[i]]) ,sep="",collapse=""))
                    fusions.names <- c(fusions.names, paste(paste(names(donorGeneList[j]), extend, sep="-"), paste(names(tmp.gene2[j]), length(tmp.gene2[[j]]), sep="-"), sep=":"))
	            }
           }
           fusions <- DNAStringSet(unlist(fusions))
           names(fusions) <- fusions.names
           return(fusions)
        }
     }else if(type=="gene"){
	       cat("\nNot implemented, yet\n")
	       return()
	 }
   }else if(length(chimera.tmp[[1]])==4){
	    if(type=="transcripts"){
	       grl1 <- new("GRangesList")
	       cat("\nDonor gene end does not fit the expected end on the detected fusion.\nFusion is located outside gene exons. Fusion sequence is assembled using:")
           start(grl[[1]]) <- start(grl[[1]]) - extend
           cat(paste(as.character(seqnames(grl[[1]])), sep=":",paste(start(grl[[1]]), end(grl[[1]]),sep="-")),"\n")
	       grl1[[1]] <- grl[[1]]
	       names(grl1) <- paste(as.character(seqnames(grl[[1]])), sep=":",paste(start(grl[[1]]), end(grl[[1]]),sep="-"))
	       #extracting seq
           donorGeneList <- list()
           z <- 1
           for(i in 1:length(grl1)){
	            gr1.tmp <- grl1[i]
	            donorGene <- NULL
	            for(j in 1:length(gr1.tmp[[1]])){
	                donorGene[j] <- getSeq(Hsapiens, seqnames(gr1.tmp[[1]][j]), start=start(ranges(gr1.tmp[[1]][j])), end=end(ranges(gr1.tmp[[1]][j])), as.character=TRUE)
	            }
	            donorGeneList[[z]] <- paste(donorGene, collapse="")
	            z <- z+1
           }
           names(donorGeneList) <- names(grl1)
#
           grl2 <- new("GRangesList")
	       cat("\nAcceptor gene start does not fit the expected start on the detected fusion\nFusion is located outside gene exons. Fusion sequence is assembled using:")
           end(grl[[2]]) <- end(grl[[2]]) + extend
           cat(paste(as.character(seqnames(grl[[2]])), sep=":",paste(start(grl[[2]]), end(grl[[2]]),sep="-")),"\n")
	       grl2[[1]] <- grl[[2]]
	       names(grl2) <- paste(as.character(seqnames(grl[[2]])), sep=":",paste(start(grl[[2]]), end(grl[[2]]),sep="-"))
	       acceptorGeneList <- list()
	       z <- 1
	       for(i in 1:length(grl2)){
		        gr2.tmp <- grl2[i]
		        acceptorGene <- NULL
		        for(j in 1:length(gr2.tmp[[1]])){
		           acceptorGene[j] <- getSeq(Hsapiens, seqnames(gr2.tmp[[1]][j]), start=start(ranges(gr2.tmp[[1]][j])), end=end(ranges(gr2.tmp[[1]][j])), as.character=TRUE)
		        }
		        acceptorGeneList[[z]] <- paste(acceptorGene, collapse="")
		        z <- z+1
	       }
	       names(acceptorGeneList) <- names(grl2)
	       #building fusions
	       fusions <- list()
	       names.fusions <- NULL
	       z <- 1
	       for(i in 1:length(donorGeneList)){
	        	 for(j in 1:length(acceptorGeneList)){
		              fusions[z] <- paste(donorGeneList[[i]], acceptorGeneList[[j]], sep="")
		              names.fusions[z] <- paste(names(donorGeneList)[i], names(acceptorGeneList)[j], sep=":")
		              z <- z + 1
		         }      
	       }
	       names(fusions) <- names.fusions
	       fusions <- DNAStringSet(unlist(fusions))
	       return(fusions)
         }else if(type=="gene"){
	           cat("\nNot implemented, yet\n")
	           return()
	   }
   }
}


