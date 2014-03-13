fusionPeptides <- function(chimeraSeq.output, annotation="hsUCSC"){
    if(annotation != "hsUCSC"){
    	cat("\nAnnotation not implemented, yet")
		return()
    }
	ann <- as.list(org.Hs.egUCSCKG)
	symbol <- as.list(org.Hs.egSYMBOL)
	f.name <- names(chimeraSeq.output)
	f.name <- strsplit(f.name, ":|\\.")
	f.name <- c(f.name[[1]][1],f.name[[1]][3])
	id1 <- lapply(ann, function(x,y){
		 tmp <- grep(y,x)
		 if(length(tmp)>0)
			return(x[tmp])
		 else
			 return(NA)
	 },y=f.name[1])
	id1.sym <- as.character(symbol[which(names(symbol)==names(id1[!is.na(id1)]))])
	if(length(id1.sym)==0)
		id1.sym <- NA 
	id1 <- as.character(id1[!is.na(id1)])
	id2 <- lapply(ann, function(x,y){
		 tmp <- grep(y,x)
		 if(length(tmp)>0)
			return(x[tmp])
		 else
			 return(NA)
	 },y=f.name[2])
	id2.sym <- as.character(symbol[which(names(symbol)==names(id2[!is.na(id2)]))])
	if(length(id2.sym)==0)
		id2.sym <- NA 
	id2 <- as.character(id2[!is.na(id2)])
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
	cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
	trs <- c(id1,id2)
	trs.sym <- c(id1.sym, id2.sym) 
	if(length(which(trs %in% names(cds_by_tx)))==1){
		cat("\nOnly ",paste(trs.sym[which(trs %in% names(cds_by_tx))],trs[which(trs %in% names(cds_by_tx))],sep="/")," is coding\n")
		return()
	}else if(length(which(trs %in% names(cds_by_tx)))==0){
		cat("\nNone of the two genes is coding\n")
		return()
	}
	#cds_seqs <- extractTranscriptsFromGenome(BSgenome.Hsapiens.UCSC.hg19, cds_by_tx[which(names(cds_by_tx) %in% trs)]) 
	cds_seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19, cds_by_tx[which(names(cds_by_tx) %in% trs)])
	proteome <- translate(cds_seqs) 
	names(proteome) <- names(cds_seqs) 
	p1 <- proteome[which(names(proteome)==id1)]
	p2 <- proteome[which(names(proteome)==id2)]
	fusion <- chimeraSeq.output[[1]]
	fusion.p <- list()
	fusion.p[[1]] <- translate(fusion[1:length(fusion)])
	fusion.p[[2]] <- translate(fusion[2:length(fusion)])
	fusion.p[[3]] <- translate(fusion[3:length(fusion)])
	fusion.p <- AAStringSet(fusion.p)
	alignment.p1 <- lapply(fusion.p, function(x,y){pairwiseAlignment(pattern=y, subject=x, type="local")},y=p1)
	scores.p1 <- sapply(alignment.p1, function(x)score(x))
	alignment.p1 <- alignment.p1[which(scores.p1==max(scores.p1))]
	alignment.p2 <- lapply(fusion.p, function(x,y){pairwiseAlignment(pattern=y, subject=x, type="local")},y=p2)
	scores.p2 <- sapply(alignment.p2, function(x)score(x))
	alignment.p2 <- alignment.p2[which(scores.p2==max(scores.p2))]
    if(length(which(scores.p1==max(scores.p1)))==1 & length(which(scores.p2==max(scores.p2)))==1){
	    if(which(scores.p1==max(scores.p1)) == which(scores.p2==max(scores.p2))){
		   pattern1 <- as.character(pattern(alignment.p1[[1]]))
		   fp1 <- matchPattern(pattern=pattern1, subject=fusion.p[[which(scores.p1==max(scores.p1))]])
		   pattern2 <- as.character(pattern(alignment.p2[[1]]))
		   fp2 <- matchPattern(pattern=pattern2, subject=fusion.p[[which(scores.p2==max(scores.p2))]])
		   if(length(as.character(fp1))==1 & length(as.character(fp2))==1){
		       cat("\nfused proteins are in frame\n")
    		   fusionp <- AAStringSet(c(paste(as.character(fp1), as.character(fp2), sep=""),as.character(fp1),as.character(fp2),as.character(p1),as.character(p2)))
	           names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")
		   }else if(length(as.character(fp1))==0 & length(as.character(fp2))==1){
		   	  cat("\n",id1.sym," does not provide any good match pattern with fusion protein frame: ",which(scores.p1==max(scores.p1)),"\n")
   		      fusionp <- AAStringSet(c("","",as.character(fp2),as.character(p1),as.character(p2)))
              names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")			  
			  #if this happen it is possible some issue in the fusion transcript reconstruction
		   }else if(length(as.character(fp1))==1 & length(as.character(fp2))==0){
		   	  cat("\n",id2.sym," does not provide any good match pattern with fusion protein frame: ",which(scores.p2==max(scores.p2)),"\n")
			  #if this happen it is possible some issue in the fusion transcript reconstruction
   		      fusionp <- AAStringSet(c("",as.character(fp1),"",as.character(p1),as.character(p2)))
              names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")
		   }
		}else{
 		   pattern1 <- as.character(pattern(alignment.p1[[1]]))
 		   fp1 <- matchPattern(pattern=pattern1, subject=fusion.p[[which(scores.p1==max(scores.p1))]])
 		   pattern2 <- as.character(pattern(alignment.p2[[1]]))
 		   fp2 <- matchPattern(pattern=pattern2, subject=fusion.p[[which(scores.p2==max(scores.p2))]])
		   if(length(as.character(fp1))==1 & length(as.character(fp2))==1){
			  cat("\nfused proteins are not in frame\n")
		      fusionp <- AAStringSet(c("",as.character(fp1),as.character(fp2),as.character(p1),as.character(p2)))
		      names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")
		   }else if(length(as.character(fp1))==0 & length(as.character(fp2))==1){
		   	  cat("\n",id1.sym," does not provide any good match pattern with fusion protein frame: ",which(scores.p1==max(scores.p1)),"\n")
   		      fusionp <- AAStringSet(c("","",as.character(fp2),as.character(p1),as.character(p2)))
              names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")			  
			  #if this happen it is possible some issue in the fusion transcript reconstruction
		   }else if(length(as.character(fp1))==1 & length(as.character(fp2))==0){
		   	  cat("\n",id2.sym," does not provide any good match pattern with fusion protein frame: ",which(scores.p2==max(scores.p2)),"\n")
   		      fusionp <- AAStringSet(c("",as.character(fp1),"",as.character(p1),as.character(p2)))
              names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")
			  #if this happen it is possible some issue in the fusion transcript reconstruction
		   }		
		}
	}else if(length(which(scores.p1==max(scores.p1)))==1){
 	 	   cat("\nfused proteins are not in frame\np2 is not providing a unique good alignment to the fusion\nIt is possible that fusion involves non-coding UTR.\n")
  		   pattern1 <- as.character(pattern(alignment.p1[[1]]))
  		   fp1 <- matchPattern(pattern=pattern1, subject=fusion.p[[which(scores.p1==max(scores.p1))]])
 		   fusionp <- AAStringSet(c("",as.character(fp1),"",as.character(p1),as.character(p2)))
 		   names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")		
	}else if(length(which(scores.p2==max(scores.p2)))==1){
 	 	   cat("\nfused proteins are not in frame.\np1 is not providing a unique good alignment to the fusion.\nIt is possible that fusion involves non-coding UTR.\n")
 		   pattern2 <- as.character(pattern(alignment.p2[[1]]))
 		   fp2 <- matchPattern(pattern=pattern2, subject=fusion.p[[which(scores.p2==max(scores.p2))]])
		   fusionp <- AAStringSet(c("","",as.character(fp2),as.character(p1),as.character(p2)))
		   names(fusionp) <- c("fusion","p1pep","p2pep","p1","p2")
		   output <- list(fusion.peptides=fusionp, fusion.transcript=chimeraSeq.output)					
	} 
	return(fusionp)
}
