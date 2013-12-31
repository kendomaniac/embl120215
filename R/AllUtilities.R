
###
filterList <- function(x,type=c("supporting.reads","fusion.names", "intronic", "annotated.genes", "read.through"),query, parallel=FALSE){
       if(type=="fusion.names"){
	       if(!is.character(query)){
		        cat("\nfiltering by fusion names needs to pass to the method vector of character names\n")
		        return()
		   }
		   if(parallel){
			     require(BiocParallel) || stop("\nMission BiocParallel library\n")
		         p <- MulticoreParam()
		         tmp <- fusionName(x, parallel=T)
			     loc <- which(tmp %in% query)
		   }else{
			     tmp <- fusionName(x)
			     loc <- which(tmp %in% query)
		   }
       }
       if(type=="supporting.reads"){
	       if(!is.numeric(query)){
		        cat("\nfiltering by supporting reads needs to pass to the method a numerical reads threshold\n")
		        return()
		   }
		   if(parallel){
			     require(BiocParallel) || stop("\nMission BiocParallel library\n")
		         p <- MulticoreParam()
		         tmp <- supportingReads(x, fusion.reads="spanning", parallel=T)
			       loc <- which(tmp >= query)
		   }else{
		       tmp <- supportingReads(x, fusion.reads="spanning")
		       loc <- which(tmp >= query)
		   }		
       }
       if(type=="intronic"){
	      cat("\n")
	        if(parallel){
			     require(BiocParallel) || stop("\nMission BiocParallel library\n")
		         p <- MulticoreParam()
		         tmp <- bplapply(x,.detectIntronic, BPPARAM=p)
		         tmp <- as.numeric(unlist(tmp))
                 loc <- which(tmp == 0)
        
		    }else{
                 tmp <- sapply(x,.detectIntronic)
                 loc <- which(tmp == 0)
            }
       }
       if(type=="annotated.genes"){
	        if(parallel){
			     require(BiocParallel) || stop("\nMission BiocParallel library\n")
		         p <- MulticoreParam()
		         tmp <- fusionName(x, parallel=T)
			     tmp.rm <- grep("chr", tmp)
			     loc <- setdiff(seq(1, length(tmp)),tmp.rm)
		    }else{
			     tmp <- fusionName(x)
			     tmp.rm <- grep("chr", tmp)
			     loc <- setdiff(seq(1, length(tmp)),tmp.rm)
		    }
		    
       }
       if(type=="read.through"){
	        if(parallel){
			     require(BiocParallel) || stop("\nMission BiocParallel library\n")
		         p <- MulticoreParam()
		         tmp <- fusionName(x, parallel=T)
			     tmp.rm <- grep("chr", tmp)
			     nsr <- bplapply(tmp.rm, function(x){
			              if(x[1]==x[2]){
				              return(1)
			              }else{return(0)}
			     }, BPPARAM=p)
			     nsr <- as.numeric(unlist(nsr))
			     loc <- which(nsr == 0)
		    }else{
			     tmp <- fusionName(x)
			     tmp.rm <- strsplit(tmp, ":")
			     nsr <- sapply(tmp.rm, function(x){
			              if(x[1]==x[2]){
				              return(1)
			              }else{return(0)}
			     })
			     nsr <- as.numeric(unlist(nsr))
			     loc <- which(nsr == 0)
		    }
		    
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





supportingReads <- function(list, fusion.reads=c("all","spanning"), parallel=FALSE){
	tmp <- list
	if(parallel){
	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
         p <- MulticoreParam()
         if(fusion.reads=="all"){
                nsr <- bplapply(tmp, function(x) x@fusionInfo$RescuedCount, BPPARAM=p)
				length.nsr <- sapply(nsr, length)
                nsr[which(length.nsr == 0)] <- 0
                nsr <- as.numeric(unlist(nsr))
         }else if(fusion.reads=="spanning"){
	            nsr <- bplapply(tmp, function(x) x@fusionInfo$SeedCount, BPPARAM=p)
				length.nsr <- sapply(nsr, length)
                nsr[which(length.nsr == 0)] <- 0
	            nsr <- as.numeric(unlist(nsr))
	     }
    }else{
	    if(fusion.reads=="all"){
                nsr <- sapply(tmp, function(x) x@fusionInfo$RescuedCount)
				length.nsr <- sapply(nsr, length)
                nsr[which(length.nsr == 0)] <- 0
				nsr <- as.numeric(unlist(nsr))
	    }else if(fusion.reads=="spanning"){ 
			    nsr <- sapply(tmp, function(x) x@fusionInfo$SeedCount)
				length.nsr <- sapply(nsr, length)
                nsr[which(length.nsr == 0)] <- 0
				nsr <- as.numeric(unlist(nsr))
	    }
    }
	return(nsr)	
}
###
fusionName <- function(list, parallel=FALSE){ 
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
starReads <- function(fusion.report, parallel=FALSE){
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
###
picardInstallation <- function(){
	    mydir <- getwd()
		picardDirLocation  <- paste(path.package("chimera", quiet = FALSE), "/picard", sep="")
		tmp.info <- dir(paste(path.package("chimera", quiet = FALSE)))
		if(length(grep("picard", tmp.info))>0){
		   unlink(picardDirLocation, recursive = T, force = T)
	    }
		dir.create(picardDirLocation, showWarnings = TRUE, recursive = FALSE)
		setwd(picardDirLocation)
		cat("\nBegin downloads of picard.....\n")
		download.file("http://sourceforge.net/projects/picard/files/picard-tools/1.92/picard-tools-1.92.zip/download", "picard.zip", mode="wb")
        cat("\nThe picard downloaded version is picard-tools-1.92\n")
        system(paste("unzip picard.zip", sep=""))
        system("cp -fR ./picard-tools-1.92/* .")
        system("rm -fR ./picard-tools-1.92/")
        system("chmod +x *")
        setwd(mydir)
        return()
}
###
validateSamFile <- function(input, output, mode=c("VERBOSE", "SUMMARY"), max.output="100"){
	tmp.info <- dir(paste(path.package("chimera", quiet = FALSE)))
	if(length(grep("picard", tmp.info))==0){
	   cat("\nIt seems that picard tools are not installed on chimera package path. \nPlease use picardInstallation() for picard tools installation\n")
       return()
    }
    picardDirLocation  <- paste(path.package("chimera", quiet = FALSE), "/picard", sep="")
	system(paste("java -jar ",picardDirLocation,"/ValidateSamFile.jar IGNORE_WARNINGS=true MAX_OUTPUT=",max.output," INPUT=",input," OUTPUT=",output,sep=""), wait=F)
    cat("SAM validation is running in background.")
	return(output)
}
##
removingErrorLine <- function(line.number, filenameIn, filenameOut){
	myline <- paste(" \'",line.number,"d\' " ,sep="")
	system(paste("sed ", myline, filenameIn," > ", filenameOut, sep=""), wait=F)
	cat("\nError line removal is running in background.\n")
}
##
filterSamReads <- function(input, output, filter=c("includeAligned","excludeAligned")){
	tmp.info <- dir(paste(path.package("chimera", quiet = FALSE)))
	if(length(grep("picard", tmp.info))==0){
	   cat("\nIt seems that picard tools are not installed on chimera package path. \nPlease use picardInstallation() for picard tools installation\n")
       return()
    }
    picardDirLocation  <- paste(path.package("chimera", quiet = FALSE), "/picard", sep="")
    tmp <- paste("tmp",gsub("[' '| :]","-", date()),sep="_")
#	system(paste("java -jar ",picardDirLocation,"/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=queryname INPUT=",input," OUTPUT=",paste(tmp,".sam",sep=""),sep=""), wait=T)
#	system(paste("java -jar ",picardDirLocation,"/FilterSamReads.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=unsorted FILTER=",filter," INPUT=",paste(tmp,".sam",sep="")," OUTPUT=",output,sep=""), wait=T)
    system(paste("java -jar ",picardDirLocation,"/FilterSamReads.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=unsorted FILTER=",filter," INPUT=",input," OUTPUT=",output,sep=""), wait=F)
 
   cat("SAM filtering is running in background.")
	return(output)
}
##
prettyPrint <- function(list, filename, fusion.reads=c("all","spanning")){
	all <- supportingReads(list, fusion.reads, parallel=FALSE)
    fusions.des <- sapply(list, function(x){
        chr.1 <- as.character(seqnames(fusionGRL(x)$gene1))
        chr.2 <- as.character(seqnames(fusionGRL(x)$gene2))
        end.1 <- end(fusionGRL(x)$gene1)
        strand.1 <- strand(fusionGRL(x)$gene1)
        start.2 <- start(fusionGRL(x)$gene2)
        strand.2 <- strand(fusionGRL(x)$gene2)
        gene1 <- elementMetadata(fusionGRL(x)$gene1)$KnownGene
        gene2 <- elementMetadata(fusionGRL(x)$gene2)$KnownGene
        trs1 <- elementMetadata(fusionGRL(x)$gene1)$KnownTranscript
        trs2 <- elementMetadata(fusionGRL(x)$gene2)$KnownTranscript
        junction1 <- elementMetadata(fusionGRL(x)$gene1)$FusionJunctionSequence
        junction2 <- elementMetadata(fusionGRL(x)$gene2)$FusionJunctionSequence
        junction <- paste(junction1, junction2, sep="")
        tmp <- paste(gene1, chr.1, end.1, strand.1, trs1, gene2, chr.2, start.2, strand.2, trs2, junction, sep="|")
    })
    fusions.des <- strsplit(fusions.des, "\\|")
    fusions.des <- as.data.frame(fusions.des)
    dimnames(fusions.des)[[1]] <-   c("gene1","chr.gene1","breakpoint.gene1", "strand.gene1","transcripts.gene1","gene2","chr.gene2","breakpoint.gene2","strand.gene2","transcripts.gene2","fusion.breakpoint")
    fusions.des <- t(fusions.des)
    fusions.des <- cbind(fusions.des, all)
    dimnames(fusions.des)[[2]] <- c("gene1","chr.gene1","breakpoint.gene1", "strand.gene1","transcripts.gene1","gene2","chr.gene2","breakpoint.gene2","strand.gene2","transcripts.gene2","fusion.breakpoint","supporting.reads")
    write.table(fusions.des, filename, sep="\t", row.names=F)
    return(paste("Fusion information is saved in ", filename, sep=""))
}

##
chimeraSeqSet <- function(fusions, parallel=F){
   if(parallel){
	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
         p <- MulticoreParam()
         trs <- bplapply(fusions, function(x){
                tmp <- chimeraSeqs(x, type="transcripts")
   			 tmp.names <- names(tmp)
   			 tmp <- tmp[[1]]
   			 return(list(seq=tmp, names=tmp.names))
             }, 
		 BPPARAM=p)
         trs.names <- sapply(trs, function(x)x[2])
         trs <- sapply(trs, function(x)x[1])
         trs <-DNAStringSet(trs)
         names(trs) <- as.character(trs.names)
   }else{
      trs <- lapply(fusions, function(x){
             tmp <- chimeraSeqs(x, type="transcripts")
			 tmp.names <- names(tmp)
			 tmp <- tmp[[1]]
			 return(list(seq=tmp, names=tmp.names))
          }
      )
      trs.names <- sapply(trs, function(x)x[2])
      trs <- sapply(trs, function(x)x[1])
      trs <-DNAStringSet(trs)
      names(trs) <- as.character(trs.names)
   }
   return(trs)
}





