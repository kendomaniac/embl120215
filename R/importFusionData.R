#functions:
#importFusionData
#.fmImport imports fusions from FusionMap
#.ffImport imports fusions from FusionFinder
#.fhImport imports fusions from FusionHunter
#.msImport imports fusions from MapSplice
#.thfImport imports fusions from TopHat-fusion
#.thfPostImport imports fusions from TopHat-fusion-post
#.dfImport imports fusions from deFuse
#.starImport imports fusions from star
#.starFset creat fset from star data
#.chimeraSeqs generates the fusion nucleotide sequence
#.rsImport importsusions from subjunc, reportAllJunctions=TRUE, of Rsubread package output for import is file with extension .fusion.txt
#the ideal situationis having mapped with STAR all reads and mapping with subjunc only the unmapped, option in STAR --outReadsUnmapped Fastx
#tophatInstallation
#tophatRun
#plotCoverage
#filterList
#supportingReads
#fusionName

.geneLevelAnnotation <- function(genome=c("hg19","mm9","hg38","mm10")){
	if(genome=="hg19"){
		TxDb.Mmusculus.UCSC.mm9.knownGene <- NULL
		TxDb.Mmusculus.UCSC.mm10.knownGene <- NULL
		TxDb.Hsapiens.UCSC.hg38.knownGene <- NULL
		Mus.musculus <- NULL
		Homo.sapiens <- NULL
		BSgenome.Hsapiens.NCBI.GRCh38 <- NULL
		
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
		genes.gr <- genes(txdb)
		simbols.eg <- select(Homo.sapiens,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
		if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
			     elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
				 return(genes.gr)
		}else{
			cat("\nGENEID in Homo.sapiens are not aligned with GENEID in TxDb.Hsapiens.UCSC.hg19.knownGene\n")
			return()
		}
	}else if(genome=="mm9"){
		require(TxDb.Mmusculus.UCSC.mm9.knownGene)||stop("\nMissing TxDb.Mmusculus.UCSC.mm9.knownGene library\n")
		txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
		genes.gr <- genes(txdb)
		require(Mus.musculus)||stop("\nMissing Mus.musculus library\n")
		simbols.eg <- select(Mus.musculus,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
		if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
			     elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
				 return(genes.gr)
		}else{
			cat("\nGENEID in Mus.musculus are not aligned with GENEID in TxDb.Mmusculus.UCSC.mm9.knownGene\n")
			return()
		}
	}else if(genome=="mm10"){
		require(TxDb.Mmusculus.UCSC.mm10.knownGene)||stop("\nMissing TxDb.Mmusculus.UCSC.mm10.knownGene library\n")
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
		genes.gr <- genes(txdb)
		require(Mus.musculus)||stop("\nMissing Mus.musculus library\n")
		simbols.eg <- select(Mus.musculus,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
		if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
			     elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
				 return(genes.gr)
		}else{
			cat("\nGENEID in Mus.musculus are not aligned with GENEID in TxDb.Mmusculus.UCSC.mm10.knownGene\n")
			return()
		}
	}else if(genome=="hg38"){
		require(TxDb.Hsapiens.UCSC.hg38.knownGene)||stop("\nMissing TxDb.Hsapiens.UCSC.hg38.knownGene library\n")
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
		genes.gr <- genes(txdb)
		simbols.eg <- select(Homo.sapiens,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
		if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
			     elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
				 return(genes.gr)
		}else{
			cat("\nGENEID in Mus.musculus are not aligned with GENEID in TxDb.Hsapiens.UCSC.hg38.knownGene\n")
			return()
		}
	} 
}

importFusionData <- function(format, filename, ...)
{
    switch(format,
           "bellerophontes" = .bfImport(filename),
           "defuse" = .dfImport(filename),
           "fusionfinder" = .ffImport(filename),
           "fusionhunter" = .fhImport(filename),
	       "mapsplice" = .msImport(filename, ...),
	       "tophat-fusion" = .thfImport(filename, ...),
		   "tophat-fusion-post" = .thfImport(filename, ...),
	       "chimerascan" = .csImport(filename, ...),
           "fusionmap" = .fmImport(filename, ...),
           "star" = .starImport(filename, ...),
		   "rsubread" = .rsImport(filename, ...),
		   "fusioncatcher" = .fcImport(filename, ...)
    )
}
#functions
#import from subjunc of Rsubread package
#tmp <- importFusionData("rsubread","spike4_subread.bam.fusion.txt", org="hs", min.support=10)

.rsImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10"), min.distance=700000, min.support=10, parallel=FALSE){
	     Mmusculus <- NULL
         if(parallel){ 
           require(BiocParallel) || stop("\nMission BiocParallel library\n")
           p <- MulticoreParam()
        }
	    report <- read.table(fusion.report, sep="\t", header=T, comment.char="")
		good1 <- report[which(as.character(report$X.Chr)!=as.character(report$Chr)),]
		good2 <- report[which(as.character(report$X.Chr)==as.character(report$Chr)),]
		good2 <- good2[which(abs(good2$Location - good2$Location.1)>= min.distance),]
		report <- rbind(good1,good2)
		
		#removing chrM
	    if(length(which(as.character(report$X.Chr)=="chrM"))>0){
	         cat("\nchrM is removed from fusion acceptor\n")
			 report <- report[setdiff(seq(1:dim(report)[1]),which(as.character(report$X.Chr)=="chrM")),]
	    }
	    if(length(which(as.character(report$Chr)=="chrM"))>0){
	         cat("\nchrM is removed from fusion donor\n")
			 report <- report[setdiff(seq(1:dim(report)[1]),which(as.character(report$Chr)=="chrM")),]
	    }
	    #removing non canonical chrs1
	    chr.g1.l <- sapply(as.character(report$X.Chr), nchar)
	    if(length(which(as.numeric(chr.g1.l) > 5)) >0){
	             removed <- as.character(unique(report$X.Chr[which(as.numeric(chr.g1.l) > 5)]))
	             cat("\nThe following chrs were removed from fusion acceptor:\n",removed,"\n")
				 report <- report[which(as.numeric(chr.g1.l) <= 5),]
	    }
	    
	    chr.g2.l <- sapply(as.character(report$Chr), nchar)
	    if(length(which(as.numeric(chr.g2.l) > 5)) >0){
	             removed <- as.character(unique(report$Chr[which(as.numeric(chr.g2.l) > 5)]))
	             cat("\nThe following chrs were removed from fusion donor:\n",removed,"\n")
				 report <- report[which(as.numeric(chr.g2.l) <= 5),]
	    }
	    #removing non canonical chrs2
#	    if(length(setdiff(seq(1, dim(report)[1]),grep("chr",as.character(report$X.Chr))))>0){
#			 removed <- as.character(unique(report$X.Chr[setdiff(seq(1, dim(report)[1]),grep("chr",as.character(report$X.Chr))]))
#	         cat("\nThe following chrs were removed from fusion donor:\n",removed,"\n")
#			 report <- report[grep("chr",as.character(report$X.Chr), ]
#	    }
#	    if(length(setdiff(seq(1, dim(report)[1]),grep("chr",as.character(report$Chr))))>0){
#		      removed <- as.character(unique(report$Chr[setdiff(seq(1, dim(report)[1]),grep("chr",as.character(report$Chr))]))
#              cat("\nThe following chrs were removed from fusion acceptor:\n",removed,"\n")
#			  report <- report[grep("chr",as.character(report$Chr), ]
#	    }
	    		
		report <- report[which(report$nSupport >= min.support),]
		cat(paste("\n",dim(report)[1]," detected fusions\n",sep=""))
		report <- apply(report, 1, function(x){
			tmp <- list(c(as.character(unlist(x[1])), as.numeric(x[2]), as.character(unlist(x[3])), as.numeric(x[4]), as.character(unlist(x[5])), as.numeric(x[6])))
		})
#		fusionreads.loc <- new("GAlignments")
	    #loading annotation
		grHs <- .geneLevelAnnotation(genome=org)
		
		#creating object
	 #   fusionreads.loc <- fusionreads.loc[[1]]
	    cat("\n") 
	    .fusionInfo <- function(x,y){
			cat(".")
			x <- x[[1]]
            grG1 <-  GRanges(seqnames = as.character(unlist(x[1])), ranges = IRanges(start = (as.numeric(x[2]) - 30), end= as.numeric(x[2])))
            grG2 <-  GRanges(seqnames = as.character(unlist(x[3])), ranges = IRanges(start = as.numeric(x[4]), end= (as.numeric(x[4])+30)))
            
			tmpG1 <- findOverlaps(grG1, y, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG1)){
                g1 <- elementMetadata(grHs[tmpG1])$symbol            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            
			tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG2)){
                g2 <- elementMetadata(grHs[tmpG2])$symbol            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            if(org=="hg19"){
   		     fs.1 <- as.character(getSeq(Hsapiens, grG1))
   		     fs.2 <- as.character(getSeq(Hsapiens, grG2))		
            }else if(org=="mm9"){
   	         require(BSgenome.Mmusculus.UCSC.mm9) || stop("\nMissing BSgenome.Mmusculus.UCSC.mm9 library\n")
   		     fs.1 <- as.character(getSeq(Mmusculus, grG1))
   		     fs.2 <- as.character(getSeq(Mmusculus, grG2))		
            }else if(org=="mm10"){
   	         require(BSgenome.Mmusculus.UCSC.mm10) || stop("\nMissing BSgenome.Mmusculus.UCSC.mm9 library\n")
   		     fs.1 <- as.character(getSeq(Mmusculus, grG1))
   		     fs.2 <- as.character(getSeq(Mmusculus, grG2))		
            }else if(org=="hg38"){
	 		    ################################################################################
	 		    #da sistemare nelle funzioni giuste per il momento mi serve solo come remind
	 		    require(BSgenome.Hsapiens.NCBI.GRCh38)||stop("\nMissing BSgenome.Hsapiens.NCBI.GRCh38 library\n")
	 		    Hsapiens <- BSgenome.Hsapiens.NCBI.GRCh38
	 		    seqlevelsStyle(Hsapiens) <- "UCSC" ## convert to UCSC style
	 		    ################################################################################
   		        fs.1 <- as.character(getSeq(Hsapiens, grG1))
   		        fs.2 <- as.character(getSeq(Hsapiens, grG2))		
            }
            
   	        gr1 <- GRanges(seqnames = as.character(unlist(x[1])),
   		            ranges = IRanges(start = (as.numeric(x[2]) - 30), end= as.numeric(x[2])),
   		            KnownGene = as.character(g1),
   		            KnownTranscript =  "",
   		            KnownExonNumber = "",
   		            KnownTranscriptStrand = "",
   		            FusionJunctionSequence =  fs.1)
   	        gr2 <- GRanges(seqnames = as.character(unlist(x[3])),
   				   ranges = IRanges(start = as.numeric(x[4]), end= (as.numeric(x[4])+30)),
   				   KnownGene = as.character(g2),
   				   KnownTranscript =  "",
   				   KnownExonNumber = "",
   				   KnownTranscriptStrand = "",
   				   FusionJunctionSequence =  fs.2)
   		    grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
   		    fusionData <- new("list", fusionTool="FusionMap", 
   		                                 UniqueCuttingPositionCount="", 
   		                                 SeedCount=x[6], 
   		                                 RescuedCount=x[6], 
   		                                 SplicePattern="",
   		                                 FusionGene=paste(g1,g2, sep="->"),
   		                                 frameShift=""
   		    )
   		    return(new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet")))		                          
	   }
	   if(parallel){           	
	       fusionList <- bplapply(report,.fusionInfo, grHs, BPPARAM=p)
	   }else{
		   fusionList <- lapply(report,.fusionInfo, grHs)
	   }
	   cat("\n") 
	   return(fusionList)
}


#FusionMap import
.fmImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	    report <- read.table(fusion.report, sep="\t", header=T)
#		fusionreads.loc <- new("GAlignments")
	    #loading annotation
        grHs <- .geneLevelAnnotation(genome=org)		#creating object
	 #  fusionreads.loc <- fusionreads.loc[[1]]            	
	    fusionList <- list()
	for(i in 1:dim(report)[1]){
	  if(as.character(report$FusionJunctionSequence[i])!=""){	
		 strand1 <- NULL
		 strand2 <- NULL
		 if(report$Strand[i] == as.character("++")) {strand1 <- "+"; strand2 <- "+"}
		 if(report$Strand[i] == as.character("+-")) {strand1 <- "+"; strand2 <- "-"}
		 if(report$Strand[i] == as.character("-+")) {strand1 <- "-"; strand2 <- "+"}
		 if(report$Strand[i] == as.character("--")) {strand1 <- "-"; strand2 <- "-"}

		 fs.tmp <- strsplit(as.character(report$FusionJunctionSequence[i]),"[a-z]")
		 fs.1 <- fs.tmp[[1]][1]
		 fs.tmp <- strsplit(as.character(report$FusionJunctionSequence[i]),"[A-Z]")
		 fs.2 <- fs.tmp[[1]][length(fs.tmp[[1]])]
		 
		 #detecting genes involved in fusions
         grG1 <-  GRanges(seqnames = paste("chr",report$Chromosome1[i],sep=""),
		            ranges = IRanges(start = (report$Position1[i] - nchar(fs.1)), end= report$Position1[i]),
		            strand = strand1)
         grG2 <-  GRanges(seqnames = paste("chr",report$Chromosome2[i],sep=""),
					ranges = IRanges(start = report$Position2[i], end= (report$Position2[i] + nchar(fs.2))),
				    strand = strand2)
					
         tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
         if(!is.na(tmpG1)){
             g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
         }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
		 
         tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
         if(!is.na(tmpG2)){
			 g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
         }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

         
		 if(length(as.character(report$KnownTranscript1[i])) == 0){
			tmpT1 <- NULL
		 } else{
		   	tmpT1 <- as.character(report$KnownTranscript1[i])	
		 } 
		 if(length(as.character(report$KnownTranscript2[i])) == 0){
			tmpT2 <- NULL
		 } else{
		   	tmpT2 <- as.character(report$KnownTranscript2[i])	
		 } 

		 if(length(as.character(report$KnownExonNumber1[i])) == 0){
			tmpEn1 <- NULL
		 } else{
		   	tmpEn1 <- as.character(report$KnownExonNumber1[i])	
		 } 
		 if(length(as.character(report$KnownExonNumber2[i])) == 0){
			tmpEn2 <- NULL
		 } else{
		   	tmpEn2 <- as.character(report$KnownExonNumber2[i])	
		 } 

		 if(length(as.character(report$KnownTranscriptStrand1[i])) == 0){
			tmpTs1 <- NULL
		 } else{
		   	tmpTs1 <- as.character(report$KnownTranscriptStrand1[i])	
		 } 
		 if(length(as.character(report$KnownTranscriptStrand2[i])) == 0){
			tmpTs2 <- NULL
		 } else{
		   	tmpTs2 <- as.character(report$KnownTranscriptStrand2[i])	
		 } 


	     gr1 <- GRanges(seqnames = paste("chr",report$Chromosome1[i],sep=""),
		            ranges = IRanges(start = (report$Position1[i] - nchar(fs.1)), end= report$Position1[i]),
		            strand = strand1,
		            KnownGene = as.character(g1),
		            KnownTranscript =  tmpT1,
		            KnownExonNumber = tmpEn1,
		            KnownTranscriptStrand = tmpTs1,
		            FusionJunctionSequence =  fs.1)
	     gr2 <- GRanges(seqnames = paste("chr",report$Chromosome2[i],sep=""),
				   ranges = IRanges(start = report$Position2[i], end= (report$Position2[i] + nchar(fs.2))),
				   strand = strand2,
				   KnownGene = as.character(g2),
				   KnownTranscript =  tmpT2,
				   KnownExonNumber = tmpEn2,
				   KnownTranscriptStrand = tmpTs2,
				   FusionJunctionSequence =  fs.2)
		grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
		fusionData <- new("list", fusionTool="FusionMap", 
		                                 UniqueCuttingPositionCount=report[i,2], 
		                                 SeedCount=report[i,3], 
		                                 RescuedCount=report[i,4], 
		                                 SplicePattern=as.character(report$SplicePattern[i]),
		                                 FusionGene=as.character(report$FusionGene[i]),
		                                 frameShift=as.character(report$FrameShift[i])
		)
		fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))		             
     }
	}
	   return(fusionList)
}

#FusionHunter import
.fhImport <- function(fusion.report){
	    con <- file(fusion.report, "r", blocking = FALSE)
		tmp <- readLines(con)
		close(con)
		tmp <- tmp[setdiff(seq(1, length(tmp)),grep("#",tmp))]
		for(i in 1:10){
			tmp <- gsub("  ", " ", tmp)
		#	nchar(tmp[1])
			i <- i + 1
		} 
		tmp <- gsub(" ", "\t", tmp)
		tmp <- tmp[setdiff(seq(1, length(tmp)),grep("---",tmp))]
		writeLines(tmp, paste(fusion.report,"tmp", sep="."), sep="\n")
	    report <- read.table(paste(fusion.report,"tmp", sep="."), sep="\t", header=F)
	    unlink(tmp, force =T)
#	    fusionreads.loc <- new("GAlignments") 
	    #strand           	
        strand1 <- "*"
        strand2 <- "*"
        #FusionGene
        fg <- paste(report[,7], report[,9], sep="->")
        fg.u <- unique(fg)
        #KnownGene1
        g.tmp <- strsplit(fg.u, "->")
        g1 <- sapply(g.tmp, function(x)x[1])
		#KnownGene2
        g2 <- sapply(g.tmp, function(x)x[2])
        #KnownTranscript1
        t1 <- NA
        #KnownTranscript2
        t2 <- NA
        #KnownExonNumber1
        e1 <- NA
        #KnownExonNumber2
        e2 <- NA
        #KnownTranscriptStrand1
        ts1 <- NA
        #KnownTranscriptStrand2
        ts2 <- NA
        #FusionJunctionSequence1-2
        #SeedCount1-2
        fs1 <- list()
        seed1 <- list()
        fs2 <- list()
        seed2 <- list()
        #fusion loc
        start1 <- list()
        end1 <- list()
        chr1 <- list()
        start2 <- list()
        end2 <- list()
        chr2 <- list()
        for(i in 1:length(fg.u)){
	        pos.tmp <- as.character(report[which(fg == fg.u[i]),5])
	        seed1[[i]] <- length(pos.tmp)
	        report.tmp <- report[which(fg == fg.u[i]),]
	        pos.tmp1 <- strsplit(pos.tmp, ":")
	        chr1.tmp <- sapply(pos.tmp1, function(x)x[1])
	        pos.tmp2 <- sapply(pos.tmp1, function(x)x[2])
	        pos.tmp3 <- strsplit(pos.tmp2, "-")
	        pos.start <- sapply(pos.tmp3, function(x)x[1])
	        pos.end <- sapply(pos.tmp3, function(x)x[2])
	        fs1.tmp <- as.character(report.tmp[which(pos.start == min(pos.start)),3])
	        fs1.tmp <- fs1.tmp[which(nchar(fs1.tmp) == max(nchar(fs1.tmp)))]
            start1.tmp <-  pos.start[which(nchar(fs1.tmp) == max(nchar(fs1.tmp)))]
            end1.tmp <- pos.end[which(nchar(fs1.tmp) == max(nchar(fs1.tmp)))]
            chr1.tmp <- chr1.tmp[which(nchar(fs1.tmp) == max(nchar(fs1.tmp)))]
	        if(length(fs1.tmp) > 1) {
		         fs1.tmp <- fs1.tmp[1]
		         start1.tmp <- start1.tmp[1]
		         end1.tmp <- end1.tmp[1]
		         chr1.tmp <- chr1.tmp[1]
		    }
	        fs1[[i]] <- fs1.tmp
            start1[[i]] <- start1.tmp
            end1[[i]] <- end1.tmp
            chr1[[i]] <- chr1.tmp
	        #
	        pos.tmp <- as.character(report[which(fg == fg.u[i]),6])
	        seed2[[i]] <- length(pos.tmp)
	        pos.tmp1 <- strsplit(pos.tmp, ":")
	        chr2.tmp <- sapply(pos.tmp2, function(x)x[1])
	        pos.tmp2 <- sapply(pos.tmp1, function(x)x[2])
	        pos.tmp3 <- strsplit(pos.tmp2, "-")
	        pos.start <- sapply(pos.tmp3, function(x)x[1])
	        pos.end <- sapply(pos.tmp3, function(x)x[2])
	        fs2.tmp <- as.character(report.tmp[which(pos.end == max(pos.end)),4])
	        fs2.tmp <- fs2.tmp[which(nchar(fs2.tmp) == max(nchar(fs2.tmp)))]
	        start2.tmp <-  pos.start[which(nchar(fs2.tmp) == max(nchar(fs2.tmp)))]
            end2.tmp <- pos.end[which(nchar(fs2.tmp) == max(nchar(fs2.tmp)))]
            chr2.tmp <- chr2.tmp[which(nchar(fs2.tmp) == max(nchar(fs2.tmp)))]
	        if(length(fs2.tmp) > 1) {
		         fs2.tmp <- fs2.tmp[1]
		         start2.tmp <- start2.tmp[1]
		         end2.tmp <- end2.tmp[1]
		         chr2.tmp <- chr2.tmp[1]
		    }
	        fs2[[i]] <- fs2.tmp
	        start2[[i]] <- start2.tmp
            end2[[i]] <- end2.tmp
            chr2[[i]] <- chr2.tmp
           
        }
        #UniqueCuttingPositionCount
        cut <- 0
        #RescuedCount
        rc <- 0
        #SplicePattern
        sp <- "-"
        #frameShift
        fs <- "-"
        #creating objects
        fusionList <- list()
        for(i in 1:length(fg.u)){
		     gr1 <- GRanges(seqnames = chr1[[i]],
			            ranges = IRanges(start = as.numeric(start1[[i]]), end= as.numeric(end1[[i]])),
			            strand = strand1,
			            KnownGene = g1[i],
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1[[i]])
		     gr2 <- GRanges(seqnames = chr2[[i]],
					   ranges = IRanges(start = as.numeric(start2[[i]]), end= as.numeric(end2[[i]])),
					   strand = strand2,
					   KnownGene = g2[i],
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2[[i]])
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="FusionHunter", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1[[i]], 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=fg.u[i],
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)
}
#tophat-fusion import
.thfImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	    report <- read.table(fusion.report, sep="\t", header=F)
	    if(dim(report)[2]!=23){
		    cat("\nThe data structure does not seem to be that of tophat-fusion\n")
		    return()
	    }
	    names(report) <- c("chrs.fusion","gene1.fusion.loc", "gene2.fusion.loc", "strand", "fusion.reads","encompassing.reads","spanning.reads",
		"contraddicting.reads","g1.nts.fusion", "g2.nts.fusion","unused1","separator1","unused2","separator2","seq1","separator3","seq2","separator4",
		"coverage1", "separator5","coverage2","separator6","unused2")
	#    fusionreads.loc <- new("GAlignments") 
        #KnownGene1
        g1 <- NA
		#KnownGene2
        g2 <- NA
        #KnownTranscript1
        t1 <- NA
        #KnownTranscript2
        t2 <- NA
        #KnownExonNumber1
        e1 <- NA
        #KnownExonNumber2
        e2 <- NA
        #KnownTranscriptStrand1
        ts1 <- NA
        #KnownTranscriptStrand2
        ts2 <- NA
        #FusionJunctionSequence1-2
        #SeedCount1-2
        #UniqueCuttingPositionCount
        cut <- 0
        #SplicePattern
        sp <- "-"
        #frameShift
        fs <- "-"
        #creating objects
        fusionList <- list()
        #loading annotation
        grHs <- .geneLevelAnnotation(genome=org)        
		for(i in 1:dim(report)[1]){
	     #   cat("\n",i)
	        strand1 <- NULL
		    strand2 <- NULL
		    if(report$strand[i] == as.character("ff")) {strand1 <- "+"; strand2 <- "+"}
		    if(report$strand[i] == as.character("fr")) {strand1 <- "+"; strand2 <- "-"}
		    if(report$strand[i] == as.character("rf")) {strand1 <- "-"; strand2 <- "+"}
		    if(report$strand[i] == as.character("rr")) {strand1 <- "-"; strand2 <- "-"}
            #spanning reads		 
	        seed1 <- report$spanning.reads[i]
	        #RescuedCount encompassing.reads
	        rc <- report$encompassing.reads[i]
	        pos.tmp <- strsplit(as.character(report$chrs.fusion[i]), "-")
	        #junction sequence1
	        fs1.tmp <- strsplit(as.character(report$seq1[i]), " ")
	        fs1 <- fs1.tmp[[1]][1]
	        #defining start1
	        coverage1.tmp <- strsplit(as.character(report$coverage1[i]), " ")
            start1 <- report$gene1.fusion.loc[i] - length(coverage1.tmp[[1]])
            end1 <- report$gene1.fusion.loc[i]
            chr1 <- pos.tmp[[1]][1]
	        #junction2
	        fs2.tmp <- strsplit(as.character(report$seq2[i]), " ")
		    fs2 <- fs2.tmp[[1]][1]
		    #defining start2
		    coverage2.tmp <- strsplit(as.character(report$coverage2[i]), " ")
	        start2 <- report$gene2.fusion.loc[i]
            end2 <- report$gene2.fusion.loc[i] + length(coverage2.tmp[[1]])
            chr2 <- pos.tmp[[1]][2]
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1), strand = strand1)
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2), strand = strand2)
            tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG1)){
                g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG2)){
                g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			         #   strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
				#	   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="TopHat-fusion", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)
}

#tophat-fusion import
.thfPostImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	    report <- read.table(fusion.report, sep="\t", header=F)
	    if(dim(report)[2]!=11){
		    cat("\nThe data structure does not seem to be that of tophat-fusion-post file\n")
		    return()
	    }
	    names(report) <- c("SampleName","gene1", "gene1.chr", "gene1.fusion.loc", "gene2","gene2.chr","gene2.fusion.loc",
		"spanning","encompassing", "encompassing.w.spanning","other")
	#    fusionreads.loc <- new("GAlignments") 
        #KnownGene1
        g1 <- NA
		#KnownGene2
        g2 <- NA
        #KnownTranscript1
        t1 <- NA
        #KnownTranscript2
        t2 <- NA
        #KnownExonNumber1
        e1 <- NA
        #KnownExonNumber2
        e2 <- NA
        #KnownTranscriptStrand1
        ts1 <- NA
        #KnownTranscriptStrand2
        ts2 <- NA
        #FusionJunctionSequence1-2
        #SeedCount1-2
        #UniqueCuttingPositionCount
        cut <- 0
        #SplicePattern
        sp <- "-"
        #frameShift
        fs <- "-"
        #creating objects
        fusionList <- list()
        #loading annotation
        grHs <- .geneLevelAnnotation(genome=org) 
        for(i in 1:dim(report)[1]){
	     #   cat("\n",i)
	        strand1 <- NULL
		    strand2 <- NULL
            #spanning reads		 
	        seed1 <- report$encompassing.w.spanning[i]
	        #RescuedCount encompassing.reads
	        rc <- report$encompassing[i]
	        #pos.tmp <- strsplit(as.character(report$chrs.fusion[i]), "-")
	        #junction sequence1
	        fs1 <- ""
	        #defining start1
	        #coverage1.tmp <- strsplit(as.character(report$coverage1[i]), " ")
            start1 <- report$gene1.fusion.loc[i] - 30
            end1 <- report$gene1.fusion.loc[i]
            chr1 <- as.character(report$gene1.chr[i])
	        #junction2
		    fs2 <- ""
		    #defining start2
		    #coverage2.tmp <- strsplit(as.character(report$coverage2[i]), " ")
	        start2 <- report$gene2.fusion.loc[i]
            end2 <- report$gene2.fusion.loc[i] + 30
            chr2 <- as.character(report$gene2.chr[i])
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1))
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2))
            tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG1)){
                g1 <- elementMetadata(grHs[tmpG1])$symbol            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG2)){
                g2 <- elementMetadata(grHs[tmpG2])$symbol            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			         #   strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
				#	   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="TopHat-fusion", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)
}


.msImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	    flankcase.description <- seq(0,6)
	    names(flankcase.description) <- c("others","ATAC","GTAT","CTGC","GCAG","GTAG","CTAC")
	    report <- read.table(fusion.report, sep="\t", header=F, skip=1)
	    names(report) <- c("chrom","donerEnd","acceptorStart","name","coverage","strand","itemRgb","blockCount","blockSizes","blockStarts","entropy","flank.case","flank.string","fusion.junction", "min.mismatch","max.mismatch","average.mismatch","maximal.of.minimal.doner.site.length","maximal.of.minimal.acceptor.site.length","minimal.anchor.difference")
#	    fusionreads.loc <- new("GAlignments") 
		#KnownGene1
	    g1 <- NA
	    #KnownGene2
		g2 <- NA
		#KnownTranscript1
		t1 <- NA
		#KnownTranscript2
		t2 <- NA
		#KnownExonNumber1
		e1 <- NA
		#KnownExonNumber2
		e2 <- NA
		#KnownTranscriptStrand1
		ts1 <- NA
		#KnownTranscriptStrand2
		ts2 <- NA
		#FusionJunctionSequence1-2
		#SeedCount1-2
		#UniqueCuttingPositionCount
		cut <- 0
		#RescuedCount
        rc <- 0
		#frameShift
		fs <- "-"
		#creating objects
		fusionList <- list()
		#loading annotation
        grHs <- .geneLevelAnnotation(genome=org)		
		for(i in 1:dim(report)[1]){
			#SplicePattern
			sp <- names(flankcase.description[which(flankcase.description == report$flank.case[i])])
            #strand
	        strand1 <- NULL
		    strand2 <- NULL
		    if(report$strand[i] == as.character("++")) {strand1 <- "+"; strand2 <- "+"}
		    if(report$strand[i] == as.character("+-")) {strand1 <- "+"; strand2 <- "-"}
		    if(report$strand[i] == as.character("-+")) {strand1 <- "-"; strand2 <- "+"}
		    if(report$strand[i] == as.character("--")) {strand1 <- "-"; strand2 <- "-"}
            #spanning reads		 
	        seed1 <- report$coverage[i]
	        fs1 <- substr(as.character(report$fusion.junction[i]), 1, 60)
	        #defining start1
            start1 <- report$donerEnd[i] - 60
            end1 <- report$donerEnd[i]
            chr.tmp <- strsplit(as.character(report$chrom[i]),"_")
            chr1 <- chr.tmp[[1]][1]
	        #junction2
		    fs2 <- substr(as.character(report$fusion.junction[i]), 61, 120)
		    #defining start2
	        start2 <- report$acceptorStart[i]
            end2 <- report$acceptorStart[i] + 60
            chr2 <- chr.tmp[[1]][2]
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1), strand = strand1)
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2), strand = strand2)
            tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG1)){
                g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG2)){
                g2 <- elementMetadata(grHs[tmpG2])$symbol           	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			            strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
					   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="mapSplice", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)
		
	
}	

.dfImport <- function(fusion.report){
	    report <- read.table(fusion.report, sep="\t", header=T)
#	    fusionreads.loc <- new("GAlignments") 
		#KnownGene1
	    g1 <- NA
	    #KnownGene2
		g2 <- NA
		#KnownTranscript1
		t1 <- NA
		#KnownTranscript2
		t2 <- NA
		#KnownExonNumber1
		e1 <- NA
		#KnownExonNumber2
		e2 <- NA
		#KnownTranscriptStrand1
		ts1 <- NA
		#KnownTranscriptStrand2
		ts2 <- NA
		#FusionJunctionSequence1-2
		#SeedCount1-2
		#UniqueCuttingPositionCount
		cut <- 0
		#RescuedCount
        rc <- 0
		#frameShift
		fs <- "-"
		#SplicePattern
		sp <- "-"
		#creating objects
		fusionList <- list()
		for(i in 1:dim(report)[1]){
            #strand
	        strand1 <- as.character(report$gene_strand1[i])
		    strand2 <- as.character(report$gene_strand2[i])
            #spanning reads		 
	        seed1 <- report$splitr_count[i]
	        fs.tmp <- strsplit(as.character(report$splitr_sequence[i]),'\\|')
	        fs1 <- fs.tmp[[1]][1]
	        #defining start1
            start1 <- report$genomic_break_pos1[i] - nchar(fs1)
            end1 <- report$genomic_break_pos1[i]
            chr.tmp <- strsplit(as.character(report$chrom[i]),"_")
            chr1 <- paste("chr",as.character(report$gene_chromosome1[i]), sep="")
	        #junction2
		    fs2 <- fs.tmp[[1]][2]
		    #defining start2
	        start2 <- report$genomic_break_pos2[i]
            end2 <- report$genomic_break_pos2[i] + nchar(fs2)
            chr2 <-  paste("chr",as.character(report$gene_chromosome2[i]), sep="")
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1), strand = strand1)
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2), strand = strand2)
            if(!is.na(report$gene_name1[i])){
                g1 <- as.character(report$gene_name1[i])            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            if(!is.na(report$gene_name2[i])){
                g2 <- as.character(report$gene_name2[i]) 	            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			            strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
					   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="deFuse", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)
		
	
}	

.ffImport <- function(fusion.report){
	    report <- read.table(fusion.report, sep="\t", header=T)
#	    fusionreads.loc <- new("GAlignments") 
		#KnownGene1
	    g1 <- NA
	    #KnownGene2
		g2 <- NA
		#KnownTranscript1
		t1 <- NA
		#KnownTranscript2
		t2 <- NA
		#KnownExonNumber1
		e1 <- NA
		#KnownExonNumber2
		e2 <- NA
		#KnownTranscriptStrand1
		ts1 <- NA
		#KnownTranscriptStrand2
		ts2 <- NA
		#FusionJunctionSequence1-2
		#SeedCount1-2
		#UniqueCuttingPositionCount
		cut <- 0
		#RescuedCount
        rc <- 0
		#frameShift
		fs <- "-"
		#SplicePattern
		sp <- "-"
		#SpliceJunction1
		fs1 <- "-"
		#SpliceJunction1
		fs2 <- "-"
		#creating objects
		fusionList <- list()
		for(i in 1:dim(report)[1]){
            #strand
	        if(report$G1_str[i] == -1) strand1 <- "-" else strand1 <- "+"
		    if(report$G2_str[i] == -1) strand2 <- "-" else strand2 <- "+"
            #spanning reads		 
	        seed1 <- report$totalreads[i]
	        #defining start1
	        loc1.tmp <- strsplit(as.character(report$G1_block[i]),"-")
            start1 <- as.numeric(loc1.tmp[[1]][1])
            end1 <- as.numeric(loc1.tmp[[1]][2])
            chr1 <- sub("chromosome_","chr",as.character(report$G1_chromosome[i]))
		    #defining start2
		    loc2.tmp <- strsplit(as.character(report$G2_block[i]),"-")
	        start2 <- as.numeric(loc2.tmp[[1]][1])
            end2 <- as.numeric(loc2.tmp[[1]][2])
            chr2 <-  paste("chr",as.character(report$gene_chromosome2[i]), sep="")
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1), strand = strand1)
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2), strand = strand2)
            tmp.G1 <- strsplit(as.character(report$G1_Ensembl_HGNC_ID[i]), "\\(")
            G1 <- gsub(" ","",tmp.G1[[1]][1])
            tmpG1 <- sub("\\)","",tmp.G1[[1]][2])
            if(nchar(tmpG1) > 1){
                g1 <- tmpG1           	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            tmp.G2 <- strsplit(as.character(report$G2_Ensembl_HGNC_ID[i]), "\\(")
            G2 <- gsub(" ","",tmp.G2[[1]][1])
            tmpG2 <- sub("\\)","",tmp.G2[[1]][2])
            if(nchar(tmpG1) > 1){
                g2 <- tmpG2 	            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			            strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  paste(G1,as.character(report$G1_exon[i]),sep=":"),
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
					   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  paste(G2,as.character(report$G2_exon[i]),sep=":"),
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="FusionFinder", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
	   return(fusionList)	
}	


#bellerophontes import
.bfImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	    report <- read.table(fusion.report, sep="\t", header=F, comment.char = "", fill=T)
	    tmp1 <- report[seq(1, dim(report)[1], by=2),]
	    names(tmp1) <- c("g1","strand1","g2","strand2", "Chromosome1", "start1","end1","Chromosome2", "start2","end2","supporting.reads")
        tmp2 <- as.character(report[seq(2, dim(report)[1], by=2),1])
        report <- cbind(tmp1, "fusion"=tmp2)
        tmp <- apply(report[,5:11], 1,function(x) {
	            tmp.x <- paste(as.character(x), collapse=":", sep="")
	            tmp.x <- gsub(" ", "", tmp.x)
	    })
		
        
#	    fusionreads.loc <- new("GAlignments")
	     
	    #KnownExonNumber1
		e1 <- NA
		#KnownExonNumber2
		e2 <- NA
		#KnownTranscriptStrand1
		ts1 <- NA
		#KnownTranscriptStrand2
		ts2 <- NA
		#UniqueCuttingPositionCount
		cut <- 0
		#RescuedCount
        rc <- 0
		#frameShift
		fs <- "-"
		#SplicePattern
		sp <- "-"
		#loading annotation
        grHs <- .geneLevelAnnotation(genome=org)		#creating object
		fusionList <- list()
		for(i in 1:dim(report)[1]){
			 #    cat(i)
			 #    cat("\n")
			     sup.reads.tmp <- strsplit(as.character(report$supporting.reads[i]), "reads#: ")
			     sup.reads <- as.numeric(sup.reads.tmp[[1]][2])
				 strand1 <- as.character(report$strand1[i])
				 strand2 <- as.character(report$strand2[i])
				 fs.tmp <- strsplit(as.character(report$fusion[i]),"\\|")
				 fs.1 <- fs.tmp[[1]][1]
				 fs.2 <- fs.tmp[[1]][2]
				 #detecting genes involved in fusions
				 if(report$end2[i] > report$start2[i] && report$end1[i] > report$start1[i]){
		                    grG1 <-  GRanges(seqnames = as.character(report$Chromosome1[i]),
				                     ranges = IRanges(start = report$start1[i], end= report$end1[i]),
				                     strand = strand1)
		                    grG2 <-  GRanges(seqnames = as.character(report$Chromosome2[i]),
							           ranges = IRanges(start = report$start2[i], end= report$end2[i]), 
							           strand = strand2)
				  } else if(report$end2[i] < report$start2[i] && report$end1[i] > report$start1[i]){
					         grG2 <-  GRanges(seqnames = as.character(report$Chromosome1[i]),
							            ranges = IRanges(start = report$start1[i], end= report$end1[i]),
							            strand = strand1)
					         grG1 <-  GRanges(seqnames = as.character(report$Chromosome2[i]),
										ranges = IRanges(start = report$end2[i], end= report$start2[i]), 
										strand = strand2)
				}else if(report$end1[i] < report$start1[i] && report$end2[i] > report$start2[i]){
		        	         grG2 <-  GRanges(seqnames = as.character(report$Chromosome1[i]),
						              ranges = IRanges(start = report$end1[i], end= report$start1[i]),
			                          strand = strand1)
				             grG1 <-  GRanges(seqnames = as.character(report$Chromosome2[i]),
				                      ranges = IRanges(start = report$start2[i], end= report$end2[i]), 
	    				              strand = strand2)
				 }
		         tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
		         if(!is.na(tmpG1)){
		             g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
		         }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
		         tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
		         if(!is.na(tmpG2)){
		             g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
		         }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}
                 if(report$end2[i] > report$start2[i] && report$end1[i] > report$start1[i]){
			        gr1 <- GRanges(seqnames = as.character(report$Chromosome1[i]),
				            ranges = IRanges(start = report$start1[i], end= report$end1[i]),
				            strand = strand1,
				            KnownGene = as.character(g1),
				            KnownTranscript =  as.character(report$g1[i]),
				            KnownExonNumber = e1,
				            KnownTranscriptStrand = ts1,
				            FusionJunctionSequence =  fs.1)
			        gr2 <- GRanges(seqnames = as.character(report$Chromosome2[i]),
						   	ranges = IRanges(start = report$start2[i], end= report$end2[i]), 
							strand = strand2,
						   KnownGene = as.character(g2),
						   KnownTranscript =  as.character(report$g2[i]),
						   KnownExonNumber = e2,
						   KnownTranscriptStrand = ts2,
						   FusionJunctionSequence =  fs.2)
				} else if(report$end2[i] < report$start2[i] && report$end1[i] > report$start1[i]){
					        gr2 <- GRanges(seqnames = as.character(report$Chromosome1[i]),
						            ranges = IRanges(start = report$start1[i], end= report$end1[i]),
						            strand = strand1,
						            KnownGene = as.character(g1),
						            KnownTranscript =  as.character(report$g1[i]),
						            KnownExonNumber = e1,
						            KnownTranscriptStrand = ts1,
						            FusionJunctionSequence =  fs.1)
					        gr1 <- GRanges(seqnames = as.character(report$Chromosome2[i]),
								   	ranges = IRanges(start = report$end2[i], end= report$start2[i]), 
									strand = strand2,
								   KnownGene = as.character(g2),
								   KnownTranscript =  as.character(report$g2[i]),
								   KnownExonNumber = e2,
								   KnownTranscriptStrand = ts2,
								   FusionJunctionSequence =  fs.2)
			  }else if(report$end1[i] < report$start1[i] && report$end2[i] > report$start2[i]){
			               gr2 <- GRanges(seqnames = as.character(report$Chromosome1[i]),
				                  ranges = IRanges(start = report$end1[i], end= report$start1[i]),
				                  strand = strand1,
				                  KnownGene = as.character(g1),
				                  KnownTranscript =  as.character(report$g1[i]),
				                  KnownExonNumber = e1,
				                  KnownTranscriptStrand = ts1,
				                  FusionJunctionSequence =  fs.1)
			               gr1 <- GRanges(seqnames = as.character(report$Chromosome2[i]),
						       	  ranges = IRanges(start = report$start2[i], end= report$end2[i]), 
							      strand = strand2,
						          KnownGene = as.character(g2),
						          KnownTranscript =  as.character(report$g2[i]),
						          KnownExonNumber = e2,
						          KnownTranscriptStrand = ts2,
						          FusionJunctionSequence =  fs.2)
			    }				
				grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
				fusionData <- new("list", fusionTool="bellerophontes", 
				                                 UniqueCuttingPositionCount=cut, 
				                                 SeedCount=sup.reads, 
				                                 RescuedCount=rc, 
				                                 SplicePattern=sp,
				                                 FusionGene=paste(as.character(report$g1[i]),as.character(report$g2[i]), sep="->"),
				                                 frameShift=fs
				)
				fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))		             
			   }
			   return(fusionList)
}
#ChimeraScann import
.csImport <- function(fusion.report, min.support=0, org=c("hg19","mm9","hg38","mm10")){
	    Mmusculus <- NULL	    
	    report <- read.table(fusion.report, sep="\t", header=F, quote ="")
	    names(report) <- c("chrom5p", "start5p", "end5p", "chrom3p", "start3p", "end3p", "chimera_cluster_id", "score", "strand5p", "strand3p", "transcript_ids_5p", "transcript_ids_3p", "genes5p", "genes3p", "type", "distance", "total_frags", "spanning_frags", "unique_alignment_positions", "isoform_fraction_5p", "isoform_fraction_3p", "breakpoint_spanning_reads", "chimera_ids")
		report <- report[which(as.numeric(report$spanning_frags) >= min.support),]
        cat(paste("\n",dim(report)[1]," detected fusions\n",sep=""))
		if(dim(report)[1]==0){
			cat(paste("\n",dim(report)[1]," fusions supported by at least ",min.support," spanning reads",sep=""))
			return()
		}
	#	fusionreads.loc <- new("GAlignments")
	#loading annotation
        grHs <- .geneLevelAnnotation(genome=org) #
	    fusionList <- list()
	    for(i in 1:dim(report)[1]){
		 strand1 <- as.character(report$strand5p[i])
		 strand2 <- as.character(report$strand3p[i])
		 #detecting genes involved in fusions		
         grG1 <-  GRanges(seqnames = as.character(report$chrom5p[i]), ranges = IRanges(start = (as.numeric(report$end5p[i]) - 30), end= as.numeric(report$end5p[i])), strand = strand1)
         grG2 <-  GRanges(seqnames = as.character(report$chrom3p[i]), ranges = IRanges(start = as.numeric(report$start3p[i]), end= (as.numeric(report$start3p[i]) + 30)), strand = strand2)
         tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
         if(!is.na(tmpG1)){
             g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
         }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
         tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
         if(!is.na(tmpG2)){
             g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
         }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}
         if(org=="hg19"){
	         junctionG1 <-  GRanges(seqnames = as.character(report$chrom5p[i]), ranges = IRanges(start = as.numeric(report$start5p[i]), end= as.numeric(report$end5p[i])), strand = strand1)
	         junctionG2 <-  GRanges(seqnames = as.character(report$chrom3p[i]), ranges = IRanges(start = as.numeric(report$start3p[i]), end= as.numeric(report$end3p[i])), strand = strand2)
		     fs.1 <- as.character(getSeq(Hsapiens, junctionG1))
		     fs.2 <- as.character(getSeq(Hsapiens, junctionG2))		
         }else if(org=="mm9"){
	         require(BSgenome.Mmusculus.UCSC.mm9) || stop("\nMissing BSgenome.Mmusculus.UCSC.mm9 library\n")
	         junctionG1 <-  GRanges(seqnames = as.character(report$chrom5p[i]), ranges = IRanges(start = as.numeric(report$start5p[i]), end= as.numeric(report$end5p[i])), strand = strand1)
	         junctionG2 <-  GRanges(seqnames = as.character(report$chrom3p[i]), ranges = IRanges(start = as.numeric(report$start3p[i]), end= as.numeric(report$end3p[i])), strand = strand2) 
		     fs.1 <- as.character(getSeq(Mmusculus, junctionG1))
		     fs.2 <- as.character(getSeq(Mmusculus, junctionG2))		
         }
		 if(length(as.character(report$transcript_ids_5p[i])) == 0){
			tmpT1 <- NULL
		 } else{
		   	tmpT1 <- as.character(report$transcript_ids_5p[i])	
		 } 
		 if(length(as.character(report$transcript_ids_3p[i])) == 0){
			tmpT2 <- NULL
		 } else{
		   	tmpT2 <- as.character(report$transcript_ids_3p[i])	
		 } 
         #exon number
			tmpEn1 <- ""
			tmpEn2 <- ""
		# transcript strand
			tmpTs1 <- ""
			tmpTs2 <- ""

	     gr1 <- GRanges(seqnames = as.character(report$chrom5p[i]), ranges = IRanges(start = (as.numeric(report$end5p[i]) - 30), end= as.numeric(report$end5p[i])),
		            strand = strand1,
		            KnownGene = as.character(g1),
		            KnownTranscript =  tmpT1,
		            KnownExonNumber = tmpEn1,
		            KnownTranscriptStrand = tmpTs1,
		            FusionJunctionSequence =  fs.1)
	     gr2 <- GRanges(seqnames = as.character(report$chrom3p[i]), ranges = IRanges(start = as.numeric(report$start3p[i]), end= (as.numeric(report$start3p[i]) + 30)),
				   strand = strand2,
				   KnownGene = as.character(g2),
				   KnownTranscript =  tmpT2,
				   KnownExonNumber = tmpEn2,
				   KnownTranscriptStrand = tmpTs2,
				   FusionJunctionSequence =  fs.2)
		grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
		fusionData <- new("list", fusionTool="FusionMap", 
		                                 UniqueCuttingPositionCount=report$unique_alignment_positions[i], 
		                                 SeedCount=report$spanning_frags[i], 
		                                 RescuedCount=report$total_frags[i], 
		                                 SplicePattern="",
		                                 FusionGene=paste(as.character(report$genes5p[i]),report$genes3p[i], sep=":"),
		                                 frameShift=""
		)
		fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))		             
	   }
	   return(fusionList)
}


.buildFusion <- function(type=c("donor.end","acceptor.start"), fusion.grl, tx.id){
    eg.lst <- list("tx_id" = tx.id)
    eg.trs.e <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene, vals=eg.lst, columns=c("tx_id","exon_id","exon_rank"))
    #handling the 5' end of the fusion
    if(type=="donor.end"){
	    donor.intron <- NA
        exons.tx <- elementMetadata(eg.trs.e)$tx_id
        exons.tx <- as.list(exons.tx)
        exons.rank <- elementMetadata(eg.trs.e)$exon_rank
        exons.idx<- sapply(exons.tx,function(x,y){grep(y, x)},y=tx.id)
        rank.e <- NULL
        for(i in 1:length(exons.idx)){
           rank.e[i] <- exons.rank[[i]][exons.idx[i]]
        }
        fusion.pos <- findOverlaps(fusion.grl[[1]],  eg.trs.e, type = "any", select = "first", ignore.strand = T)
        if(!is.na(fusion.pos)){
          end(eg.trs.e[fusion.pos]) <- end(fusion.grl[[1]])
          if(unique(as.character(strand(eg.trs.e))) == "-"){
             eg.trs.e <- eg.trs.e[fusion.pos:length(eg.trs.e)]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             eg.trs.rnk <- order(eg.trs.rnk,decreasing=T)
             donor.seq <- NULL
             for(i in eg.trs.rnk){
                donor.seq <- c(donor.seq, as.character(eg.trs.seq[i]))
                donor.seq <- paste(donor.seq, collapse="")
             }             
          }else{
             eg.trs.e <- eg.trs.e[1:fusion.pos]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             donor.seq <- NULL
             for(i in eg.trs.rnk){
                donor.seq <- c(donor.seq, as.character(eg.trs.seq[i]))
                donor.seq <- paste(donor.seq, collapse="")
             }
          }
        }else{
	        #first intron start between exFIRST and exFIRST+1
	        #last intron between exLAST-1 and exLAST
            start.i <- end(eg.trs.e)
            end.i <- start(eg.trs.e[2:length(eg.trs.e)])
            #adding a fake intron at the end to keep the same organization of the exons
            end.i <- c(end.i, start.i[length(start.i)])
            names.introns <- paste("I",paste(rank.e, (c(rank.e[2:length(rank.e)],rank.e[length(rank.e)])),sep="-"),sep="")
            eg.trs.i <- GRanges(seqnames =seqnames(eg.trs.e),ranges = IRanges(start=start.i, end= end.i, names = names.introns), strand = strand(eg.trs.e))
            fusion.pos <- findOverlaps(fusion.grl[[1]],  eg.trs.i, type = "any", select = "first", ignore.strand = T)
            donor.intron <- names(eg.trs.i[fusion.pos])
            my.intron <- sub("I","",names(eg.trs.i[fusion.pos]))
            my.intron.lst <- strsplit(my.intron,"-")
            my.intron <- as.numeric(my.intron.lst[[1]][1])
            end(eg.trs.e[my.intron]) <- end(fusion.grl[[1]])
            if(unique(as.character(strand(eg.trs.e))) == "-"){
             eg.trs.e <- eg.trs.e[fusion.pos:length(eg.trs.e)]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             eg.trs.rnk <- order(eg.trs.rnk,decreasing=T)
             donor.seq <- NULL
             for(i in eg.trs.rnk){
                donor.seq <- c(donor.seq, as.character(eg.trs.seq[i]))
                donor.seq <- paste(donor.seq, collapse="")
             }             
          }else{
             eg.trs.e <- eg.trs.e[1:fusion.pos]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             donor.seq <- NULL
             for(i in eg.trs.rnk){
                donor.seq <- c(donor.seq, as.character(eg.trs.seq[i]))
                donor.seq <- paste(donor.seq, collapse="")
             }
          }
        }
        donor.seq <- DNAString(donor.seq)
        return(list(seq=donor.seq, intron.location=donor.intron))
    }else if(type=="acceptor.start"){
	    acceptor.intron <- NA
        exons.tx <- elementMetadata(eg.trs.e)$tx_id
        exons.tx <- as.list(exons.tx)
        exons.rank <- elementMetadata(eg.trs.e)$exon_rank
        exons.idx<- sapply(exons.tx,function(x,y){grep(y, x)},y=tx.id)
        rank.e <- NULL
        for(i in 1:length(exons.idx)){
           rank.e[i] <- exons.rank[[i]][exons.idx[i]]
        }
        fusion.pos <- findOverlaps(fusion.grl[[2]],  eg.trs.e, type = "any", select = "first", ignore.strand = T)
        if(!is.na(fusion.pos)){
          start(eg.trs.e[fusion.pos]) <- start(fusion.grl[[2]])
          if(unique(as.character(strand(eg.trs.e))) == "-"){
             eg.trs.e <- eg.trs.e[1:fusion.pos]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             eg.trs.rnk <- order(eg.trs.rnk,decreasing=T)
             accept.seq <- NULL
             for(i in eg.trs.rnk){
                accept.seq <- c(accept.seq, as.character(eg.trs.seq[i]))
                accept.seq <- paste(accept.seq, collapse="")
             }             
          }else{
             eg.trs.e <- eg.trs.e[fusion.pos:length(eg.trs.e)]
             eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
             eg.trs.rnk <- seq(1,length(eg.trs.seq))
             accept.seq <- NULL
             for(i in eg.trs.rnk){
                accept.seq <- c(accept.seq, as.character(eg.trs.seq[i]))
                accept.seq <- paste(accept.seq, collapse="")
             }
          }
        }else{
             #first intron start between exFIRST and exFIRST+1
	        #last intron between exLAST-1 and exLAST
            start.i <- end(eg.trs.e)
            start.i <- start.i + 1
            end.i <- start(eg.trs.e[2:length(eg.trs.e)])
            #adding a fake intron at the end to keep the same organization of the exons
            end.i <- c(end.i, start.i[length(start.i)])
            end.i <- end.i - 1
            names.introns <- paste("I",paste(rank.e, (c(rank.e[2:length(rank.e)],rank.e[length(rank.e)])),sep="-"),sep="")
            eg.trs.i <- GRanges(seqnames =seqnames(eg.trs.e),ranges = IRanges(start=start.i, end= end.i, names = names.introns), strand = strand(eg.trs.e))
            fusion.pos <- findOverlaps(fusion.grl[[2]],  eg.trs.i, type = "any", select = "first", ignore.strand = T)
            acceptor.intron <- names(eg.trs.i[fusion.pos])
            my.intron <- sub("I","",names(eg.trs.i[fusion.pos]))
            my.intron.lst <- strsplit(my.intron,"-")
            my.intron <- as.numeric(my.intron.lst[[1]][2])
            #changing the exon end with the intron end
            end(eg.trs.e[my.intron]) <- end(eg.trs.i[my.intron])
            start(eg.trs.e[my.intron]) <- start(fusion.grl[[2]])
            if(unique(as.character(strand(eg.trs.e))) == "-"){
               #right!
               eg.trs.e <- eg.trs.e[1:fusion.pos]
               eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
               eg.trs.rnk <- seq(1,length(eg.trs.seq))
               eg.trs.rnk <- order(eg.trs.rnk,decreasing=T)
               accept.seq <- NULL
               for(i in eg.trs.rnk){
                  accept.seq <- c(accept.seq, as.character(eg.trs.seq[i]))
                  accept.seq <- paste(accept.seq, collapse="")
               }             
             }else{
               eg.trs.e <- eg.trs.e[fusion.pos:length(eg.trs.e)]
               eg.trs.seq <- getSeq(Hsapiens, eg.trs.e)
               eg.trs.rnk <- seq(1,length(eg.trs.seq))
               accept.seq <- NULL
               for(i in eg.trs.rnk){
                  accept.seq <- c(accept.seq, as.character(eg.trs.seq[i]))
                  accept.seq <- paste(accept.seq, collapse="")
               }
            }
          }
        accept.seq <- DNAString(accept.seq)
        return(list(seq=accept.seq, intron.location=acceptor.intron))
    }
}
#############
.starImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10"), min.support=10){	        
    tmp.file <- gsub(" ","_",date()) 
    tmp.file <- sub("__","_",tmp.file)
    tmp.file <- gsub(":","-",tmp.file) 
	
    .C("StarParser",2,c("program",fusion.report,tmp.file), PACKAGE = 'chimera')

   report <- read.table(tmp.file, sep="\t", header=TRUE)
    names(report) <- c("gene1.chr","gene1.fusion.loc","gene1.strand","gene2.chr","gene2.fusion.loc","gene2.strand","n.spanning","n.encompassing")
	#removing chrM
    if(length(which(as.character(report$gene1.chr)=="chrM"))>0){
         cat("\nchrM is removed from fusion acceptor\n")
    }
    report <- report[setdiff(seq(1:dim(report)[1]),which(as.character(report$gene1.chr)=="chrM")),]
    if(length(which(as.character(report$gene2.chr)=="chrM"))>0){
         cat("\nchrM is removed from fusion donor\n")
    }
    report <- report[setdiff(seq(1:dim(report)[1]),which(as.character(report$gene2.chr)=="chrM")),]
    #removing non canonical chrs1
    chr.g1.l <- sapply(as.character(report$gene1.chr), nchar)
    if(length(which(as.numeric(chr.g1.l) > 5)) >0){
             removed <- as.character(unique(report$gene1.chr[which(as.numeric(chr.g1.l) > 5)]))
             cat("\nThe following chrs were removed from fusion acceptor:\n",removed,"\n")
    }
    report <- report[which(as.numeric(chr.g1.l) <= 5),]
    chr.g2.l <- sapply(as.character(report$gene2.chr), nchar)
    if(length(which(as.numeric(chr.g2.l) > 5)) >0){
             removed <- as.character(unique(report$gene2.chr[which(as.numeric(chr.g2.l) > 5)]))
             cat("\nThe following chrs were removed from fusion donor:\n",removed,"\n")
    }
	
	#just a reminder of the star file structure
    #the encompassing chimeric reads are marked with -1 in column junction.type, 
    #while spanning reads are marked with 0 for non-GT/AG introns, 1-GT/AG, 2-CT/AC.
	
    report <- report[which(report$n.spanning >= min.support),]
    if(sum(report$n.spanning) == 0){
		cat("\nThe input file does not have any spanning read.\nYour fusion lacking of spanning reads are most probably artifacts\nThe analysis of fusions lacking spanning reads is not supported.\n")
        return()
	}
	cat(paste("\nImporting ",dim(report)[1]," fusions\n", sep=""))
    #loading annotation
    
    grHs <- .geneLevelAnnotation(genome=org)

	#    fusionreads.loc <- new("GAlignments") 
        #KnownGene1
        g1 <- NA
		#KnownGene2
        g2 <- NA
        #KnownTranscript1
        t1 <- NA
        #KnownTranscript2
        t2 <- NA
        #KnownExonNumber1
        e1 <- NA
        #KnownExonNumber2
        e2 <- NA
        #KnownTranscriptStrand1
        ts1 <- NA
        #KnownTranscriptStrand2
        ts2 <- NA
        #FusionJunctionSequence1-2
        #SeedCount1-2
        #UniqueCuttingPositionCount
        cut <- 0
        #SplicePattern
        sp <- "-"
        #frameShift
        fs <- "-"
        #creating objects
        fusionList <- list()
	    cat("\n")
        for(i in 1:dim(report)[1]){
			cat(".")
	        strand1 <- as.character(report$gene1.strand[i])
		    strand2 <- as.character(report$gene2.strand[i])
            #spanning reads		 
	        seed1 <- report$n.spanning[i]
	        #RescuedCount encompassing.reads
	        rc <- report$n.encompassing[i]
	        #pos.tmp <- strsplit(as.character(report$chrs.fusion[i]), "-")
	        #junction sequence1
	        fs1 <- ""
	        #defining start1
	        #coverage1.tmp <- strsplit(as.character(report$coverage1[i]), " ")
            start1 <- report$gene1.fusion.loc[i] - 30
            end1 <- report$gene1.fusion.loc[i]
            chr1 <- as.character(report$gene1.chr[i])
	        #junction2
		    fs2 <- ""
		    #defining start2
		    #coverage2.tmp <- strsplit(as.character(report$coverage2[i]), " ")
	        start2 <- report$gene2.fusion.loc[i]
            end2 <- report$gene2.fusion.loc[i] + 30
            chr2 <- as.character(report$gene2.chr[i])
            #defining gene1
            grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1))
            grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2))
            tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG1)){
                g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
            }else{g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")}
            tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
            if(!is.na(tmpG2)){
                g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
            }else{g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")}

            gr1 <- GRanges(seqnames = chr1,
			            ranges = IRanges(start = start1, end= end1),
			            strand = strand1,
			            KnownGene = as.character(g1),
			            KnownTranscript =  t1,
			            KnownExonNumber = e1,
			            KnownTranscriptStrand = ts1,
			            FusionJunctionSequence =  fs1)
		     gr2 <- GRanges(seqnames = chr2,
					   ranges = IRanges(start = start2, end= end2),
					   strand = strand2,
					   KnownGene = as.character(g2),
					   KnownTranscript =  t2,
					   KnownExonNumber = e2,
					   KnownTranscriptStrand = ts2,
					   FusionJunctionSequence =  fs2)
			grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
			fusionData <- new("list", fusionTool="Star", 
			                                 UniqueCuttingPositionCount=cut, 
			                                 SeedCount=seed1, 
			                                 RescuedCount=rc, 
			                                 SplicePattern=sp,
			                                 FusionGene=paste(g1,g2, sep="->"),
			                                 frameShift=fs
			)
			fusionList[[i]] <- new("fSet",fusionInfo=fusionData,fusionLoc=grl, fusionRNA=new("DNAStringSet"))
        }
		cat("\n")
	   return(fusionList)
}
####fusionCatcher
.fcImport <- function(fusion.report, org=c("hg19","mm9","hg38","mm10")){
	report <- read.table(fusion.report, sep="\t", header=TRUE)
	
	#KnownTranscript1
	t1 <- NA
	#KnownTranscript2
	t2 <- NA
	#KnownExonNumber1
	e1 <- NA
	#KnownExonNumber2
	e2 <- NA
	#KnownTranscriptStrand1
	ts1 <- NA
	#KnownTranscriptStrand2
	ts2 <- NA

    grHs <- .geneLevelAnnotation(genome=org)	
	
	fusionList <- list()
	for(i in 1:dim(report)[1]){
		fs.tmp <- strsplit(as.character(report$Fusion_sequence[i]), split="\\*")
		fs.1 <- fs.tmp[[1]][1]
		fs.2 <- fs.tmp[[1]][2]

		gene1_pos <- unlist(strsplit(as.character(report$Fusion_point_for_gene_1.5end_fusion_partner.[i]), split=":"))
		gene2_pos <- unlist(strsplit(as.character(report$Fusion_point_for_gene_2.3end_fusion_partner.[i]), split=":"))
		grG1 <- GRanges(seqnames = paste("chr", gene1_pos[1], sep=""), ranges = IRanges(start = (as.numeric(gene1_pos[2]) - nchar(fs.1)), end = as.numeric(gene1_pos[2])), strand = gene1_pos[3])
		grG2 <- GRanges(seqnames = paste("chr", gene2_pos[1], sep=""), ranges = IRanges(start = as.numeric(gene2_pos[2]), end = (as.numeric(gene2_pos[2]) + nchar(fs.2))), strand = gene2_pos[3])
		tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
		if(!is.na(tmpG1)){
			g1 <- elementMetadata(grHs[tmpG1])$symbol
		}else{
			g1 <- paste(seqnames(grG1), paste(start(grG1), end(grG1), sep="-"), sep=":")
		}
		tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand=T)
		if(!is.na(tmpG2)){
			g2 <- elementMetadata(grHs[tmpG2])$symbol
		}else{
			g2 <- paste(seqnames(grG2), paste(start(grG2), end(grG2), sep="-"), sep=":")
		}

		gr1 <- GRanges(seqnames = paste("chr", gene1_pos[1], sep=""),
			ranges = IRanges(start = (as.numeric(gene1_pos[2]) - nchar(fs.1)), end = as.numeric(gene1_pos[2])),
			strand = gene1_pos[3],
			KnownGene = as.character(g1),
			KnownTranscript = t1,
			KnownExonNumber = e1,
			KnownTranscriptStrand = ts1,
			FusionJunctionSequence = fs.1)
		gr2 <- GRanges(seqnames = paste("chr", gene2_pos[1], sep=""),
			ranges = IRanges(start = as.numeric(gene2_pos[2]), end = (as.numeric(gene2_pos[2]) + nchar(fs.2))),
			strand = gene2_pos[3],
			KnownGene = as.character(g2),
			KnownTranscript = t2,
			KnownExonNumber = e2,
			KnownTranscriptStrand = ts2,
			FusionJunctionSequence = fs.2)
		grl <- GRangesList("gene1" = gr1, "gene2" = gr2)

		fusionData <- new("list", fusionTool="FusionCatcher",
			UniqueCuttingPositionCount = NA,
			SeedCount = report$Spanning_unique_reads[i],
			RescuedCount = report$Spanning_pairs[i],
			SplicePattern = NA,
			FusionGene = paste(report$Gene_1_symbol.5end_fusion_partner.[i], "->", report$Gene_2_symbol.3end_fusion_partner., sep=""),
			frameShift = NA)
		fusionList[[i]] <- new("fSet", fusionInfo=fusionData, fusionLoc=grl, fusionRNA=new("DNAStringSet"))

	}
    return(fusionList)
}



