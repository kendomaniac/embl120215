

###
filterList <- function(x,type=c("spanning.reads","fusion.names", "intronic", "annotated.genes", "read.through", "oncofuse"),oncofuse.output=NULL, query=NULL, oncofuse.type=c("g5CDS", "g3CDS", "g5g3CDS", "g5exon", "g3exon", "g5g3exon", "passenger.prob", "expression.gain"), parallel=FALSE){
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
       if(type=="spanning.reads"){
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
			     tmp.rm <- strsplit(tmp, ":")
			     nsr <- sapply(tmp.rm, function(x){
			              if(x[1]==x[2]){
				              return(1)
			              }else{return(0)}
			     })
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
	   if(type=="oncofuse"){
		   if(length(oncofuse.type)!=1){
		   	   cat("\nYou need to select one option for oncofuse filtering\n")
		   }
		  if(dim(oncofuse.output)[2]!=34){
			  cat("\nIt seems that the output of oncofuse has not the correct number of columns\n")
			  return()
		  }
		  if(length(x)<1){
			  cat("\nIt seems that the list of fusions you provided is empty\n")
			  return()
		  }
          if(parallel){
		     require(BiocParallel) || stop("\nMission BiocParallel library\n")
   	         p <- MulticoreParam()
	         tmp <- fusionName(x, parallel=T)
	      }else{
		     tmp <- fusionName(x)
	      }
		  of.fn <- paste(as.character(oncofuse.output$X5_FPG_GENE_NAME),as.character(oncofuse.output$X3_FPG_GENE_NAME), sep=":")
          if(length(intersect(tmp, of.fn))==0){
          	  cat("\nIt seems that there is no overlaps between oncofuse output and the fusion names of the fusions list.\n")
			  return()
          }
		  cat("\n", "The number of fusions detected by oncofuse are ", length(intersect(tmp, of.fn)),"\n")
		  if(oncofuse.type=="g5CDS"){
		  	loc <- intersect(tmp, of.fn[which(oncofuse.output$X5_IN_CDS.=="Yes")])
			cat("\n", "The number of fusions detected by oncofuse into CDS of gene at 5' end are ", length(intersect(tmp, of.fn)),"\n")
		  }
		  if(oncofuse.type=="g3CDS"){
		  	loc <- intersect(tmp, of.fn[which(oncofuse.output$X3_IN_CDS.=="Yes")])
			cat("\n", "The number of fusions detected by oncofuse into CDS of gene at 3' end are ", length(loc),"\n")
		  }
		  if(oncofuse.type=="g5g3CDS"){
		  	loc <- intersect(tmp, of.fn[intersect(which(oncofuse.output$X5_IN_CDS.=="Yes"), which(oncofuse.output$X3_IN_CDS.=="Yes"))])
			cat("\n", "The number of fusions detected by oncofuse into CDS of both genes are ", length(loc),"\n")
		  }
		  if(oncofuse.type=="g5exon"){
		  	loc <- intersect(tmp, of.fn[which(oncofuse.output$X5_SEGMENT_TYPE=="Exon")])
			cat("\n", "The number of fusions detected by oncofuse into an exon of gene at 5' end are ", length(intersect(tmp, of.fn)),"\n")
		  }
		  if(oncofuse.type=="g3exon"){
		  	loc <- intersect(tmp, of.fn[which(oncofuse.output$X3_SEGMENT_TYPE=="Exon")])
			cat("\n", "The number of fusions detected by oncofuse into an exon of gene at 3' end are ", length(loc),"\n")
		  }
		  if(oncofuse.type=="g5g3exon"){
		  	loc <- intersect(tmp, of.fn[intersect(which(oncofuse.output$X5_SEGMENT_TYPE=="Exon"), which(oncofuse.output$X3_SEGMENT_TYPE=="Yes"))])
			cat("\n", "The number of fusions detected by oncofuse into an exon of both genes are ", length(loc),"\n")
		  }
		  if(oncofuse.type=="passenger.prob"){
		  	loc <- intersect(tmp, of.fn[which(as.numeric(oncofuse.output$P_VAL_CORR)<=query)])
			cat("\n", "The number of fusions detected by oncofuse with a passenger probability lower then ",query, " are ", length(loc),"\n")
		  }
		  if(oncofuse.type=="expression.gain"){
		  	loc <- intersect(tmp, of.fn[which(as.numeric(oncofuse.output$EXPRESSION_GAIN)>=query)])
			cat("\n", "The number of fusions detected by oncofuse with a expression gain score greater then ",query, " are ", length(loc),"\n")
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
		tmp.seq <- .onlyExons(type="donor.end", grl, tmp.tx[i])
        tmp.gene1 <- c(tmp.gene1, tmp.seq$seq)
        if(!is.na(tmp.seq$intron.location)){
	           return(0)#junction in an exon
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
		  tmp.seq <- .onlyExons(type="acceptor.start", grl, tmp.tx[i])
          tmp.gene2 <- c(tmp.gene2, tmp.seq$seq)
          if(!is.na(tmp.seq$intron.location)){
	           return(0)#junction in an exon
        }
  }
  return(1)#junction in an intron
}


#a function to identify the presence of intronic region at a junction point
.onlyExons <- function(type=c("donor.end","acceptor.start"), fusion.grl, tx.id){
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
        if(is.na(fusion.pos)){
			return(list(seq=NA, intron.location=NA))#junctionin an intron
		}else{
			return(list(seq="OK", intron.location="OK"))
		}             
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
        if(is.na(fusion.pos)){
			return(list(seq=NA, intron.location=NA))#junctionin an intron
		}else{
			return(list(seq="OK", intron.location="OK"))
		}
	}
}





supportingReads <- function(list, fusion.reads=c("encompassing","spanning"), parallel=FALSE){
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
	    if(fusion.reads=="encompassing"){
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

#oncofuse installation
oncofuseInstallation <- function(){
	    mydir <- getwd()
		oncofuseDirLocation  <- paste(path.package("chimera", quiet = FALSE), "/oncofuse", sep="")
		tmp.info <- dir(paste(path.package("chimera", quiet = FALSE)))
		if(length(grep("oncofuse", tmp.info))>0){
		   unlink(oncofuseDirLocation, recursive = T, force = T)
	    }
		dir.create(oncofuseDirLocation, showWarnings = TRUE, recursive = FALSE)
		setwd(oncofuseDirLocation)
		cat("\nBegin downloads of oncofuse.....\n")
		download.file("https://github.com/mikessh/oncofuse/releases/download/1.0.9/oncofuse-1.0.9.zip", "oncofuse.zip", mode="wb")
        cat("\nThe oncofuse downloaded version is oncofuse-v1.0.9\n")
        system(paste("unzip oncofuse.zip", sep=""))
        system("cp -fR ./oncofuse-v1.0.9/* .")
        system("rm -fR ./oncofuse-v1.0.9/")
        system("chmod +x *")
        setwd(mydir)
        return()
}
oncofuseRun <- function(listfSet, tissue=c("EPI","HEM","MES","AVG"), org=c("hg19","hg38"), threads=1, plot=FALSE){
	oncofuseDirLocation  <- paste(path.package("chimera", quiet = FALSE), "/oncofuse", sep="")
	of.input <- paste("of",gsub("[' '| :]","-", date()),sep="_")
	extract.loc <-function(fset, tissue){
	    chr.g1 <- as.character(seqnames(fusionGRL(fset)$gene1))
		chr.g2 <- as.character(seqnames(fusionGRL(fset)$gene2))
		end.g1 <- end(fusionGRL(fset)$gene1)
		start.g2 <- start(fusionGRL(fset)$gene2)
		return(list(chr.g1, end.g1, chr.g2, start.g2, tissue))
    }
	fset.of <- sapply(listfSet, extract.loc, tissue)
	fset.of <- as.data.frame(fset.of)
	fset.of <- t(fset.of)
	write.table(fset.of, of.input, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	cat("\nStart Oncofuse analysis\n")
	system(paste("java -Xmx1G -jar ",paste(oncofuseDirLocation,"/Oncofuse.jar ",sep=""), of.input," coord ",tissue," ", sub("of","of.out",of.input)," -a ",org," -p ",threads, sep=""))
    of.out <-read.table(sub("of","of.out",of.input), sep="\t", header=TRUE, fill=TRUE, na.strings="NaN", quote = "")
	cat("\nEnd Oncofuse analysis\n")
	return(of.out)
	if(plot){
	    smoothScatter(log10(as.numeric(of.out$P_VAL_CORR)), log10(as.numeric(of.out$EXPRESSION_GAIN)), xlab="log10 PASSENGER PROBABILITY", ylab="log10 EXPRESSION GAIN", pch=19, cex=0.5, col="red")
    }
}
#Bayesian probability of fusion being a passenger (class 0), given as Bonferroni-corrected P-value
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
bam2fastq<- function(bam, filename="ready4gapfiller", ref,parallel=FALSE){
    bam <- scanBam(bam)
    tmp <- which(as.character(bam[[1]]$rname) == ref)
    seq <- as.character(bam[[1]]$seq[tmp])
    qual <- as.character(bam[[1]]$qual[tmp])
    read.n <- as.character(bam[[1]]$qname[tmp])
    names(seq) <- read.n 
    names(qual) <- read.n 
    seq <- seq[order(names(seq))]
    qual <- qual[order(names(qual))]
    rname.u <- unique(names(seq))
    if(parallel){
 	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
          p <- MulticoreParam()
          paired <- bplapply(rname.u, function(x,y){
   	          if(length(which(y==x)) == 2){
   	         	return(paste(which(y==x),sep=" "))
   	          }else{
   	        	return(NA)
   	          }
          }, y=names(seq), BPPARAM=p)
		  
    }else{
       paired <- sapply(rname.u, function(x,y){
	       if(length(which(y==x)) == 2){
	         	return(paste(which(y==x),sep=" "))
	       }else{
	        	return(NA)
	       }
       }, y=names(seq))
   }
	
    paired <- paired[!is.na(paired)]
    paired <- as.data.frame(paired)
    paired1 <-as.numeric(as.character(unlist(paired[1,])))
    paired2 <-as.numeric(as.character(unlist(paired[2,])))

    fastq1 <- rep(NA, (length(paired1) * 4))
    fastq1[seq(1, length(fastq1),by=4)] <- paste("@",names(seq)[paired1], sep="")
    fastq1[seq(2, length(fastq1),by=4)] <- seq[paired1]
    fastq1[seq(3, length(fastq1),by=4)] <- rep("+", length(seq[paired1]))
    fastq1[seq(4, length(fastq1),by=4)] <- qual[paired1]
    writeLines(fastq1, paste(filename,"_R1.fastq",sep=""))

    fastq2 <- rep(NA, (length(paired2) * 4))
    fastq2[seq(1, length(fastq2),by=4)] <- paste("@",names(seq)[paired2], sep="")
    fastq2[seq(2, length(fastq2),by=4)] <- seq[paired2]
    fastq2[seq(3, length(fastq2),by=4)] <- rep("+", length(seq[paired2]))
    fastq2[seq(4, length(fastq2),by=4)] <- qual[paired2]
    writeLines(fastq2, paste(filename,"_R2.fastq",sep=""))
}
##
MHmakeRandomString <- function()
{
    n <- 1
	lenght <- 12
	randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                 lenght, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}

##
.plotCoverage <- function (x, x1, start = 1, end = length(x), col = "yellow", col1="blue", xlab = "Index", ylab = "Coverage", breack.point=NULL, ylim=NULL){
     xWindow <- as.vector(window(x, start, end))
	 x1Window <- as.vector(window(x1, start, end))
     x <- start:end
     xlim <- c(start, end)
	 if(length(ylim)==0){
        ylim <- c(0, max(x1Window))
	 }
     plot(x = start, y = 0, xlim = xlim, ylim = ylim, xlab = xlab,ylab = ylab, type = "n")
	 polygon(c(start, x, end), c(0, x1Window, 0), col = col1)
	 polygon(c(start, x, end), c(0, xWindow, 0), col = col)
	 abline(v=breack.point,lty=2, lwd=2, col="red")
}						

breakpointOverlaps <- function(fset, plot=FALSE, ylim=NULL){
 if(length(fusionRNA(fset))==0){
	 cat("\nMissinng fusion transcript information. See chimeraSeqs and addRNA functions\n")
	 return()
 }
 if(length(fusionGA(fset))==0){
	 cat("\nMissinng mapping reads information. See subreadRun and addGA functions\n")
	 return()
 }	 	
 tmp <- names(fusionRNA(fset))
 tmp <- strsplit(tmp, ":")
 tmp <- tmp[[1]][1]
 tmp <- strsplit(tmp, "-")
 #length of tr1 insert till breakpoint
 tmp <- as.numeric(tmp[[1]][2])
 subj.ga <- fusionGA(fset)
 tmp.ga <- GRanges(seqnames = as.character(seqnames(subj.ga[1])), ranges = IRanges(start = tmp, width= 1), strand = "+")
 seqlengths(tmp.ga) <- seqlengths(subj.ga)
 breakpoints <- findOverlaps(query=tmp.ga, subject=subj.ga,type = "within",select = "all", ignore.strand = TRUE)
 bp.ga <- subj.ga[subjectHits(breakpoints)]
 tmp1.cov <- coverage(subj.ga)[[1]]
 tmp.cov <- coverage(bp.ga)[[1]]
 if(plot){
      .plotCoverage(tmp.cov, tmp1.cov, breack.point=tmp, ylim=ylim)
 }
 return(bp.ga)
}     
