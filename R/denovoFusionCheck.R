gapfillerRun <- function(fusion.fa, seed1, seed2, gapfiller=NULL, seed.ins=200, seed.var=50, block.length=5, overlap=20, max.length=5000, slack=7, k=6, global.mismatch=5, perc.identity=0.6){
    alignments <- list()
    if(length(gapfiller)==0){
        chimeraDirLocation  <- path.package("chimera", quiet = FALSE)
        gapfiller.check <- grep("gapfiller",dir(chimeraDirLocation))
        if(length(gapfiller.check)==0){
	         cat("\nYoun need to provide the full path to GapFiller or you can run gapfillerInstallation function\n")
	         return()
        }	else{
			   gapfiller <- paste(path.package("chimera", quiet = FALSE),"/gapfiller/src/GapFiller",sep="")
	    }
    }
    tempdir <- tempdir()
    mydir <- getwd()
    system(paste("cp ",fusion.fa," ",seed1," ", seed2," ",tempdir, sep=""))
    setwd(tempdir)
    system(paste("cat",seed1,seed2,"> r1r2.fasta",sep=" ")) 
    system(paste(gapfiller," --no-read-cycle --output contigs.fasta --statistics contigs.stats --seed1 ",seed1,
                  " --seed2 ",seed2," --query r1r2.fasta --seed-ins ",seed.ins," --seed-var ",seed.var,
                  " --block-length ",block.length," --overlap  ",overlap," --read-length 200 ext.Thr 1 --max-length ",
                  max.length," --slack ",slack," --k ",k," --global-mismatch ",global.mismatch," --perc-identity ",perc.identity, sep=""))

#   gapfiller.stats <- read.table("contigs.stats", sep = " ", header = T)
#	gapfiller.stats <- t(gapfiller.stats)
#	cat("\n")
#	print(gapfiller.stats)
#	cat("\n")
    tmp <- readDNAStringSet("contigs.fasta", format="fasta")
    target <- readDNAStringSet(fusion.fa, format="fasta")
    alignment <- pairwiseAlignment(pattern=tmp, subject=target[1], type="local")
    alignment.view <- Views(alignment)
    tmp1 <- strsplit(names(target),":")
    tmp1.1 <- strsplit(tmp1[[1]][1],"\\.")
    tmp1.2 <- strsplit(tmp1.1[[1]][length(tmp1.1[[1]])],"-")
    threshold <- as.numeric(tmp1.2[[1]][2]) - as.numeric(tmp1.2[[1]][1])
    tmp.acceptor <- which(start(alignment.view) < threshold)
    tmp.donor <- which(end(alignment.view) > threshold)
    alignments <- alignment[intersect(tmp.acceptor, tmp.donor)]
    stat.direct <- length(intersect(tmp.acceptor, tmp.donor))
     if(max(stat.direct) > 0) {
	   cat("\nde novo alignment has overlap over the fusion break point\n")
	   setwd(mydir)
	   return(list(contigs=alignments,junction.contigs=DNAStringSet(Views(alignments)),fusion=target))
    }else{
	   	cat("\nde novo alignment has no overlap over the fusion break point\n")
	    setwd(mydir)
		return()
    }    
}
####
.gfWrap<-function(chimeraSeq.out, bam.file, parallel=FALSE){
	cat(".")
	tmp.file <- MHmakeRandomString()
    bam2fastq(bam=bam.file, filename=tmp.file, ref=names(chimeraSeq.out),parallel=parallel)
    tmp.fa <- paste(MHmakeRandomString(),".fa", sep="")
    writeXStringSet(chimeraSeq.out, tmp.fa, format="fasta")
    mylist  <- gapfillerRun(tmp.fa, seed1=paste(tmp.file,"_R1.fastq",sep=""),  
    seed2=paste(tmp.file,"_R2.fastq",sep=""), gapfiller=NULL, seed.ins=200, 
    seed.var=50, block.length=5, overlap=20, max.length=5000, 
    slack=7, k=6, global.mismatch=5, perc.identity=0.6)
	return(mylist)
}
gapfillerWrap <- function(chimeraSeqSet.out, bam, parallel=c(FALSE,TRUE)){
    if(parallel[1]){
 	     require(BiocParallel) || stop("\nMission BiocParallel library\n")
          p <- MulticoreParam()
          mylist <- bplapply(chimeraSeqSet.out, .gfWrap, bam, parallel=parallel[2], BPPARAM=p)
    }else{
        mylist <- lapply(chimeraSeqSet.out, .gfWrap, bam, parallel=parallel[2])
   }
   return(mylist)	
}

####
gapfillerInstallation <- function(os=c("mac64","unix64")){
	mydir <- getwd()
	chimeraDirLocation  <- path.package("chimera", quiet = FALSE)
	setwd(chimeraDirLocation)
	if(os=="mac64"){
		cat("\nBegin downloads of GapFiller.....\n")
		download.file("http://sourceforge.net/projects/ochguiextras/files/chimera/GapFiller/gapfiller_mac64.zip/download", "gapfiller.zip", mode="wb")
        cat("\nThe GapFiller downloaded version is GapFiller_1_0\n")
        system(paste("unzip gapfiller.zip", sep=""))
        system("chmod -fR +x gapfiller")
	}else if(os=="unix64"){
			cat("\nBegin downloads of GapFiller.....\n")
			download.file("http://sourceforge.net/projects/ochguiextras/files/chimera/GapFiller/gapfiller_unix64.zip/download", "gapfiller.zip", mode="wb")
	        cat("\nThe GapFiller downloaded version is GapFiller_1_0\n")
	        system(paste("unzip gapfiller.zip", sep=""))
	        system("chmod -fR +x gapfiller")
	}else{
		cat("\nThe selected OS is not supported\n")
	}
	unlink("gapfiller.zip")
	setwd(mydir)
	return()
}







