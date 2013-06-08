gapfillerRun <- function(fusion.fa, seed1, seed2, gapfiller=NULL, seed.ins=200, seed.var=50, block.length=6, overlap=6, max.length=5000, slack=6, k=6, global.mismatch=5, perc.identity=0.6){
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
    system(paste(gapfiller," --output contigs.fasta --statistics contigs.stats --seed1 ",seed1,
                  " --seed2 ",seed2," --query r1r2.fasta --seed-ins ",seed.ins," --seed-var ",seed.var,
                  " --block-length ",block.length," --overlap  ",overlap," --read-length 200 ext.Thr 1 --max-length ",
                  max.length," --slack ",slack," --k ",k," --global-mismatch ",global.mismatch," --perc-identity ",perc.identity, sep=""))
    gapfiller.stats <- read.table("contigs.stats", sep = " ", header = T)
	gapfiller.stats <- t(gapfiller.stats)
	dimnames(gapfiller.stats)[2] <- "statistics"
	cat("\n")
	print(gapfiller.stats)
	cat("\n")
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
    alignment <- alignment[intersect(tmp.acceptor, tmp.donor)]
    if(length(intersect(tmp.acceptor, tmp.donor) > 0)) {
	   cat("\nde novo alignment has overlap over the fusion break point\n")
	   setwd(mydir)
	   return(alignment)
    }else{
	   	cat("\nde novo alignment has no overlap over the fusion break point\n")
	    setwd(mydir)
		return()
    }    
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







