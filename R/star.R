starInstallation <- function(binDir, os=c("unix","mac")){
	mydir <- getwd()
	if(os=="unix"){
		starDirLocation  <- paste(.path.package("chimera", quiet = FALSE), "/star", sep="")
		tmp.info <- dir(paste(.path.package("chimera", quiet = FALSE)))
		if(length(grep("star", tmp.info))>0){
		unlink(starDirLocation, recursive = T, force = T)
	    }
		dir.create(starDirLocation, showWarnings = TRUE, recursive = FALSE)
		setwd(starDirLocation)
		cat("\nBegin downloads of STAR.....\n")
		download.file("ftp://ftp2.cshl.edu/gingeraslab/tracks/STARrelease/2.3.0/STAR_2.3.0e.Linux_x86_64_static.tgz", "star.tgz", mode="wb")
        cat("\nThe STAR downloaded version is 2.3.0\n")
        system("gzip -d star.tgz")
        system("tar xvf star.tar")        
        system("cp -fR ./STAR_2.3.0e.Linux_x86_64_static/* .")
        system("rm -fR ./STAR_2.3.0e.Linux_x86_64_static/")
	    system(paste("chmod +x ",starDirLocation, "/star", sep=""))
		system(paste("ln -fs ",starDirLocation,"/star ", binDir, "/STAR", sep=""))
		cat(paste("\nSTAR has been installed and soft links created in\n",binDir,"\n",sep=""))
        setwd(mydir)
        return()
	}else if(os=="mac"){
	starDirLocation  <- paste(.path.package("chimera", quiet = FALSE), "/star", sep="")
	tmp.info <- dir(paste(.path.package("chimera", quiet = FALSE)))
	if(length(grep("star", tmp.info))>0){
	unlink(starDirLocation, recursive = T, force = T)
    }
	dir.create(starDirLocation, showWarnings = TRUE, recursive = FALSE)
	setwd(starDirLocation)
	cat("\nBegin downloads of STAR.....\n")
	download.file("ftp://ftp2.cshl.edu/gingeraslab/tracks/STARrelease/2.3.0/STAR_2.3.0e.OSX_x86_64.tgz", "star.tgz", mode="wb")
    cat("\nThe STAR downloaded version is 2.3.0\n")
    system("gzip -d star.tgz")
    system("tar xvf star.tar")        
    system("cp -fR ./STAR_2.3.0e.OSX_x86_64/* .")
    system("rm -fR ./STAR_2.3.0e.OSX_x86_64/")
    system(paste("chmod +x ",starDirLocation, "/star", sep=""))
	system(paste("ln -fs ",starDirLocation,"/star ", binDir, "/STAR", sep=""))
	cat(paste("\nSTAR has been installed and soft links created in\n",binDir,"\n",sep=""))
    setwd(mydir)
    return()
	}else{
		cat("\nYou have selected a os not compatible with STAR!\n")
		setwd(mydir)
		return()
	}
	
}

starRun <- function(input1, input2, cores=1, star= "STAR", samtools="samtools", fa, alignment=c("se","pe"), chimSegmentMin=0, chimJunctionOverhangMin=20)
{
           check.star <- system("STAR 2>&1", intern=T)
#           if(length(grep("Error",check.tophat[[1]])) > 0){
	        if(length(as.character(check.star)) == 0){
	            cat("\nIt seems that STAR is not installed in your system:\n")
	            cat("please install it using the function: starInstallation\nAfter running starInstallation. Close R and the shell\n")
	            return()
           }
           #creating time base name to differentiate the various runs
	       time <- gsub(" ","_",date())
		   time <- gsub(":","-",time)
		   time <- gsub(":","-",time)
           #building chimera db
           chimera.db <- paste("chimeraDB_",time,sep="-")
	       dir.create(paste(getwd(),chimera.db, sep="/"))
	       dir.create(paste(getwd(),paste("output_",time, sep="-"), sep="/"))
	       mydir <- getwd()
	       setwd(paste(getwd(),paste("output_",time, sep="-"), sep="/"))
	       outputdir <- getwd()
           system(paste(star," --runMode genomeGenerate --genomeDir ", paste(mydir,chimera.db, sep="/"), " --genomeFastaFiles ", paste(mydir,fa, sep="/"), " --runThreadN ", cores,sep=""), wait=T)
           if(alignment=="se"){
	          star.run <- paste(star," --genomeDir ", paste(mydir,chimera.db, sep="/")," --readFilesIn ", paste(mydir,input1, sep="/")," --runThreadN ", cores, " --chimSegmentMin ",chimSegmentMin, " --chimJunctionOverhangMin ",chimJunctionOverhangMin, sep="")
              system(star.run, wait=T)
           }else if(alignment=="pe"){
	          star.run <- paste(star," --genomeDir ", paste(mydir,chimera.db, sep="/")," --readFilesIn ", paste(paste(mydir,input1, sep="/"),paste(mydir,input2, sep="/"), sep=" ")," --runThreadN ", cores, " --chimSegmentMin ",chimSegmentMin, " --chimJunctionOverhangMin ",chimJunctionOverhangMin, sep="")
              system(star.run, wait=T)
           }
           system(paste(samtools, " view -Sb ", paste(getwd(), "Aligned.out.sam", sep="/"), " > ",  paste(getwd(),"accepted_hits.bam", sep="/"), sep=""))
           setwd(mydir)
           unlink(paste(getwd(),chimera.db, sep="/"), recursive = T, force = T)
           return(outputdir)
}
