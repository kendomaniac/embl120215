subreadRun <- function(ebwt, input1, input2, outfile.prefix="accepted_hits", alignment=c("se","pe"), cores=1){
    require(Rsubread) || stop("\nMissing Rsubread package\n")
	chimera.db <- paste("chimeraDB",gsub("[' '| :]","-", date()),sep="_")
	buildindex(basename=chimera.db,reference=ebwt)
	if(alignment=="se"){
		tmp.fq <- paste("tmp",gsub("[' '| :]","-", date()),sep="_")
		system(paste("cat ",input1," ",input2," >",tmp.fq,".fastq",sep=""))
		align(index=chimera.db,readfile1=paste(tmp.fq,".fastq", sep=""),output_file=paste(outfile.prefix,".sam",sep=""), nthreads=cores)
#        filterSamReads(input=paste(outfile.prefix,".sam",sep=""), output=paste(outfile.prefix,"_mapped.sam",sep=""), filter="includeAligned")
		asBam(paste(outfile.prefix,".sam",sep=""), outfile.prefix, overwrite=TRUE)
	}else{
		align(index=chimera.db,readfile1=input1,readfile2=input2,output_file=paste(outfile.prefix,".sam",sep=""), nthreads=cores)
 #       filterSamReads(input=paste(outfile.prefix,".sam",sep=""), output=paste(outfile.prefix,"_mapped.sam",sep=""), filter="includeAligned")
		asBam(paste(outfile.prefix,".sam",sep=""), outfile.prefix, overwrite=TRUE)
	}
}











