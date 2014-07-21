subreadRun <- function(ebwt, input1, input2, outfile.prefix="accepted_hits", alignment=c("se","pe"), cores=1){
	buildindex <- NULL
	align <- NULL
    if(.Platform$OS.type!="windows"){	
	     require(Rsubread) || stop("\nMissing Rsubread package\n")
	     chimera.db <- paste("chimeraDB",gsub("[' '| :]","-", date()),sep="_")
	     buildindex(basename=chimera.db,reference=ebwt)
	     if(alignment=="se"){
		     tmp.fq <- paste("tmp",gsub("[' '| :]","-", date()),sep="_")
		     system(paste("cat ",input1," ",input2," >",tmp.fq,".fastq",sep=""))
		     align(index=chimera.db,readfile1=paste(tmp.fq,".fastq", sep=""),output_file=paste(outfile.prefix,".bam",sep=""), nthreads=cores)
		     sortBam(paste(outfile.prefix,".bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
		     indexBam(paste(outfile.prefix,"_sorted.bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
             filterBam(paste(outfile.prefix,"_sorted.bam",sep=""),paste(outfile.prefix,"_mapped.bam",sep=""), param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
	     }else{
		     align(index=chimera.db,readfile1=input1,readfile2=input2,output_file=paste(outfile.prefix,".bam",sep=""), nthreads=cores)
		#     asBam(paste(outfile.prefix,".sam",sep=""), outfile.prefix, overwrite=TRUE)
		     sortBam(paste(outfile.prefix,".bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
		     indexBam(paste(outfile.prefix,"_sorted.bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
             filterBam(paste(outfile.prefix,"_sorted.bam",sep=""),paste(outfile.prefix,"_mapped.bam",sep=""), param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
	    }
   }else{
 	   chimera.db <- paste("chimeraDB",gsub("[' '| :]","-", date()),sep="_")
	   tmp <- bowtie_build(references=ebwt, outdir=getwd(), prefix=chimera.db, force=TRUE)
	   if(alignment=="se"){
	       tmp.fq <- paste("tmp",gsub("[' '| :]","-", date()),sep="_")
	       system(paste("cat ",input1," ",input2," >",tmp.fq,".fastq",sep="")) 
		   bowtie(sequences=tmp.fq, index=chimera.db, type ="single", outfile=paste(outfile.prefix,".bam",sep=""), sam=TRUE, best=TRUE, force=TRUE)
		   asBam(paste(outfile.prefix,".sam",sep=""),outfile.prefix,overwrite=T)
	       sortBam(paste(outfile.prefix,".bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
	       indexBam(paste(outfile.prefix,"_sorted.bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
           filterBam(paste(outfile.prefix,"_sorted.bam",sep=""),paste(outfile.prefix,"_mapped.bam",sep=""), param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
	   }else{
		   bowtie(sequences=list(input1,input2), index=chimera.db, type ="paired", outfile=paste(outfile.prefix,".bam",sep=""), sam=TRUE, best=TRUE, force=TRUE)
   	       asBam(paste(outfile.prefix,".sam",sep=""), outfile.prefix, overwrite=TRUE)
   		   sortBam(paste(outfile.prefix,".bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
   		   indexBam(paste(outfile.prefix,"_sorted.bam",sep=""), paste(outfile.prefix,"_sorted",sep=""))
           filterBam(paste(outfile.prefix,"_sorted.bam",sep=""),paste(outfile.prefix,"_mapped.bam",sep=""), param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))   	
	  }
   }
}











