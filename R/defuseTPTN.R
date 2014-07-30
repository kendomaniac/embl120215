defuseTPTN <- function(){
	chimeraDirLocation  <- path.package("chimera", quiet = FALSE)
	cat("\nReading supplementary table 8 of deFuse paper McPherson et al. (2011)PLoS Comput Biol 7(5): e1001138\n") 
	tmp <- read.table(paste(chimeraDirLocation, "/examples/defuse_TP_FP.txt", sep=""), header=TRUE)
	cat(paste("creating a list of ", dim(tmp)[1], " fSet Objects\n"))
	
	.my.newfset <- function(x){
		gr1 <-
		  GRanges(seqnames = as.character(unlist(x[3])), ranges = IRanges(start = (as.numeric(unlist(x[5])) - 30), end= as.numeric(unlist(x[5]))),
		          strand = as.character(unlist(x[4])),
	              KnownGene = as.character(unlist(x[2])),
	              KnownTranscript =  "",
	              KnownExonNumber = "",
	              KnownTranscriptStrand = "",
	              FusionJunctionSequence =  ""
				  
	    )
		gr2 <-
	     GRanges(seqnames = as.character(unlist(x[7])), ranges = IRanges(start = (as.numeric(unlist(x[9])) - 30), end= as.numeric(unlist(x[9]))),
	          strand = as.character(unlist(x[8])),
              KnownGene = as.character(unlist(x[6])),
              KnownTranscript =  "",
              KnownExonNumber = "",
              KnownTranscriptStrand = "",
              FusionJunctionSequence =  ""
			  
        )
		grl <- GRangesList("gene1" = gr1, "gene2" = gr2)
		
		newfSet(fusionInfo = list(fusionTool = "deFuse",
                                       UniqueCuttingPositionCount = 0,
                                       SeedCount = as.numeric(unlist(x[10])), 
                                       RescuedCount = as.numeric(unlist(x[10]))+ as.numeric(unlist(x[11])), 
                                       SplicePattern = NA, 
                                       FusionGene = paste(as.character(unlist(x[2])), as.character(unlist(x[6])), sep="->"),
                                       frameShift = NA), 
                             fusionLoc = grl,
 		                     fusionRNA = DNAStringSet(),
 		                     fusionGA = GAlignments())
	}
	fset.lst <- list()
	for(i in 1:dim(tmp)[1]){
		 fset.lst[[i]]  <- .my.newfset(tmp[i,])
	}
	return(fset.lst)
}

