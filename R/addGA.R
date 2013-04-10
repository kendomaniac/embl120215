#insert DNAStringSet for fused trasncripts
setMethod("addGA","fSet",function(x, bam) {
	sortBam(bam, paste(bam,"_sorted", sep=""))
	indexBam(paste(bam,"_sorted.bam", sep=""))
	ga <- readGAlignmentsFromBam(paste(bam,"_sorted.bam", sep=""))
    tmp <- new("fSet", fusionInfo=x@fusionInfo, fusionLoc=x@fusionLoc, fusionRNA=x@fusionRNA, fusionGA=ga)
	return(tmp)
})
