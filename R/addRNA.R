#insert DNAStringSet for fused trasncripts
setMethod("addRNA","fSet",function(x, rna) {
    tmp <- new("fSet", fusionInfo=x@fusionInfo, fusionLoc=x@fusionLoc, fusionRNA=rna, fusionGA=x@fusionGA)
	return(tmp)
})
