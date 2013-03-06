#retrieve the chr coordinates of the fusion as GRangesList with two elements
setMethod("fusionGRL","fSet",function(x) return(x@fusionLoc))
