## necessary generics
if(is.null(getGeneric("fusionData"))) setGeneric("fusionData",function(x) standardGeneric("fusionData"))
if(is.null(getGeneric("fusionGRL"))) setGeneric("fusionGRL",function(x) standardGeneric("fusionGRL"))
if(is.null(getGeneric("fusionRNA"))) setGeneric("fusionRNA",function(x) standardGeneric("fusionRNA"))
if(is.null(getGeneric("addRNA"))) setGeneric("addRNA",function(x, rna) standardGeneric("addRNA"))
if(is.null(getGeneric("fusionGA"))) setGeneric("fusionGA",function(x) standardGeneric("fusionGA"))
if(is.null(getGeneric("addGA"))) setGeneric("addGA",function(x, bam) standardGeneric("addGA"))
