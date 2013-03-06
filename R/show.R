#return the size and class of the object
setMethod("show","list", function(object){
	cat("\nAn object of class list including ", length(object), " fusions\nThe data were acquired by ", 
object[[1]]@fusionInfo$fusionTool, " fusions detection tool.")
})
