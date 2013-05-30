#evaluate if an fSet object is an fSet object
"is.fSet"<- function(x){
	   tmp.summary <- summary(x)
	   if(tmp.summary[[2]]=="fSet") {return(TRUE)} else {return(FALSE)}
}
