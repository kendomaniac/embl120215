installing.rcasc <-
function(){
  test <- system("docker -v", intern = TRUE)
  if (length(test)==0){
	    cat("\nERROR: Docker seems not to be installed in your system\n")
	    return(FALSE)
	  }else{
	    cat(paste("\n In your system the following version of Docker is installed:\n",test,sep=""))
	    return(TRUE)
	  }
  
  require(devtools)
	install_github("kendomaniac/rCASC", ref="master")
}
