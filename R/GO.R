###plotting GO term
plotGO <- function(go.term, go.class=c("BP","MF","CC")){
  i <- go.term
  ReturnVal <- go.class
#defining the GO class
  if(ReturnVal=="BP"){
     tfG <- GOGraph(i, GOBPPARENTS)
     rootNames <- c("all", "top")
     root <- rootNames[rootNames %in% nodes(tfG)]
     if (length(root) > 0)  tfG <- removeNode(root, tfG)
     nL <- nodes(tfG)
     nLterm = unlist(getGOTerm(nL), use.names=FALSE)
     names(nLterm) = nL
     nLterm = gsub("_", " ", nLterm)
     tGfnA <- list()
     tmp <- rep("ellipse", length(nL))
     names(tmp) <- nL
     tGfnA$shape = tmp;
     tmp <- rep(FALSE, length(nL))
     names(tmp) = nL
     tGfnA$fixedsize = tmp
     names(nLterm) <- nL
     tGfnA$label = nLterm
     plot(tfG, nodeAttrs=tGfnA)
  } else if(ReturnVal=="MF"){
     tfG <- GOGraph(i, GOMFPARENTS)
     rootNames <- c("all", "top")
     root <- rootNames[rootNames %in% nodes(tfG)]
     if (length(root) > 0)  tfG <- removeNode(root, tfG)
     nL <- nodes(tfG)
     nLterm = unlist(getGOTerm(nL), use.names=FALSE)
     names(nLterm) = nL
     nLterm = gsub("_", " ", nLterm)
     tGfnA <- list()
     tmp <- rep("ellipse", length(nL))
     names(tmp) <- nL
     tGfnA$shape = tmp;
     tmp <- rep(FALSE, length(nL))
     names(tmp) = nL
     tGfnA$fixedsize = tmp
     names(nLterm) <- nL
     tGfnA$label = nLterm
     plot(tfG, nodeAttrs=tGfnA)
  } else if(ReturnVal=="CC"){
     tfG <- GOGraph(i, GOCCPARENTS)
     rootNames <- c("all", "top")
     root <- rootNames[rootNames %in% nodes(tfG)]
     if (length(root) > 0)  tfG <- removeNode(root, tfG)
     nL <- nodes(tfG)
     nLterm = unlist(getGOTerm(nL), use.names=FALSE)
     names(nLterm) = nL
     nLterm = gsub("_", " ", nLterm)
     tGfnA <- list()
     tmp <- rep("ellipse", length(nL))
     names(tmp) <- nL
     tGfnA$shape = tmp;
     tmp <- rep(FALSE, length(nL))
     names(tmp) = nL
     tGfnA$fixedsize = tmp
     names(nLterm) <- nL
     tGfnA$label = nLterm
     plot(tfG, nodeAttrs=tGfnA)
  } else{
      return()
  } 
}


###extracting DE from GO term

deInGO <- function(go.term, de.universe, org=c("Mus.musculus","Homo.sapiens")){
	if(org=="Mus.musculus"){
		   allGO <- select(Mus.musculus, keys=go.term, columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="GOID")
		}else if(org=="Homo.sapiens"){
		   allGO <- select(Homo.sapiens, keys=go.term, columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="GOID")
		}   
    deGO <- allGO[which(allGO$SYMBOL%in%names(de.universe)[which(de.universe==1)]),]
	return(deGO)
}

