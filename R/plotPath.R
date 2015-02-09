plotPath <- function(my.path, path.db=c("kegg","biocarta","nci","spike","humancyc","panther"), type=c("entrez","symbol")){
	if(path.db=="kegg"){
	    p <- kegg[[which(names(kegg)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
	}else if(path.db=="biocarta"){
	    p <- biocarta[[which(names(biocarta)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
    }else if(path.db=="nci"){
	    p <- nci[[which(names(nci)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
    }else if(path.db=="spike"){
	    p <- spike[[which(names(spike)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
    }else if(path.db=="humancyc"){
	    p <- humancyc[[which(names(humancyc)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
    }else if(path.db=="panther"){
	    p <- panther[[which(names(panther)==my.path)]]
	    p <- convertIdentifiers(p, type)
	    g <- pathwayGraph(p)
	    plot(g)
    }
}


 
 
 
 
 
 
 
 

