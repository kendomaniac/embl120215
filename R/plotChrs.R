#plotting genes in chjromosomes
library(lattice)
  .prepanel.idio <-
      function(chrs, ..., xlim)
  {
      if (missing(xlim))
          list(xlim=c(1, max(chrs[["Length"]])))
      else
          list(xlim=xlim)
  }

  .panel.idio <-
      function(x, y, chrs, ..., gtick=.2, colp="blue", colm="red")
  {
      ## helpers: y coordinates at integer values
      yy <- unique(y)
      o <- order(yy)
      yi <- match(y, yy[o])
      ## draw 'chromosomes'
      with(chrs, {
          idx <- match(yy, Chromosome)
          panel.segments(1, seq_along(yy), Length[idx][o], seq_along(yy))
      })
      ## plot genes
      panel.segments(abs(x), yi, abs(x),
                     ifelse(x<0, yi-gtick, yi+gtick),
                     ..., col=ifelse(x<0, colm, colp))
}

  .idioplot <-
      function(chrs, ..., prepanel=.prepanel.idio, panel=.panel.idio)
  {
      xyplot(..., prepanel=prepanel, panel=panel, chrs=chrs)
  }


  #de.df contains log2FC, SYMBOL, ENTREZID, GENENAME
  plotChrs <- function(de.df, org, genome, plot=c("bars","chart")){
    	if(org=="Mus.musculus"){
    		if(genome=="mm10"){
  			    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  			    genes.m <- select(Mus.musculus, keys=de.df$ENTREZID, columns=c("ENTREZID", "CHRLOC"), keytype="ENTREZID")
  		    }else if(genome=="mm9"){
  				txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  				genes.m <- select(Mus.musculus, keys=de.df$ENTREZID, columns=c("ENTREZID", "CHRLOC"), keytype="ENTREZID")
  		    }
		}else if(org=="Homo.sapiens"){
    		if(genome=="hg19"){
  			   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  			   genes.m <- select(Homo.sapiens, keys=de.df$ENTREZID, columns=c("ENTREZID", "CHRLOC"), keytype="ENTREZID")
  		    }else if(genome=="hg38"){
  				txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  				genes.m <- select(Homo.sapiens, keys=de.df$ENTREZID, columns=c("ENTREZID", "CHRLOC"), keytype="ENTREZID")
   		    }
    	}
		
  	genome.gr <- genes(txdb)
  	de.gr <-  genome.gr[which(names(genome.gr)%in%de.df$ENTREZID)]
  	#adding information to GenomicRange object
  	de.gr <- de.gr[order(elementMetadata(de.gr)[,1])]
  	de.df <- de.df[order(de.df$ENTREZID),]		
  	elementMetadata(de.gr)[,2] <- de.df$SYMBOL
  	names(elementMetadata(de.gr))[2] <- "symbol"
  	#extracting the chrs
  	chrs <- as.character(seqnames(de.gr))
  	chr.length <- seqlengths(de.gr)
  	chrs <- chrs[order(chrs)]
  	chrs <- rle(chrs)
  	chr.length <- chr.length[which(names(chr.length)%in%chrs$values)]
  	chr.length <- chr.length[order(names(chr.length))]
  	chrs.v <- chrs$lengths
  	names(chrs.v) <- chrs$values
  	chrs.v <- chrs.v[order(names(chrs.v))]
  	chrs.v <- chrs.v[order(chr.length, decreasing=TRUE)]
	
	genes.m <- data.frame(Id=genes.m$ENTREZID,
	                      Chromosome=factor(paste("chr",genes.m$CHRLOCCHR, sep=""), levels=names(chr.length)),
						  Position=as.integer(genes.m$CHRLOC))
	
	
  	##plotting genes on chrs

  	chrs.m <- data.frame(Chromosome=factor(names(chr.length),levels=names(chr.length)),
  						Length=as.integer(chr.length))
		
	
  	if(plot=="bars"){
  		#plotting DE in chrs in funciton of chr length
  	    barplot(chrs.v, las=2)
  	}else if(plot=="chart"){
  	    idio <- .idioplot(chrs.m, Chromosome~Position, genes.m)
		show(idio)
	}else{
			cat("\nThis is not a plot parameter\n")
	}
  	return(list(gr=de.gr, chr.freq=chrs.v))
  }



