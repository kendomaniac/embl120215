
plotCoverage <- function(fset, plot.type=c("exons","junctions"), junction.spanning= 20, fusion.only=FALSE, xlab="nts", ylab="Coverage", main="", col.box1="red", col.box2="green", ybox.lim=c(-4,-1)) {
  x <- fset
  ga <- fusionGA(x)
  ga.x <- coverage(ga)[[1]]
  fusion.transcript <- paste(elementMetadata(fusionGRL(x)$gene1)$KnownGene, elementMetadata(fusionGRL(x)$gene2)$KnownGene, sep=":")
  grl <- fusionGRL(x)
  chimera.tmp <- strsplit(fusion.transcript,":")
  if(length(chimera.tmp[[1]]) > 2){
		     cat("\nAt list one of the fusion element is not an annotated gene.\nPlot of fusion regions for this specific situation is not implemented, yet.\n")
	         return()
  }
  chr.sym <- as.list(org.Hs.egSYMBOL)
  g1 <- chimera.tmp[[1]][1]
  eg1 <- names(chr.sym[which(chr.sym == g1)])
  g2 <- chimera.tmp[[1]][2]
  eg2 <- names(chr.sym[which(chr.sym == g2)])	 
  eg.lst <- list(gene_id=eg1)
  eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, filter=eg.lst, columns=c("tx_id", "tx_name"))
  fusion.trs <- findOverlaps(grl[[1]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
  eg.trs.n <- eg.trs.n[fusion.trs]
  tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
  tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)
  
  eg.lst <- list("tx_id" = tmp.tx)
  eg.trs.e <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene, filter=eg.lst, columns=c("tx_id","exon_id","exon_rank"))
  #handling the 5' end of the fusion
  donor.intron <- NA
  exons.tx <- elementMetadata(eg.trs.e)$tx_id
  exons.tx <- as.list(exons.tx)
  exons.rank <- elementMetadata(eg.trs.e)$exon_rank
  exons.idx<- sapply(exons.tx,function(x,y){grep(y, x)},y= tmp.tx)
  rank.e <- NULL
  for(i in 1:length(exons.idx)){
         rank.e[i] <- exons.rank[[i]][exons.idx[i]]
  }
  donor.exon <- findOverlaps(grl[[1]],  eg.trs.e, type = "any", select = "first", ignore.strand = T)
  tr1 <- eg.trs.e
  end(tr1[donor.exon]) <- end(grl[[1]])
  if(unique(as.character(strand(tr1))) == "-"){
     tr1 <- tr1[donor.exon:length(tr1)]
  }else{
	 tr1 <- tr1[1:donor.exon]
  }  
  new.start1 <- start(tr1) - min(start(tr1))
  new.end1 <- end(tr1) - min(start(tr1))
  new.start1[which(new.start1 == 0)] <- 1       
  tr1.new <- GRanges(seqnames =  elementMetadata(tr1)$exon_id, 
	                     ranges = IRanges(start = new.start1, end= new.end1),strand=strand(tr1), 
	                     rank=elementMetadata(tr1)$rank, 
	                     transcript.id=elementMetadata(tr1)$transcript.id,
	                     exon.id=elementMetadata(tr1)$exon.id)
  myranges1 <- end(tr1.new) - start(tr1.new)
  myranges1 <- c(1,myranges1)
  myranges1.tmp <- myranges1
  for(i in 1:length(myranges1.tmp)){
	          tmp.sum <- sum(myranges1.tmp[1:i])
	           myranges1[i] <- tmp.sum
  }

  ##g2
  eg.lst <- list(gene_id=eg2)
  eg.trs.n <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, filter=eg.lst, columns=c("tx_id", "tx_name"))
  fusion.trs <- findOverlaps(grl[[2]],  eg.trs.n, type = "any", select = "first", ignore.strand = T)
  eg.trs.n <- eg.trs.n[fusion.trs]
  tmp.tx <- as.character(elementMetadata(eg.trs.n)$tx_id)
  tmp.name <- as.character(elementMetadata(eg.trs.n)$tx_name)

  eg.lst <- list("tx_id" = tmp.tx)
  eg.trs.e <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene, filter=eg.lst, columns=c("tx_id","exon_id","exon_rank"))
  #handling the 3' end of the fusion
  acceptor.intron <- NA
  exons.tx <- elementMetadata(eg.trs.e)$tx_id
  exons.tx <- as.list(exons.tx)
  exons.rank <- elementMetadata(eg.trs.e)$exon_rank
  exons.idx<- sapply(exons.tx,function(x,y){grep(y, x)},y=tmp.tx)
  rank.e <- NULL
  for(i in 1:length(exons.idx)){
      rank.e[i] <- exons.rank[[i]][exons.idx[i]]
  }
  acceptor.exon <- findOverlaps(grl[[2]],  eg.trs.e, type = "any", select = "first", ignore.strand = T)
  tr2 <- eg.trs.e
  start(tr2[acceptor.exon]) <- start(grl[[2]])
  if(unique(as.character(strand(tr2))) == "-"){
	  #right!
       tr2 <- tr2[1:acceptor.exon]
  }else{
	   tr2 <- tr2[acceptor.exon:length(tr2)]
  }
  new.start2 <- start(tr2) - min(start(tr2))
  new.end2 <- end(tr2) - min(start(tr2))
  new.start2[which(new.start2 == 0)] <- 1       
  tr2.new <- GRanges(seqnames =  elementMetadata(tr2)$exon_id, 
        ranges = IRanges(start = new.start2, end= new.end2),strand=strand(tr2), 
        rank=elementMetadata(tr2)$rank, 
        transcript.id=elementMetadata(tr2)$transcript.id,
        exon.id=elementMetadata(tr2)$exon.id)
  myranges2 <- end(tr2.new) - start(tr2.new)
  myranges2.tmp <- myranges2
  for(i in 1:length(myranges2.tmp)){
        tmp.sum <- sum(myranges2.tmp[1:i])
         myranges2[i] <- tmp.sum
  }
  #qualcosa che non mi torna controllare
  myranges2 <- c((myranges1[length(myranges1)] +1), (myranges1[length(myranges1)] + myranges2))
#coverage
           start.j1 <- NULL
	       end.j1 <- NULL
           all.coverage <- new("Rle",values=0, lengths=length(ga.x))
           start <- 1
		   end <- length(ga.x)
		   if(length(myranges1) > 2){
	          for(i in 2:(length(myranges1)-1)){
	             start.j1 <- myranges1[i] - junction.spanning
	             end.j1 <- myranges1[i] + junction.spanning
	             all.coverage[start.j1:end.j1] <- ga.x[start.j1:end.j1]
	         }
           }else{
	           i <- 2 
	           start.j1 <- myranges1[i] - junction.spanning
	           end.j1 <- myranges1[i] + junction.spanning
	           all.coverage[start.j1:end.j1] <- ga.x[start.j1:end.j1]
            }
	        if(fusion.only){
		       fusion.coverage <- new("Rle",values=0, lengths=length(ga.x))
		       start.f <- myranges1[length(myranges1)] - junction.spanning
		       end.f <- myranges1[length(myranges1)] + junction.spanning
	           fusion.coverage[start.f:end.f] <- ga.x[start.f:end.f]
	        }
		    start.f <- myranges1[length(myranges1)] - junction.spanning
		    end.f <- myranges1[length(myranges1)] + junction.spanning
	        all.coverage[start.f:end.f] <- ga.x[start.f:end.f]
	  	    start.j2 <- NULL
		    end.j2 <- NULL
		    if(length(myranges2) > 2){
	          for(i in 2:(length(myranges2)-1)){
	             start.j2 <- myranges2[i] - junction.spanning
	             end.j2 <- myranges2[i] + junction.spanning
		         all.coverage[start.j2:end.j2] <- ga.x[start.j2:end.j2]
	          }
	        }else{
		         i <- 2
	             start.j2 <- myranges2[i] - junction.spanning
	             end.j2 <- myranges2[i] + junction.spanning
		         all.coverage[start.j2:end.j2] <- ga.x[start.j2:end.j2]	
	        }
	        if(fusion.only){
		            xWindow <- as(window(fusion.coverage, start.f, end.f), "vector")
		            x <- start.f:end.f
		            xlim <- c(start.f, end.f)
		            ylim <- c(0, max(xWindow))
		            plot(x = start.f, y = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, type = "n")
		            polygon(c(start.f, x, end.f), c(0, xWindow, 0), col = "yellow")
		            abline(v=(start.f + junction.spanning), lty=3, col="red")
	       }else{
	            xWindow <- as(window(ga.x, start, end), "vector")
	            x <- start:end
	            xlim <- c(start, end)
	            ylim <- c(0, max(xWindow))
	            jWindow <- as(window(all.coverage, start, end), "vector")
	            plot(x = start, y = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, type = "n")
	            if(plot.type=="exons"){
	                polygon(c(start, x, end), c(0, xWindow, 0), col = "blue")
	            }else if(plot.type=="junctions"){
	                polygon(c(start, seq(from=start, to=end, length=length(jWindow)), end), c(0, jWindow, 0), col = "yellow")
	            }
	            for(i in 1:(length(myranges1)-1)){
	                rect(xleft=myranges1[i], ybottom=ybox.lim[1], xright=myranges1[i+1], ytop=ybox.lim[2], col=col.box1) 
	            }
	            for(i in 1:(length(myranges2)-1)){
	              rect(xleft=myranges2[i], ybottom=ybox.lim[1], xright=myranges2[i+1], ytop=ybox.lim[2], col=col.box2) 
	            }
	      }
}
