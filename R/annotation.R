#subsidiary function to associate ucsc to EG it uses only a subset of UCSC id, i.e. the first one of the ordered set of transcripts ids associated to a gene.
.ucsc2eg <- function(org=c("hg19","mm9"), ucsc.ids){
	if(org=="hg19"){
	    genes.gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	    simbols.eg <- select(Homo.sapiens,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
	    if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
	         elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
	    }else{
		     cat("\nGENEID in Homo.sapiens are not aligned with GENEID in TxDb.Hsapiens.UCSC.hg19.knownGene\n")
		     return()
	    }
	    #extracting ucsc linked to entrez geneid
	    tmp.x <- select(Homo.sapiens,keys=elementMetadata(genes.gr)$gene_id,columns="UCSCKG",keytype="GENEID")
	}else if(org=="mm9"){
	    genes.gr <- genes(TxDb.Mmusculus.UCSC.mm9.knownGene)
		#trs.gr <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)
	    simbols.eg <- select(Mus.musculus,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
	    if(identical(simbols.eg$GENEID, elementMetadata(genes.gr)$gene_id)){
	         elementMetadata(genes.gr)$symbol <- simbols.eg$SYMBOL
	    }else{
		     cat("\nGENEID in Mus.musculus are not aligned with GENEID in TxDb.Mmusculus.UCSC.mm9.knownGene\n")
		     return()
	    }
	    #extracting ucsc linked to entrez geneid instead of being assocaited to UCSCKG in Mm are associated to TXNAME
	    tmp.x <- select(Mus.musculus,keys=elementMetadata(genes.gr)$gene_id,columns="TXNAME",keytype="GENEID")     	
	}
	#only ucsc associated to entrez geneid  
	tmp.x <- tmp.x[which(tmp.x[,2]%in%ucsc.ids),]
	tmp.x <- tmp.x[which(tmp.x[,1]%in%elementMetadata(genes.gr)$gene_id),]
	tmp.x <- tmp.x[setdiff(seq(1,dim(tmp.x)[1]), which(duplicated(tmp.x[,1]))),]
	tmp.x <- tmp.x[order(tmp.x[,1]),]
	#
	genes.gr <- genes.gr[which(elementMetadata(genes.gr)$gene_id %in% tmp.x[,1])]
	genes.gr <- genes.gr[order(elementMetadata(genes.gr)$gene_id)]
	if(identical(tmp.x[,1], elementMetadata(genes.gr)$gene_id)){
	     elementMetadata(genes.gr)$UCSC <- tmp.x[,2]
	}else{
		cat("\nGENEID in Homo.sapiens are not aligned with UCSC in Homo.sapiens\n")
		return()
	}
	return(genes.gr)
}


annotatingGenes <- function(filename, org=c("hg19","mm9"), truncating.expected.counts=FALSE){
	tmp <- read.table(filename, sep="\t", header=T)
	tmp.ucsc <- as.character(tmp[,2])
	sorting <- function(x){
		x <- strsplit(x, ",")
		x <- x[[1]]
		x <- x[order(x)]
		return(x[1])
	}
	tmp.ucsc <- sapply(tmp.ucsc, sorting)
	rownames(tmp) <- tmp.ucsc
	tmp.ucsc <- sapply(tmp.ucsc, function(x)x[1])
	annotation <- .ucsc2eg(org=org, ucsc.ids=as.character(tmp.ucsc))
	annotation <- annotation[which(elementMetadata(annotation)$UCSC %in% rownames(tmp))]
    tmp <- tmp[which(rownames(tmp) %in% elementMetadata(annotation)$UCSC),]
	tmp <- tmp[order(rownames(tmp)),]
	annotation <- annotation[order(elementMetadata(annotation)$UCSC)]
	if(identical(rownames(tmp), elementMetadata(annotation)$UCSC)){
	     tmp <- cbind(elementMetadata(annotation)$gene_id, elementMetadata(annotation)$symbol, tmp)
		 names(tmp) <- c("Entrez_GeneID", "Symbol", "USCS_gene_cluster","collapsed_transcripts", "length", "effective_length", "expected_count", "TPM", "FPKM")
		 if(truncating.expected.counts){
		      tmp$expected_count <- as.integer(trunc(tmp$expected_count))
	     }
		 return(as.data.frame(tmp))
	}else{
		cat("\nUCSC in count file are not aligned with UCSC in annotation\n")
		return()
	}
}
#dir <- dir()
#dir <- dir[grep(".genes.results", dir)]
#tmp <- annotatingGenes(filename=dir[1], org="mm9")