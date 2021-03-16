##################################################################################################
####                         consensusMIBC Subtyping Classifer                                ####
##################################################################################################
MIBCSubtyping <- function(countMat, geneID = "SYMBOL", normalized = FALSE){
  suppressMessages(if(!require(consensusMIBC)){devtools::install_github("cit-bioinfo/consensusMIBC", build_vignettes = TRUE); require(consensusMIBC})
  suppressMessages(if(!require(dplyr)){BiocManager::install("dplyr");require(dplyr)})
  
  if(!normalized){countMat.norm <- log2(countNormlization(round(countMat)) + 1) %>% as.data.frame()}else{countMat.norm = countMat}
  ####
  if(geneID == "ENSEMBL"){
    countMat.norm <- geneLocalAnno(exprMatrix = countMat.norm, gtf = "gencode.gene.info.v22.tsv")
    countMat.norm <- dupRemove(countMat.norm[,-ncol(countMat.norm)])
    geneAnno <- geneIDannotation(rownames(countMat.norm))
    geneAnno <- as.data.frame(unique(na.omit(geneAnno)[,c("gene_id","symbol")]))
    rownames(geneAnno) <- geneAnno$symbol
    countMat.norm$EntrezID <- as.character(geneAnno[rownames(countMat.norm),"gene_id"])
    countMat.norm <- dupRemove(na.omit(countMat.norm), geneID = "EntrezID") %>% as.data.frame()
    res <- getConsensusClass(countMat.norm)
  }
  ####
  if(geneID == "SYMBOL"){ 
    geneAnno <- geneIDannotation(rownames(countMat.norm))
    geneAnno <- as.data.frame(unique(na.omit(geneAnno)[,c("gene_id","symbol")]))
    rownames(geneAnno) <- geneAnno$symbol
    countMat.norm$EntrezID <- as.character(geneAnno[rownames(countMat.norm),"gene_id"])
    countMat.norm <- dupRemove(na.omit(countMat.norm), geneID = "EntrezID") %>% as.data.frame()
    res <- getConsensusClass(countMat.norm)
  }
  ####
  if(geneID == "ENTREZ"){
    res <- getConsensusClass(countMat.norm)
  }
  ####
  return(res)
}
