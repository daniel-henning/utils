###############################################################
####     RNAseq count normalization by DESeq2 or edgeR.    ####
###############################################################
normalizeGeneCounts <- function(counts, TxDb, method) {
  require("GenomicFeatures")
  length <- width(genes(TxDb))/10 ^ 3 #gene length per kb
  names(length) <- substring(rownames(as.data.frame(genes(TxDb))),1,15)
  counts <- counts[intersect(names(length), rownames(counts)),]
  if (method == "CPM") {
    normalized <-
      as.data.frame(apply(counts,2, function(x) (x/sum(x)) * 10 ^ 6))
  } else if (method == "RPKM") {
      rate <- counts / lengths
      normalized <- rate / sum(counts) * 1e6
  } else if (method == "TPM") {
      rate <- counts / length[rownames(counts)]
      normalized <- rate / sum(rate) * 1e6
  } else {
    stop("method must be CPM, RPKM or TPM")
  }
  return(normalized)
}



############################################################################
####      Batch effect remove of Micro-Array or RNAseq expression       ####
############################################################################
BatchRemove <- function(dataMat, batch){
  library("sva")
  # subset to common genes andbatch correct using ComBat
  dataMat <- dataMat[which(apply(dataMat, 1, sum) > 0),]
  mod <- data.frame(Intercept = rep(1, ncol(dataMat)))
  rownames(mod) <- colnames(dataMat)
  whichbatch <- as.factor(batch)
  combatout <- ComBat(as.matrix(dataMat), whichbatch, mod = mod)
  return(combatout)
}
