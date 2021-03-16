############################################
############################################
#### date : 8/12/2020
#### author : hening
#### Caclulate EMT score of tumor samples

EMT.Score.Calc <- function(exp.dat){
  
  # Computes the sample's EMT(Epithelial–mesenchymal transition) score.
  #
  # Args:
  #   exp.dat: data.frame of genes' expression from RNAseq or Microarray, samples in column and gene in row.
  #
  # Returns:
  #   The sample EMTscore in vector.
  
  # function: normalized to standard deviations from the median.
  median.norm <- function(x) (x - median(x)) / sd(x)
  
  # EMT signature genes from Gibbons et al's article.
  Gibbons.mesenchymal.genes <- c("VIM", "CDH2", "FOXC2", "SNAI1", "SNAI2", "TWIST1", "FN1", "ITGB6", "MMP2", "MMP3", "MMP9", "SOX10", "GCLC")
  
  Gibbons.epithelial.genes  <- c("CDH1", "DSP", "OCLN")
  

  Byers.mesenchymal.genes   <- c("ADAM12", "ADAMTS12", "ADAMTS2", "AEBP1", "ANGPTL2", "ANTXR1", "AXL", "BNC2", "CALD1", "CDH2", "CMTM3", "CNRIP1",
                                 "COL10A1", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL8A1", "DACT1", 
                                 "EMP3", "FAP", "FBN1", "FN1", "FSTL1", "GPC6", "GYPC", "HTRA1", "INHBA", "ITGA11", "LOXL2", "LRRC15", "MMP2", 
                                 "MSRB3", "NAP1L3", "NID2", "OLFML2B", "PCOLCE", "PDGFRB", "PMP22", "POSTN", "SPARC", "SPOCK1", "SULF1", "SYT11",
                                 "THBS2", "VCAN", "VIM", "ZEB2")
  
  Byers.epithelial.genes    <- c("AP1G1", "ATP8B1", "CDH1", "CDS1", "CGN", "CLDN4", "CNOT1", "CTNND1", "DYNC1LI2", "ERBB3", "ESRP1", "ESRP2", 
                                 "F11R", "GALNT3", "GRHL2", "HOOK1", "IRF6", "MAP7", "MARVELD2", "MARVELD3", "MYO5B", "OCLN", "PRSS8", "SPINT1")
  
  # Error handling
  if (!all(Gibbons.mesenchymal.genes %in% rownames(exp.dat)) || !all(Gibbons.epithelial.genes %in% rownames(exp.dat))) {
    genes <- c(setdiff(Gibbons.mesenchymal.genes, rownames(exp.dat)), setdiff(Gibbons.mesenchymal.genes, rownames(exp.dat)))
    stop("Some EMT genes of Gibbons do not exist : ", genes, " .")
  }
  
  if (!all(Byers.mesenchymal.genes %in% rownames(exp.dat)) || !all(Byers.epithelial.genes %in% rownames(exp.dat))) {
    genes <- c(setdiff(Byers.mesenchymal.genes, rownames(exp.dat)), setdiff(Byers.epithelial.genes, rownames(exp.dat)))
    stop("Some EMT genes of Byers do not exist : ", genes, " .")
  }
  
  ##############################################################
  #### Gibbons et al's EMT score calculate formula.
  # By a simple addition of the mesenchymal-associated gene values and 
  #  a corresponding subtraction of the epithelial-associated gene values,
  #  each tumor profile would be given a summary score for EMT phenotype. 
  
  EMT.score.Gibbons <- apply(apply(exp.dat[Gibbons.mesenchymal.genes,], 1, median.norm ), 1, sum) - 
    apply(apply(exp.dat[Gibbons.epithelial.genes,], 1, median.norm ), 1, sum)
  
  ##############################################################
  #### Byers et al's EMT score calculate formula.
  # EMT score is calculated as the mean expression of the M markers subtracted by the mean expression of the E markers, 
  #  a higher EMT score means a more mesenchymal phenotype. 
  EMT.score.Byers <- apply(apply(exp.dat[Byers.mesenchymal.genes,], 1, median.norm ), 1, sum)/length(Byers.mesenchymal.genes) - 
    apply(apply(exp.dat[Byers.epithelial.genes,], 1, median.norm ), 1, sum)/length(Byers.epithelial.genes)
  
  ##############################################################
  #### Combining two EMT score.
  EMT.score <- data.frame(Gibbons.EMT.score = EMT.score.Gibbons[colnames(exp.dat)],
                          Byers.EMT.score = EMT.score.Byers[colnames(exp.dat)],
                          row.names = colnames(exp.dat))
  return(EMT.score)
}
