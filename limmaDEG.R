#############################################################
####         DEGenes expression analysis by limma        ####    
#############################################################
limmaDEG <- function(normExpr, sample.info, type = "array"){
  library(limma)
  library(edgeR)
  if (type == "array"){
    group_list <- factor(sample.info$Group)
    design <- model.matrix(~group_list)
    colnames(design) <- levels(group_list)
    rownames(design) <- rownames(sample.info)
    ####
    fit <- lmFit(normExpr, design)
    fit <- eBayes(fit, trend=TRUE)
    output <- topTable(fit, coef=2,n=Inf)
  }else if(type == "RNAseq"){
    ####
    counts <- round(normExpr)
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    ####
    group_list <- factor(sample.info$Group)
    design <- model.matrix(~group_list)
    colnames(design) <- levels(group_list)
    rownames(design) <- rownames(sample.info)
    ####
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    output <- topTable(fit, coef=2,n=Inf)
  }
  return(output)
}
