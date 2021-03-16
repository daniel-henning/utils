####wilcox test for identifing gene expression pattern of sepecific cell type.

group <- (benjamin.tumor.metabolic$label)

exprMat <- benjamin.tumor.metabolic[,-5338]

wilcox.pattern <- function(exprMat, group, pvalue = 0.000001){
  uni.group <- unique(group)
  gene.pattern <- list();for (i in 1:length(uni.group)) { gene.pattern[[`i`]] <- list()}; names(gene.pattern) <- uni.group
  gene.pattern.2 <- list();for (i in 1:length(uni.group)) { gene.pattern.2[[`i`]] <- c("NA")}; names(gene.pattern.2) <- uni.group
  for (i in 1:length(uni.group)) { gene.pattern[[`i`]] <- gene.pattern.2[setdiff(1:length(uni.group),i)] }; rm(gene.pattern.2)
  ####
  for (g in uni.group) {
    for (s in setdiff(uni.group,g)) {
      for (i in 1:ncol(exprMat)) {
        expr <- exprMat[,i]
        w.res <- wilcox.test(x = expr[which(group == g)], y = expr[which(group == s)], paired = FALSE)
        if(!is.na(w.res$p.value) && w.res$p.value < pvalue) { gene.pattern[[`g`]][[`s`]] =  c(gene.pattern[[`g`]][[`s`]], "TRUE")}else{
          gene.pattern[[`g`]][[`s`]] =  c(gene.pattern[[`g`]][[`s`]], "FALSE")}
      }
    }
  }
  gene.pattern.2 <- list();for (i in 1:length(uni.group)) { gene.pattern.2[[`i`]] <- c('NA')}; names(gene.pattern.2) <- uni.group
  for(cell0 in gene.pattern){
    sp.genes <- c(); cell.type <- setdiff(uni.group,names(cell0))
    for(cell in cell0){if(length(sp.genes)==0){sp.genes<-which(cell[-1]=="TRUE")}else{sp.genes<-intersect(sp.genes,which(cell[-1]=="TRUE"))}}
    gene.pattern.2[[`cell.type`]] <- colnames(exprMat)[sp.genes]
  }
  return(list(diff.gene = gene.pattern, uni.gene = gene.pattern.2))
}


gene.pattern <- wilcox.pattern(exprMat = exprMat, group = group)


ann_colors = list(Cell.Type = c(Mal = "#e21a1c", Macrophage = "#fddf8a", Endo. = "#e6ab02", T.CD4 = "#1b9e77", 
                                CAF = "#d95f02", T.CD8 = "#1b9e77", T.cell = "#1b9e77", NK = "#756fb3", B.cell = "#66a61e"))

genes <- unique(c(gene.pattern$uni.gene$T.cell, gene.pattern$uni.gene$T.CD8, gene.pattern$uni.gene$T.CD4))

pheatmap(log2(exprMat[,genes] + 1), annotation_row = annotation_row, 
         show_rownames = FALSE, show_colnames = TRUE, annotation_colors = ann_colors,
         col = colorRampPalette(c('steelblue','snow','firebrick3'))(120))


genes.new <- c()

mal.idx <- which(group %in% c("Mal"))

for (gene in genes) {
  fc <- mean((exprMat[mal.idx,gene]))/mean((exprMat[-mal.idx,gene]))
  if (fc > 4 | fc< 0.25){
    genes.new = c(genes.new,gene)
  }
}


idx <- which(apply(exprMat[,genes],2,function(x)quantile(x,0.5)>1))


