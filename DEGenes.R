################
DEGenes <- function(exprMatrix, sample.info, designMatrix = NULL, outdir = "./",
                    SampleName = "SampleName", Group = "Group", package = "DESeq2"){
  ####################################
  ### Turning warnings off for the sake of a cleaner aoutput
  ####################################
  oldw <- getOption("warn")
  options( warn = -1 )
  ####
  source("D://Scripts/utils/functions.R")
  suppressMessages(if(!require(easypackages)){BiocManager::install('easypackages')})
  packages(c("DESeq2", "ggplot2", "org.Hs.eg.db", "edgeR", "limma", "ggpubr",
             "pheatmap", "RColorBrewer", "EnhancedVolcano", "BiocParallel"))
  #starting parallel
  register(MulticoreParam(8))

  ##load sample RNA-seq count data matrix
  if(!is.object(exprMatrix)){rawdata <- loadExprData(exprMatrix)
  }else{rawdata = exprMatrix;cat("the input exprMatrix being checked...\n")}

  ##load sample information data matrix
  if(!is.object(sample.info)){sample.info <- loadSampleInfo(sample.info)}

  rownames(sample.info) <- sample.info[,SampleName]

  # introduce batch to DEGenes analysis.
  if("batch" %in% colnames(sample.info)){ batch <- sample.info$batch }else{ batch <- NA }

  ##load Homo_sapiens.gene.GRCh38.93.gtf as gene transform model.
  gtf = 'D://Database/HomoGeneInfo/gencode.gene.info.v22.tsv'

  ## create output directory if it isn't exist.
  if(!dir.exists(outdir)){dir.create(outdir)}else{cat("outpur directory is already exist...\n")}

  ## make the design matrix for next step analysis.
  group <- as.factor(sample.info[,Group])
  names(group) <- rownames(sample.info)
  factors = levels(group)

  ASY = toTable(org.Hs.egSYMBOL2EG)
  gl  = geneIDannotation(ASY$gene_id)
  #### 判断是否有比较矩阵。
  if(!is.object(designMatrix)){ designMatrix = read.table(designMatrix, header = TRUE)}

  exprSet <- as.matrix(round(rawdata)[,names(group)])

  exprSet <- exprSet[-which(apply(exprSet,1,sum) == 0 ),]

  colData <- data.frame(row.names = colnames(exprSet), group)
  if(!is.na(batch)){colData$batch = batch}
  ## choose DEGene analysis package (DESeq2 or edgeR).
  if (package == 'DESeq2'){
    if(!dir.exists(paste(outdir, '/DESeq2', sep = ''))){dir.create(paste(outdir, '/DESeq2', sep = ''))}
    ## create dds object from raw count matrix.
    if(!is.na(batch)){
      dds  <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, design = ~ batch + group)
      dds2 <- DESeq(dds)
    }else{
      dds  <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, design = ~ group)
      dds2 <- DESeq(dds)
    }

    ####

    if(ncol(dds2) > 50){
      vsd  <- vst(dds, blind = FALSE)
      if(is.object(batch)){assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)}
      rld.data <- as.data.frame(vsd@assays@.xData$data[[1]], row.names = rownames(dds2))
    }else{
      rld  <- rlog(dds, blind=FALSE)
      if(is.object(batch)){assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)}
      rld.data <- as.data.frame(rld@assays@.xData$data[[1]], row.names = rownames(dds2))
    }
    #resultsNames(dds2)

    ################################################################################################
    df = data.table::melt(rld.data)
    p <- ggplot(data = df, aes(x = variable, y = value)) +
      geom_boxplot(aes(fill=variable)) + theme(legend.position='none') +
      ylab("rlog normalized gene expression") + xlab("samples")
    ggsave(plot = p, filename = paste(outdir, '/All genes expression.pdf',sep = ''),
           device = "pdf",units = "in", width = 12, height = 8)
    rm(df)
    #### pca plot.
    ggsave(plot = plotPCA(rld, c("batch","group"))+ ggtitle("PCA of input RNAseq expression data"),
           filename =  paste(outdir, '/PCA of Sample\'s relationship.pdf',sep = ''),
           device = "pdf", units = "in", width = 7, height = 7)
    #### Heatmap of the sample-to-sample distances
    sampleDists <- dist(t(rld.data))
    sampleDistMatrix <- as.matrix(sampleDists)
    colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
    anno_col = data.frame(Group = sample.info[rownames(sampleDistMatrix),Group], row.names = rownames(sampleDistMatrix))
    pdf(paste(outdir, '/heatmap of Sample\'s correlation.pdf',sep = ''), width = 12, height = 12)
    pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,col=colors,
             annotation_col = anno_col, annotation_legend = TRUE,
             border_color="grey", cellwidth = 25, cellheight = 25)
    dev.off()
    ################################################################################################
    ##
    if(is.null(designMatrix)){
      for (i in 1:l ) { for (j in 1:l) {if(j > i){
        res <- results(dds2, contrast=c("group",factors[i],factors[j]))
        res <- res[which(res$baseMean >= 4),]
        #resOrdered <- res[order(res$pvalue),]
        res =  as.data.frame(res)
        selected <- rownames(res)[which(abs(res$log2FoldChange) > 1 & res$padj < 0.05)]
        ####
        res <- res[order(res$padj),]

        res <- meanCount(DE_Results = res,i = i,j = j, colData = colData,factors = factors, exprSet = exprSet)
        #### if input geneID is ensembl_id, then transform it to symbol.
        if (length(which(rownames(res)[1:100] %in% gl$ensembl_id)) > 1){
          res <- geneLocalAnno(gtf = gtf, exprMatrix = res)
        }
        #### write the final DESeq2 results to file.
        write.csv(res,paste(outdir, '/DESeq2/',factors[i],'_vs_',factors[j],'.csv',sep = ''), quote = FALSE)

        ####
        h.idx <- names(group[which(group %in% c(factors[i],factors[j]))])

        #### heatmap of differential genes.
        ggsave(plot = pheatmap(rld.data[selected,h.idx], scale = 'row',
                               annotation_col = colData, clustering_method = "complete",
                               show_colnames = T, show_rownames = F, cluster_cols = T, cluster_rows = T,
                               col = colorRampPalette(c('steelblue','snow','firebrick3'))(120), border_color = NA),
               filename = paste(outdir, '/DESeq2/', factors[i],'_vs_', factors[j],'.DEGenes.heatmap.png',sep = ''),
               units = "in", width = 8, height = 8)

        #### volcano plot.
        if (length(which(rownames(res)[1:100] %in% gl$ensembl_id)) > 1){ lab = res$geneSymbol }else{ lab = rownames(res)}
        p <- EnhancedVolcano(res, lab = lab,
                             title = "Volcano plot of DESeq2 analysis results", subtitle = '',
                             caption = "FC cutoff, 2; p-value cutoff, 10e-2",
                             x = "log2FoldChange", y = "padj",  legendVisible = FALSE,
                             pCutoff = 10e-2, FCcutoff = 2, ylim = c(0,15), xlim = c(-12,12))
        ggsave(plot = p, filename = paste(outdir, '/DESeq2/', factors[i],'_vs_', factors[j],'.DEGenes.volcano.png',sep = ''),
               device = "png", units = 'in', width = 8, height = 8)
      }
      }
      }
    }else{
      for (i in 1:dim(designMatrix)[1]){
        res <- results(dds2, contrast=c("group",as.character(designMatrix[i,1]),as.character(designMatrix[i,2])))
        res <- res[which(res$baseMean >= 4),]
        #resOrdered <- res[order(res$pvalue),]
        res =  as.data.frame(res)
        selected <- rownames(res)[which(abs(res$log2FoldChange) > 1 & res$padj < 0.05)]
        ####
        res <- res[order(res$padj),]

        res <- meanCount(DE_Results = res, colData = colData,
                         factors = as.character(t(designMatrix[i,])), exprSet = exprSet)
        #### if input geneID is ensembl_id, then transform it to symbol.
        if (length(which(rownames(res)[1:100] %in% gl$ensembl_id)) > 1){
          res <- geneLocalAnno(gtf = gtf, exprMatrix = res)}
        #### write the final DESeq2 results to file.
        write.csv(res,paste(outdir, '/DESeq2/',designMatrix[i,1],'_vs_',designMatrix[i,2],'.csv',sep = ''), quote = FALSE)

        ####
        h.idx <- names(group[which(group %in% c(as.character(designMatrix[i,1]),as.character(designMatrix[i,2])))])

        #### heatmap of differential genes.
        ggsave(plot = pheatmap(rld.data[selected,h.idx], scale = 'row',
                               annotation_col = colData, clustering_method = "complete",
                               show_colnames = T, show_rownames = F, cluster_cols = T, cluster_rows = T,
                               col = colorRampPalette(c('steelblue','snow','firebrick3'))(120), border_color = NA), filename = paste(outdir, '/DESeq2/', designMatrix[i,1],'_vs_',designMatrix[i,2],'.DEGenes.heatmap.png',sep = ''),
               units = "in", width = 8, height = 8)

        #### volcano plot.
        if (length(which(rownames(res)[1:100] %in% gl$ensembl_id)) > 1){ lab = res$geneSymbol }else{ lab = rownames(res)}
        p <- EnhancedVolcano(res, lab = lab,
                             title = "Volcano plot of DESeq2 analysis results", subtitle = '',
                             caption = "FC cutoff, 2; p-value cutoff, 10e-2",
                             x = "log2FoldChange", y = "padj",  legendVisible = FALSE,
                             pCutoff = 10e-2, FCcutoff = 2, ylim = c(0,15), xlim = c(-12,12))
        ggsave(plot = p, filename = paste(outdir, '/DESeq2/', designMatrix[i,1],'_vs_',designMatrix[i,2],'.DEGenes.volcano.png',sep = ''),
               device = "png", units = 'in', width = 8, height = 8)
      }
    }
    save(list = c("exprMatrix", "dds", "dds2", "designMatrix", "rld.data"),
         file = paste(outdir,'DEGene_Analysis_DESeq2.RData', sep = "/"))

    # } else if (package == 'edgeR'){
    #   if(!dir.exists(paste(outdir, '/edgeR', sep = ''))){dir.create(paste(outdir, '/edgeR', sep = ''))}
    #
    #   if (!is.null(numOfRep) && numOfRep < 2){
    #     if(!dir.exists(paste(outdir, '/edgeR', sep = ''))){dir.create(paste(outdir, '/edgeR', sep = ''))}
    #     ########################edgeR########################
    #     y <- DGEList(counts=rawdata,genes=rownames(rawdata),group = group)
    #
    #     y <- calcNormFactors(y)
    #
    #     bcv <- 0.4
    #
    #     for ( i in 1:l ) { for (j in 1:l) {if( j > i){
    #       et <- exactTest(y, dispersion=bcv^2, pair = c(i,j))
    #       topTags(et)
    #       summary(de <- decideTestsDGE(et))
    #       ####????????
    #       DE <- topTags(et, length(et$genes[,1]))$table
    #
    #       DE <- meanCount(DE_Results = DE,i = i,j = j, colData = colData,factors = factors)
    #
    #       write.csv(DE[,-3],paste(outdir, '/edgeR/',factors[j],'_vs_',factors[i],'.csv',sep = ''), quote = FALSE,row.names = FALSE)
    #     } } }
    #   } else {
    #     if(!dir.exists(paste(outdir, '/edgeR', sep = ''))){dir.create(paste(outdir, '/edgeR', sep = ''))}
    #
    #     design <- model.matrix(~0+group)
    #
    #     y<-DGEList(exprSet,genes=rownames(rawdata),group = group)
    #
    #     keep <- rowSums(cpm(y)>1) >= 2
    #
    #     #an exon is only retained if it is expressed at a count-per-million (CPM) above 1 in at least 2 samples.
    #     y <-y[keep, ]
    #     y <-DGEList(counts=y$counts,genes=y$genes,group = group)
    #     y <-calcNormFactors(y)
    #     y$samples
    #     ## outlier and relationship
    #     y <-plotMDS(y)
    #     y <- estimateDisp(y,design = design)
    #     y <- estimateCommonDisp(y)
    #     y <- estimateTagwiseDisp(y)
    #     ## dispersion
    #     y<-estimateGLMCommonDisp(y,method="deviance", robust=TRUE, subset=NULL)
    #     y<-estimateGLMTrendedDisp(y)
    #     y<-estimateGLMTagwiseDisp(y)
    #     ## tagwise dispersions
    #     ##Testing for DE genes??quasi-likelihood F-tests,
    #     #fit <- glmQLFit(y,design,robust=TRUE)
    #     #qlf <- glmQLFTest(fit)
    #     #topTags(qlf)
    #
    #     ## Testing for DE genes likelihood ratio tests
    #     fit<-glmFit(y, design)
    #     lrt<-glmLRT(fit)
    #     topTags(lrt)
    #     for ( i in 1:l ) { for (j in 1:l) {if( j > i){
    #       et <- exactTest(y, pair = c(i,j))
    #       summary(de <- decideTestsDGE(et))
    #       ####????????
    #       DE <- topTags(et, length(et$genes[,1]))$table
    #       DE <- meanCount(DE_Results = DE,i = i,j = j, colData = colData,factors = factors)
    #       DE <- geneLocalAnno(gtf, DE)
    #       write.csv(DE[,-3],paste(outdir, '/edgeR/',factors[j],'_vs_',factors[i],'.csv',sep = ''), quote = FALSE,row.names = FALSE)
    #     } } }
    #   }
    # }
  }
}


