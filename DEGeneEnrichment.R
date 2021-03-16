DEGeneEnrichment <- function(DEgenes = DEgenes,
                             logFC = 1 , 
                             pcutoff = 0.05,
                             qcutoff = 0.05, 
                             organism = "hsa",
                             outdir = "./"){
  ####################################
  ### Turning warnings off for the sake of a cleaner output
  ####################################
  oldw <- getOption("warn")
  options( warn = -1 )

  source("D://Scripts/utils/functions.R")
  suppressMessages(if(!require(ggplot2)){BiocManager::install('ggplot2')})
  suppressMessages(if(!require(pathview)){BiocManager::install('pathview')})
  suppressMessages(if(!require(clusterProfiler)){BiocManager::install('clusterProfiler')})
  
  ####load data####
  if(!is.object(DEgenes)){ genes = readFile(DEgenes, sep = ','); rownames(genes) <- as.factor(genes[,1])}else{genes = DEgenes}

  if (is.null(pcutoff)){
    gene      = genes[which(abs(genes$log2FoldChange) > logFC & genes$padj < 0.05),1]
    gene_up   = genes[which(genes$log2FoldChange >  logFC & genes$padj < 0.05),1]
    gene_down = genes[which(genes$log2FoldChange < -logFC & genes$padj < 0.05),1]
    #genes     = genes[which(genes$log2FoldChange < -1 & genes$padj < 0.05),]
  }else{
    gene      = genes[which(abs(genes$log2FoldChange) >logFC & genes$pvalue < pcutoff),1]
    gene_up   = genes[which(genes$log2FoldChange >  logFC & genes$pvalue < pcutoff),1]
    gene_down = genes[which(genes$log2FoldChange < -logFC & genes$pvalue < pcutoff),1]
    #genes     = genes[which(abs(genes$log2FoldChange) >1 & genes$padj < pcutoff),]
  }

  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(paste(outdir,'/GO',sep = ''))){dir.create(paste(outdir,'/GO',sep = ''))}
  if(!dir.exists(paste(outdir,'/KEGG',sep = ''))){dir.create(paste(outdir,'/KEGG',sep = ''))}
  if(!dir.exists(paste(outdir,'/KEGG/Pathview',sep = ''))){dir.create(paste(outdir,'/KEGG/Pathview',sep = ''))}

  ####convert gene symbol to entrezID####
  eg <- geneIDannotation(gene, name = "symbol")
  gene <- as.character(unique(na.omit(eg$gene_id)))
  ####
  gcSample <- list()
  gcSample$exprUp   <- as.integer(as.character(na.omit(geneIDannotation(gene_up)$gene_id)))
  gcSample$exprDown <- as.integer(as.character(na.omit(geneIDannotation(gene_down)$gene_id)))
  ####
  tryCatch({
    ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG", organism = organism)
    ggsave(paste(outdir, '/KEGG/KEGG_compare_cluster.png', sep = ''),dotplot(ck), width = 9, height = 9, units = "in")
    write.csv(as.data.frame(ck),paste(outdir, '/KEGG/KEGG_compare_cluster.csv',sep = ''),
              quote = FALSE, row.names = FALSE)
    ####
    cg <- compareCluster(geneCluster = gcSample, fun='enrichGO', OrgDb='org.Hs.eg.db')
    ggsave(paste(outdir, '/GO/GO_compare_cluster.png', sep = ''),dotplot(cg), width = 9, height = 9, units = "in")
    write.csv(as.data.frame(cg),paste(outdir, '/GO/GO_compare_cluster.csv',sep = ''), quote = FALSE, row.names = FALSE)
  }, warning = function(war) {
    # warning handler picks up where error was generated
    print(paste("MY_WARNING:  ",war))

  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("MY_ERROR:  ", err))
  }
  )
  
  

  ############geneList contains three features,numeric vector: fold change or other type of numerical variable.
  eg <- geneIDannotation(genes[,1])
  RT.ENS = length(which(genes[1:100,1] %in% eg$ensembl_id))
  RT.SYM = length(which(genes[1:100,1] %in% eg$symbol))
  RT.ENT = length(which(genes[1:100,1] %in% eg$gene_id))
  if (RT.ENS > 20){
    geneList <- genes[eg$ensembl_id,3]
  }else if (RT.SYM > 20){
    geneList <- genes[eg$symbol,3]
  }else if (RT.ENT > 20){
    geneList <- genes[eg$gene_id,3]
  }
  names(geneList) = as.character(eg$gene_id)
  ####sort the geneList by fold-change.
  #geneList = 2^geneList
  geneList <- sort(geneList,decreasing = TRUE)

  #########################################KEGG analysis###########################################
  ekk <- enrichKEGG(gene, organism = organism, keyType  ='kegg', pvalueCutoff = pcutoff, 
                    pAdjustMethod  = "BH", minGSSize = 10, maxGSSize = 500)
  ####
  ekk2 <- setReadable(ekk, 'org.Hs.eg.db', 'ENTREZID')

  if(dim(ekk2)[1] > 1){
    ggsave(paste(outdir, '/KEGG/KEGG_Enrichment_cnetplot.png', sep = ''), plot = cnetplot(ekk2, foldChange=geneList),
           width = 8, height = 8, units = "in")
    ####
    ggsave(paste(outdir, '/KEGG/KEGG_Enrichment_emapplot.png', sep = ''), plot = emapplot(ekk2),
           width = 8, height = 8, units = "in")
    ####
    ggsave(paste(outdir, '/KEGG/KEGG_Enrichment_dot_plot.png', sep = ''), plot = dotplot(ekk2),
           width = 8, height = 8, units = "in")
    ####
    ggsave(paste(outdir, '/KEGG/KEGG_Enrichment_bar_plot.png', sep = ''), plot = barplot(ekk2),
           width = 8, height = 6, units = "in")
    ####
    ggsave(paste(outdir, '/KEGG/KEGG_Enrichment_heat_plot.png', sep = ''),
           plot = heatplot(ekk2, foldChange=geneList),
           width = 15, height = 8, units = "in")
  }

  ####
  write.csv(as.data.frame(ekk),paste(outdir, "/KEGG/KEGG_Enrichment.csv", sep = ''), row.names =F)
  #pathview of all enrichment KEGG pathways.
  library("pathview")
  ek = as.data.frame(ekk)
  if (dim(ek)[1] > 1 & dim(ek)[1] < 50){KGBubble(ek, outdir, KEGG = TRUE, width = 25, height = 20)
  }else{ KGBubble(ek[1:50,], outdir = paste(outdir,'KEGG',sep = '/'), KEGG = TRUE, width = 25, height = 20) }

  drop.path.id <- which(ek$ID %in% c("hsa04723", "hsa05206"))
  if(length(drop.path.id) > 1){hsa_id <- ek$ID[-drop.path.id]}else{hsa_id <- ek$ID}
  for ( id in hsa_id ){
    pathview(gene.data=geneList,pathway.id=id, species="hsa",
             kegg.dir = "D://Scripts/database/pathways/KEGG/", kegg.native = TRUE,
             limit = list(gene = 4, cpd = 4), node.sum = "max")
  }
  file.remove(list.files()[grep('xml',list.files())])
  file.copy(list.files()[grep('pathview',list.files())], paste(outdir, '/KEGG/Pathview',sep = ''))
  file.remove(list.files()[grep('hsa',list.files())])

  ##########################################enrichGO###############################################
  mfgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "MF", pvalueCutoff = pcutoff, readable = TRUE)
  bpgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "BP", pvalueCutoff = pcutoff, readable = TRUE)
  ccgo <- enrichGO(gene = gene, OrgDb="org.Hs.eg.db",ont = "CC", pvalueCutoff = pcutoff, readable = TRUE)
  mf <- as.data.frame(mfgo)
  m <- c('MF')
  mf_t <- rep(m,times=length(mfgo$ID))
  mf$GO_term <- mf_t
  bp <- as.data.frame(bpgo)
  b <- c('BP')
  bp_t <- rep(b,times=length(bpgo$ID))
  bp$GO_term <- bp_t
  cc <- as.data.frame(ccgo)
  c <- c('CC')
  cc_t <- rep(c, times=length(ccgo$ID))
  cc$GO_term <- cc_t
  GO <- rbind(mf,bp,cc)
  GO2 <- rbind(mf[1:10,],bp[1:10,],cc[1:10,])
  ####
  write.csv(GO, paste(outdir, "/GO/GO_Enrichment.csv", sep = ''), row.names=FALSE, quote = FALSE)

  KGBubble(na.omit(GO), paste(outdir, 'GO',sep = '/'), KEGG = FALSE, width = 30, height = 20)

  ####
  if(dim(GO2)[1] > 1){
    ggsave(paste(outdir, '/GO/GO_MF_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(mfgo, showCategory=10)+ ggtitle("Molecular Function of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 10, units = "in")
    ####
    ggsave(paste(outdir, '/GO/GO_BP_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(bpgo, showCategory=10)+ ggtitle("Biological Process of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 10, units = "in")
    ####
    ggsave(paste(outdir, '/GO/GO_CC_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(ccgo, showCategory=10)+ ggtitle("Cellular Component of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 10, units = "in")
    ####
    ggsave(paste(outdir, '/GO/GO_MF_Enrichment_barplot.png', sep = ''),
           plot = barplot(mfgo, showCategory=10)+ ggtitle("Molecular Function of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 6, units = "in")
    ####
    ggsave(paste(outdir, '/GO/GO_BP_Enrichment_barplot.png', sep = ''),
           plot = barplot(bpgo, showCategory=10)+ ggtitle("Biological Process of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 6, units = "in")
    ####
    ggsave(paste(outdir, '/GO/GO_CC_Enrichment_barplot.png', sep = ''),
           plot = barplot(ccgo, showCategory=10)+ ggtitle("Cellular Component of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 6, units = "in")
    ####
    png(paste(outdir, "/GO/GO_MF_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
    cnetplot(mfgo,categorySize="pvalue",foldChange=geneList)
    dev.off()
    ####
    png(paste(outdir, "/GO/GO_BP_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
    cnetplot(bpgo,categorySize="pvalue",foldChange=geneList)
    dev.off()
    ####
    png(paste(outdir, "/GO/GO_CC_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
    cnetplot(ccgo,categorySize="pvalue",foldChange=geneList)
    dev.off()
    png(paste(outdir, "/GO/GO_BP_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
    emapplot(bpgo)
    dev.off()
    ####
    png(paste(outdir, "/GO/GO_MF_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
    emapplot(mfgo)
    dev.off()
    ####
    png(paste(outdir, "/GO/GO_CC_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
    emapplot(ccgo)
    dev.off()
    ####
    png(paste(outdir, "/GO/topGO_MF.png", sep = ''),width = 1200, height = 1200)
    suppressMessages(plotGOgraph(mfgo))
    dev.off()
    ####
    png(paste(outdir, "/GO/topGO_BP.png", sep = ''),width = 1200, height = 1200)
    suppressMessages(plotGOgraph(bpgo))
    dev.off()
    ####
    png(paste(outdir, "/GO/topGO_CC.png", sep = ''),width = 1200, height = 1200)
    suppressMessages(plotGOgraph(ccgo))
    dev.off()
  }
  
  #### save runing image data.
  save.image(paste(outdir,'DEGenes_Enrichment.RData'))
 }
