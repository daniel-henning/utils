genelist.enrichment <- function(gene.list = gene.list, pcutoff = 0.05, qcutoff = 0.05, outdir = "./"){
  ####载入所需安装包
  suppressWarnings(if(!require(ggplot2)){BiocManager::install('ggplot2')})
  suppressWarnings(if(!require(pathview)){BiocManager::install('pathview')})
  suppressWarnings(if(!require(data.table)){BiocManager::install('data.table')})
  suppressWarnings(if(!require(org.Hs.eg.db)){BiocManager::install('org.Hs.eg.db')})
  suppressWarnings(if(!require(clusterProfiler)){BiocManager::install('clusterProfiler')})
  suppressWarnings(source("./utils/functions.R"))
  ####symbol转entrezID
  idx <- grep('///|;',gene.list)
  if (length(idx) >= 1) { 
    sgene <- as.vector(na.omit(unlist(tstrsplit(gene.list[idx], ";|///")))) 
    gene.list <- c(gene.list[-idx],sgene)
  }
  ####creating output results directoty.
  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(paste(outdir,'/GO',sep = ''))){dir.create(paste(outdir,'/GO',sep = ''))}
  if(!dir.exists(paste(outdir,'/KEGG',sep = ''))){dir.create(paste(outdir,'/KEGG',sep = ''))}
  if(!dir.exists(paste(outdir,'/KEGG/Pathview',sep = ''))){dir.create(paste(outdir,'/KEGG/Pathview',sep = ''))}
  
  ####convert gene symbol to entrezID####
  eg   <- geneIDannotation(as.character(gene.list), name = T, org = "hsa")
  gene.list <- as.vector(na.omit(unique(eg$gene_id)))
  
  ####create a false genelist to view these kegg-pathway.
  geneList        <- c(rep(2, length(gene.list)))
  names(geneList) <- gene.list
  #########################################KEGG analysis###########################################
  ekk <- enrichKEGG(gene.list,organism="hsa",keyType='kegg', pvalueCutoff = pcutoff,
                    pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
  ####
  ekk2 <- setReadable(ekk, 'org.Hs.eg.db', 'ENTREZID')
  
  if(nrow(as.data.frame(ekk2)) > 0){
    
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
    ####
    write.csv(as.data.frame(ekk),paste(outdir, "/KEGG/KEGG_Enrichment.csv", sep = ''), row.names =F)
    #pathview of all enrichment KEGG pathways.
    library("pathview")
    ek = as.data.frame(ekk)
    if (dim(ek)[1] > 0 & dim(ek)[1] < 50){KGBubble(ek, outdir, KEGG = TRUE, width = 25, height = 20)
    }else{ KGBubble(ek[1:50,], outdir = paste(outdir,'KEGG',sep = '/'), KEGG = TRUE, width = 25, height = 20) }
    
    hsa_id <- ek$ID[-which(ek$ID %in% c("hsa04723", "hsa05206"))]
    for ( id in hsa_id ){
      pathview(gene.data=geneList,pathway.id=id, species="hsa", 
               kegg.dir = "D://Scripts/database/pathways/KEGG/", 
               limit = list(gene = 2, cpd = 2), kegg.native = TRUE)
    }
    file.remove(list.files()[grep('xml',list.files())])
    file.copy(list.files()[grep('pathview',list.files())], paste(outdir, '/KEGG/Pathview',sep = ''))
    file.remove(list.files()[grep('hsa',list.files())])
  }else{
    print("Warning: No KEGG pathways were enriched !")
  }
  ##########################################enrichGO###############################################
  mfgo <- enrichGO(gene = gene.list, OrgDb="org.Hs.eg.db",ont = "MF", pvalueCutoff = pcutoff, readable = TRUE)
  bpgo <- enrichGO(gene = gene.list, OrgDb="org.Hs.eg.db",ont = "BP", pvalueCutoff = pcutoff, readable = TRUE)
  ccgo <- enrichGO(gene = gene.list, OrgDb="org.Hs.eg.db",ont = "CC", pvalueCutoff = pcutoff, readable = TRUE)
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
 
    GO2 <- na.omit(rbind(mf[1:10,],bp[1:10,],cc[1:10,]))
    ####
    write.csv(GO, paste(outdir, "GO/GO_Enrichment.csv", sep = '/'), row.names=FALSE, quote = FALSE)
    KGBubble(GO2, paste(outdir, 'GO',sep = '/'), KEGG = FALSE, width = 25, height = 20)
    ####
    if(nrow(mf) >= 1){
      ggsave(paste(outdir, '/GO/GO_MF_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(mfgo, showCategory=10)+ ggtitle("Molecular Function of GO") +
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
      png(paste(outdir, "/GO/GO_MF_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
      cnetplot(mfgo,categorySize="pvalue",foldChange=geneList)
      dev.off()
      ####
      png(paste(outdir, "/GO/GO_MF_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
      emapplot(mfgo)
      dev.off()
      ####
      png(paste(outdir, "/GO/topGO_MF.png", sep = ''),width = 1200, height = 1200)
      suppressMessages(plotGOgraph(mfgo))
      dev.off()
      
    }
    ####
    if(nrow(bp) >= 1){
      ggsave(paste(outdir, '/GO/GO_BP_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(bpgo, showCategory=10)+ ggtitle("Biology Process of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 10, units = "in")
      ####
      ggsave(paste(outdir, '/GO/GO_BP_Enrichment_barplot.png', sep = ''),
             plot = barplot(bpgo, showCategory=10)+ ggtitle("Biology Process of GO") +
               theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                     plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
             width = 10, height = 6, units = "in")
      
      ####
      png(paste(outdir, "/GO/GO_BP_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
      cnetplot(bpgo,categorySize="pvalue",foldChange=geneList)
      dev.off()
      ####
      png(paste(outdir, "/GO/GO_BP_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
      emapplot(bpgo)
      dev.off()
      ####
      png(paste(outdir, "/GO/topGO_BP.png", sep = ''),width = 1200, height = 1200)
      suppressMessages(plotGOgraph(bpgo))
      dev.off()
    }
    ####
    if(nrow(cc) >= 1){
      ggsave(paste(outdir, '/GO/GO_CC_Enrichment_dotplot.png', sep = ''),
           plot = dotplot(ccgo, showCategory=10)+ ggtitle("Cell Circle of GO") +
             theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                   plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
           width = 10, height = 10, units = "in")
      
      ####
      ggsave(paste(outdir, '/GO/GO_CC_Enrichment_barplot.png', sep = ''),
             plot = barplot(ccgo, showCategory=10)+ ggtitle("Cell Circle of GO") +
               theme(axis.title = element_text(size=15, vjust=0.5, hjust=0.5),
                     plot.title = element_text(size=15, face = 'bold', vjust=0.5, hjust=0.5)),
             width = 10, height = 6, units = "in")
      ####
      png(paste(outdir, "/GO/GO_CC_Enrichment_cnetplot.png", sep = ''),width = 1000, height = 800)
      cnetplot(ccgo,categorySize="pvalue",foldChange=geneList)
      dev.off()
      ####
      png(paste(outdir, "/GO/GO_CC_Enrichment_emapplot.png", sep = ''),width = 1000, height = 800)
      emapplot(ccgo)
      dev.off()
      ####
      png(paste(outdir, "/GO/topGO_CC.png", sep = ''),width = 1200, height = 1200)
      suppressMessages(plotGOgraph(ccgo))
      dev.off()
    }
  ####
  # save.image(paste(outdir,'DEGenes_Enrichment.RData'))
}

