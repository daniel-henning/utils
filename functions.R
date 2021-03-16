# source("./utils/genelist.enrichment.R")
# source("./utils/DEgenes.enrichment.R")
# source("./utils/countNormlization.R")
# source("./utils/geneIDannotation.R")
# source("./utils/geneLocalAnno.R")
# source('./utils/DEG.analysis.R')
# source("./utils/BatchRemove.R")
# source("./utils/dupRemove.R")
# source("./utils/expr2gsea.R")
# source("./utils/fig2ppt.R")
# source('./utils/cytolyticCalc.R')
# source("./utils/geoSoftDataExtract.R")
# source("./utils/RD.RNA.Extraction.R")
# source("./utils/GEOsetsProcess.R")
# source("./utils/GenomeCompare.R")
# source("./utils/expr.dotplot.R")
# source("./utils/bioFuncQuant.R")
# source("./utils/scatter.plot.R")
# source("./utils/violin.plot.R")
# source("./utils/tsne.plot.R")
# source("./utils/KGBubble.R")
# source("./utils/gmt2list.R")
# source("./utils/wssplot.R")
# source('./utils/deg.limma.R')
# source("./utils/MIBCSubtyping.R")
# source("./utils/BulidProjTree.R")
# source("./utils/GEP_Analysis.R")

####FUNCTION I.
####################################################################################
####    normalize raw count expression data by deseq2 and dedup gene symbols    ####
####################################################################################
nmrc2sym <- function(rawCount, gtf = 'gencode.gene.info.v22.tsv'){
  raw.norm      <- log2(countNormlization(rawCount) + 1)
  raw.norm.anno <- geneLocalAnno(gtf = gtf, resDEG = raw.norm)
  raw.dedup     <- dupRemove(raw.norm.anno[,-ncol(raw.norm.anno)])
  return(raw.dedup[,-ncol(raw.dedup)])
}

####FUNCTION II.
####################################################################################
####define a function to read RNAseq count data from csv/txt/xlsx format files. ####
####################################################################################
loadCountMatrix <- function(input.file.name, check.names = FALSE, sep = "\t", rname = 1){
  if (length(grep('xlsx',input.file.name,value = T)) > 1){
    dat = read.xlsx(input.file.name,1)
  }else{
    dat = read.table(input.file.name,header = TRUE,sep = sep, check.names = check.names, row.names = rname)}
  return(dat)
}


####FUNCTION III.
####################################################################################
####           define a function to caculate mean count of a group.             ####
####################################################################################
meanCount <- function(DE_Results = DE_Results, i, j, colData = colData, factors = factors, exprSet = exprSet){
  index_grp1 <- rownames(colData)[which(colData$group == factors[i])]
  index_grp2 <- rownames(colData)[which(colData$group == factors[j])]
  DE_Results$grp1 <- NA
  DE_Results$grp2 <- NA
  for ( g in 1:dim(DE_Results)[1]){
    index.gene <- which(rownames(exprSet) == rownames(DE_Results)[g])
    DE_Results$grp1[g] <- mean(t(exprSet[index.gene,index_grp1]))
    DE_Results$grp2[g] <- mean(t(exprSet[index.gene,index_grp2]))
  }
  colnames(DE_Results) <- sub('grp1',paste('mean Count Of',factors[i],sep = ' '),colnames(DE_Results))
  colnames(DE_Results) <- sub('grp2',paste('mean Count Of',factors[j],sep = ' '),colnames(DE_Results))
  return(DE_Results)
}

#FUNCTION IV.
############################################################################################
####define a function to select duplication and return those elements in vector format. ####
############################################################################################
seleDup <- function(vect){
  vect.uniq    <- unique(vect)
  vect.nonuniq <- c()
  for (v in vect.uniq){
    index <- which(vect == v)
    if (length(index) > 1){
      print(paste(v,' is not uniq in this vector.', sep = ''))
      vect.nonuniq = c(vect.nonuniq, v)
    }
  }
  return(vect.nonuniq)
}

####
####Function V.
####################################################################################
####           get the gene fold-change list from limma DEG anaysis             ####
####################################################################################
geneRanks <- function(degenes){
  if (length(grep('logFC', colnames(degenes))) > 0){
    gene.ranks <- degenes$logFC
  }else if (length(grep('log2FoldChange', colnames(degenes))) > 0){
    gene.ranks <- degenes$log2FoldChange
  }
  names(gene.ranks) <- rownames(degenes)
  gene.ranks <- gene.ranks[which(gene.ranks != 'NA')]
  gene.ranks <- sort(gene.ranks,decreasing = T)
  #for (i in 1:length(gene.ranks)) { if(gene.ranks[i] > 0 ){
  #  gene.ranks[i] <- 2^abs(gene.ranks[i]) }else{gene.ranks[i] <- -2^abs(gene.ranks[i])}}
  #rm(list = c('i'))
  return(gene.ranks)
}

####Function VI.
####################################################################################
####           volcano plot of significant differential expression genes.       ####
####################################################################################
plot.volcano <- function(res,main){
  library(calibrate)
  # Make a basic volcano plot
  with(res, plot(logFC, -log10(P.Value), pch=20, main = main, xlim=c(-2.5,2)))
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(res, P.Value<.001 ), points(logFC, -log10(P.Value), pch=20, col="red"))
  with(subset(res, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
  with(subset(res, P.Value<.001 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="green"))
  # Label points with the textxy function from the calibrate plot
  with(subset(res, P.Value<.001 & abs(logFC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))
}

####Function VII.
####################################################################################
####              merge gene expression data from tcga data portal.             ####
####################################################################################
tcgaMerge <- function(rootdir = '.', sample.sheet, outdir = '.'){
  sample.sheet = read.table(sample.sheet,header = T, sep = '\t')
  #read.table(gzfile("file.text.gz")) # Requires you have ./file.text.bz2
  projects <- unique(sample.sheet$Project.ID)
  ########
  proj.expr.data <- matrix('TCGA.htseq.count', nrow = 60488, ncol = 1)
  for (project in projects) {
    samps <- sample.sheet[which(sample.sheet$Project.ID == project),]
    expr.data = list()
    sample.ID = as.character(samps$Sample.ID)
    ####
    empty.file = c()
    for (i in 1:length(samps$File.Name)) {
      f = paste(rootdir,samps[i,'File.ID'],samps[i,'File.Name'],sep = '/')
      if (file.exists(f)){
        expr.data[[i]] <- read.table(gzfile(f),header = FALSE,row.names = 1,sep = '\t')
        print(paste('The ',i,'th sample \'',samps$Sample.ID[i],'\' has been loaded...',sep = ''))
      }else{
        empty.file <- c(empty.file,i)
      }
    }
    ########
    index = c()
    for (exp in expr.data) { index <- union(index, rownames(exp))}
    proj = as.data.frame(expr.data[[1]][index,])
    rownames(proj) = index
    for (j in 2:length(expr.data)) {
      tmp = as.data.frame(expr.data[[j]][index,])
      rownames(tmp) = index
      proj <- cbind(proj,tmp) }
    if (length(empty.file) == 0){
      colnames(proj) <- sample.ID
    }else{
      colnames(proj) <- sample.ID[-(empty.file)]
    }
    rownames(proj) <- substring(rownames(proj),1,15)
    save(proj, file = paste(outdir,project,'.htseq.count.RData',sep = ''))
    proj.expr.data <- cbind(proj.expr.data, proj)
  }
  return(proj.expr.data)
}

####FUNCTION VIII.
####################################################################################
####                  function to select top varinces by MAD.                   ####
####################################################################################
mad.filter <- function(expr.matrix, top=1500, bottom=FALSE){
  if(bottom){
    expr.mad <- apply(expr.matrix,1,FUN = mad)
    idx <- order(expr.mad, decreasing = FALSE)
    return(dat = expr.matrix[idx[1:top],])
  }else{
    expr.mad <- apply(expr.matrix,1,FUN = mad)
    idx <- order(expr.mad, decreasing = TRUE)
    return(dat = expr.matrix[idx[1:top],])
  }
}


####FUNCTION IX.
####################################################################################
####                                                                            ####
####################################################################################
readFile <- function(file, sep = '\t', header = TRUE, rownm = NULL, bigFile = FALSE, check.names = FALSE, quote = ""){
  if (bigFile) {
    library(data.table)
    dat <- fread(file, header = header, sep = sep, quote = quote, check.names = check.names)
    dat <- as.data.frame(dat)
    if(typeof(rownm) == "NULL"){return(dat)
    }else{rownames(dat) <- dat[,rownm]; return(dat[,-rownm])}
  }else{return(read.table(file, header = header, sep = sep, row.names = rownm, quote = quote, check.names = check.names))}
}

####FUNCTION X.
####################################################################################
####                                                                            ####
####################################################################################
getGeneMut.TCGA <- function(maf = maf, gene = gene){
  require(data.table)
  maf <- as.data.frame(fread(maf, header = TRUE, sep = '\t'))
  mut.samp  <- unique(maf$Tumor_Sample_Barcode[which(maf$Hugo_Symbol == gene)])
  wild.samp <- setdiff(unique(maf$Tumor_Sample_Barcode), mut.samp)
  return(list(mut.samples = mut.samp, wild.samples = wild.samp))
}

####FUNCTION XI.
####Remove Duplicated Genes expression with same geneID.
exprGeneDedup <- function(dat){ #the first column is the geneID.
  dup.genes <- as.character(unique(dat[,1][which(duplicated(dat[,1]))]))
  uni.genes <- setdiff(dat[,1], dup.genes)
  dat.out <- dat[which(dat[,1] %in% uni.genes),]
  dat.out[(nrow(dat.out)+1):(nrow(dat.out)+length(dup.genes)),1] <- dup.genes
  ####
  for (gene in dup.genes) {
    expr.data <- dat[which(dat[,1] == gene),]
    for (i in 2:ncol(dat.out)) {
      dat.out[which(dat.out[,1] == gene),i] <- max(expr.data[,i])
    }
  }
  dat.out <- as.data.frame(dat.out)
  rownames(dat.out) <- 1:nrow(dat.out)
  return(dat.out)
}

####FUNCTION XII.
####################################################################################
####         make a data tree strcture for a new RNAseq-analysis project.       ####
####################################################################################
mk.proj.tree <- function(dir){
  if(!dir.exists(dir)){dir.create(dir)}
  setwd(dir = dir)
  if(!dir.exists("./data")){dir.create("./data")}
  if(!dir.exists("./data/QC/")){dir.create("./data/QC/")}
  if(!dir.exists("./data/quantification/")){dir.create("./data/quantification/")}
  if(!dir.exists("./data/gene.fusion/")){dir.create("./data/gene.fusion/")}
  if(!dir.exists("./data/DEGanalysis/")){dir.create("./data/DEGanalysis/")}
  if(!dir.exists("./data/immuneAnalysis/")){dir.create("./data/immuneAnalysis/")}
  if(!dir.exists("./data/Enrichment/")){dir.create("./data/Enrichment//")}
  if(!dir.exists("./data/Enrichment/GO/")){dir.create("./data/Enrichment/GO/")}
  if(!dir.exists("./data/Enrichment/KEGG/")){dir.create("./data/Enrichment/KEGG/")}
  if(!dir.exists("./data/Enrichment/GSEA/")){dir.create("./data/Enrichment/GSEA/")}
  if(!dir.exists("./data/Enrichment/GSEA/GO/")){dir.create("./data/Enrichment/GSEA/GO/")}
  if(!dir.exists("./data/Enrichment/GSEA/KEGG/")){dir.create("./data/Enrichment/GSEA/KEGG/")}
  if(!dir.exists("./data/Enrichment/GSEA/HALLMARK/")){dir.create("./data/Enrichment/GSEA/HALLMARK/")}
  if(!dir.exists("./data/Enrichment/GSEA/REACTOME/")){dir.create("./data/Enrichment/GSEA/REACTOME/")}
  system(paste("chmod -R 777 ",abs.dir, '/data', sep = ''))
}


#FUNCTION XiII.
############################################################################
####define a function to transform a expreesion matrix to .gct and .cls.####
############################################################################
expr2gsea <- function(exprSet, group, outdir = './', samp = 'sample'){
  require(stringr)
  Descript <- as.data.frame((rep('na', nrow(exprSet))))
  rownames(Descript) <- rownames(exprSet)
  colnames(Descript) <- 'DESCRIPT'
  gct <- cbind(Descript, exprSet)
  cls = paste(paste(ncol(exprSet)," 2 1\n#Group1 Group2\n", sep = ''),
              paste(str_c(group, collapse = " "), '\n',sep = ' '), sep = '')
  write.table(cls,paste(outdir, samp, '.cls', sep = ''),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '')
  header = paste("#1.2\n", paste(nrow(exprSet), ncol(exprSet), sep = '\t'), '\nNAME\tDESCRIPT\t', sep = '')
  header = paste(header, str_c(colnames(exprSet), collapse = '\t'), sep = '')
  write.table(header, paste(outdir, samp, '.gct', sep = ''),
              quote = FALSE, sep = '', row.names = FALSE, col.names = FALSE)
  write.table(gct, paste(outdir, samp, '.gct', sep = ''),
              quote = FALSE, sep = '\t',append = TRUE,col.names = FALSE)
}

####################################################################################
####       Functions to intersect or unoin multi vectors in a list object.      ####
####################################################################################
multi.intersect<-function(dat.li){inter<-c();for(l in 2:length(dat.li)){inter<-intersect(dat.li[[l-1]],dat.li[[l]])};return(inter)}
multi.union<-function(dat.li){uni<-c();for(l in 2:length(dat.li)){uni<-union(dat.li[[l-1]],dat.li[[l]])};return(uni)}

####################################################################################
####                                                                            ####
####################################################################################
loadSampleInfo <- function(sampleInfo){sample.info=read.table(sampleInfo, header = TRUE, sep = '\t', check.names = FALSE);return(sample.info)}
loadExprData   <- function(exprData){expr=read.table(exprData, header = T, sep = '\t', row.names = 1, check.names = FALSE);return(expr)}

####################################################################################
####                            数据标准化处理算法                              ####
####################################################################################
min.max.norm <- function(x){return((x-min(x))/(max(x) - min(x)))}
zscore <- function(x){return((x-mean(x))/sd(x))}


####################################################################################
####                              基因表达分组                                  ####
####################################################################################
strata.gene <- function(exp.dat, gene, samples, top = 0.75, bottom = 0.25, med.label = NA){
  return(ifelse(t(exp.dat[gene,samples]) > quantile(t(exp.dat[gene,samples]),top),"High", ifelse(t(exp.dat[gene,samples]) <= quantile(t(exp.dat[gene,samples]),bottom),"Low",med.label))[,1])
}
