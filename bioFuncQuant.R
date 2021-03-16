#######################################################################################
####                          通路或免疫细胞等的 ssgsea 分析                         ####
#######################################################################################
bioFuncQuant <- function(exprData, gmt){
  if(!require(GSVA)){BiocManager::install('GSVA')}
  # ssgsea analysis is based on GSVA R package.
  if(!is.list(gmt)){genesets <- gmt2list(gmt = gmt)}else{genesets = gmt}
  suppressMessages(res <- as.data.frame((t(gsva(as.matrix(exprData), genesets, method = "ssgsea", parallel.sz = 40)))))
  return(res)
  rm(progressBar);rm(iSample);rm(nSamples)
}


############################################################################
####          define a function to run ssGSEA(GSVA package).            ####
############################################################################
gmt2list <- function(gmt){
  gmtlist <- list()
  open(con <- file(gmt))
  i = 1
  term = c()
  line <- readLines(con, n = 1)
  while (length(line) > 0){
    ele <- stringr::str_split(line, '\t')[[1]]
    gmtlist[[i]] <- ele[3:length(ele)]
    term <- c(term, stringr::str_split(line, '\t')[[1]][1])
    i = i + 1
    line <- readLines(con, n = 1)
  }
  close(con)
  names(gmtlist) <- term
  return(gmtlist)
}



############################################################################
####                                                                    ####
############################################################################
kegg2gmt <- function(species = 'hsa', outdir = '.'){
  if(!require(gage)){BiocManager::install('gage');require(gage)}
  if(!require(stringr)){BiocManager::install('stringr');require(stringr)}
  kegg.gmt <- kegg.gsets(species = species, id.type = 'entrez', check.new = TRUE)
  headers  <- names(kegg.gmt$kg.sets)
  con_str  <- str_c(kegg.gmt$kg.sets[[headers[1]]], collapse='\t', sep='\t')
  con_str    = paste(headers[1],'\t\t',con_str,sep = '')
  write.table(con_str, file = paste(species,'_kegg.gmt',sep = ''),
              quote = FALSE, sep = '\t',append = FALSE,col.names = FALSE, row.names = FALSE)
  for (i in 2:length(headers)) {
    con_str = str_c(kegg.gmt$kg.sets[[headers[i]]], collapse = '\t', sep = '\t')
    con_str = paste(headers[i],'\t\t',con_str,sep = '')
    write.table(con_str, file = paste(species,'_kegg.gmt',sep = ''),
                quote = FALSE, sep = '\t',append = TRUE,col.names = FALSE, row.names = FALSE)
  }
}

