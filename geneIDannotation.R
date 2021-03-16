####FUNCTION I.
###############################################################################
####define a function to annotate all kinds of input gene ID/symbol types. ####
###############################################################################
geneIDannotation <- function(geneLists=c(1,2,9), name = T, org = 'hsa'){
  ## input ID type : So far I just accept entrezID or symbol
  ## default, we will annotate the entrezID and symbol, chromosone location and gene name
  if ( org == 'hsa'){
    if(!require("org.Hs.eg.db")){BiocManager::install('org.Hs.eg.db')}
    all_EG=mappedkeys(org.Hs.egSYMBOL)
    all_ES=mappedRkeys(org.Hs.egENSEMBL)
    all_SY=mappedRkeys(org.Hs.egSYMBOL)
    EG2name=toTable(org.Hs.egGENENAME)
    EG2Symbol=toTable(org.Hs.egSYMBOL)
    EG2Ensembl=toTable(org.Hs.egENSEMBL)
    if(length(geneLists[1:10] %in% all_SY)>1){
      inputType='symbol'
      geneLists=data.frame(symbol=geneLists)
      results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
      results=merge(results,EG2Ensembl,by='gene_id',all.x=T)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }else if(length(geneLists[1:10] %in% all_EG)>1){
      inputType='entrezID'
      geneLists=data.frame(gene_id=geneLists)
      results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
      results=merge(results,EG2Ensembl,by='gene_id',all.x=T)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }else if(length(geneLists[1:10] %in% all_ES) > 1){
      inputType='ensemblID'
      geneLists=data.frame(ensembl_id=geneLists)
      results=merge(geneLists,EG2Ensembl,by='ensembl_id')
      results=merge(results,EG2Symbol,by='gene_id')
      results=merge(results,EG2name,by='gene_id')
    }else{
      return(NA)
    }
    return(results)
  }else if ( org == 'mmu'){
    if(!require("org.Mm.eg.db")){BiocManager::install('org.Mm.eg.db')}
    all_EG=mappedkeys(org.Mm.egSYMBOL)
    all_ES=mappedRkeys(org.Mm.egENSEMBL)
    all_SY=mappedRkeys(org.Mm.egSYMBOL)
    EG2name=toTable(org.Mm.egGENENAME)
    EG2Symbol=toTable(org.Mm.egSYMBOL)
    EG2Ensembl=toTable(org.Mm.egENSEMBL)
    if((geneLists[10] %in% all_SY)){
      inputType='symbol'
      geneLists=data.frame(symbol=geneLists)
      results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
      results=merge(results,EG2Ensembl,by='gene_id',all.x=T)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }else if((geneLists[10] %in% all_EG)){
      inputType='entrezID'
      geneLists=data.frame(gene_id=geneLists)
      results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
      results=merge(results,EG2Ensembl,by='gene_id',all.x=T)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }else if(length(geneLists[10] %in% all_ES) >= 1){
      inputType='ensemblID'
      geneLists=data.frame(ensembl_id=geneLists)
      results=merge(geneLists,EG2Ensembl,by='ensembl_id')
      results=merge(results,EG2Symbol,by='gene_id')
      results=merge(results,EG2name,by='gene_id')
    }else{
      return(NA)
    }
    return(results)
  }
}
