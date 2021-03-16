cytolyticCalc <- function(exprMatrix){
  #######################################################################################################
  #### the immune cytolytic activity was calculated as the geometric mean of the genes GZMA and PRF1, 
  #### 《Molecular and genetic properties of tumors associated with local immune cytolytic activity》
  #### 《The expression and prognostic impact of immune cytolytic activity-related markers 
  ####  in human malignancies: a comprehensive meta-analysis》
  #######################################################################################################
  
  if("GZMA" %in% rownames(exprMatrix) & "PRF1" %in% rownames(exprMatrix)){
    
    gzma <- exprMatrix["GZMA",]
    
    prf1 <- exprMatrix["PRF1",]
    
    cytolytic <- as.data.frame(t(sqrt(gzma*prf1)))
    
    colnames(cytolytic) <- "cytolytic activity"
    
    return(cytolytic)
    
  }else{
    
    cat("Your input data do not have gene expression value of GZMA or PRF1.")
    
  }
  
}
