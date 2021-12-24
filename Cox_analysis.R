#### FUNCTION TO RUN SINGLE VARIABLE COXPH ANALYSIS.
uni_coxph <- function(dat, variables, time = "time", status = "status"){
  if(!require(survival)){install.packages("survival");require(survival)}
  if(!require(survminer)){install.packages("survminer");require(survminer)}
  
  y <- Surv(time=dat[,time],event=dat[,status]==1)

  Uni_cox_model<- function(x){
    FML <- as.formula(paste0 ("y~",x))
    cox<- coxph(FML,data=dat)
    cox1<-summary(cox)
    HR <- round(cox1$coefficients[,2],3)
    PValue <- round(cox1$coefficients[,5],4)
    CI5 <-round(cox1$conf.int[,3],3)
    CI95 <-round(cox1$conf.int[,4],3)
    Uni_cox_model<- data.frame( "variables" = x, 'HR' = HR, 'CI.95' = paste(CI5,"-",CI95), 'p' = PValue)
    return(Uni_cox_model)
  }
  #Runing Uni-cox analysis.
  Uni_cox <- ldply(lapply(variables, Uni_cox_model))
  
  return(Uni_cox)
}



#### FUNCTION TO RUN MULTI-VARIABLES COXPH ANALYSIS.
mul_coxph <- function(dat, covariates, time = "time", status = "status"){
  if(!require(survival)){install.packages("survival");require(survival)}
  if(!require(survminer)){install.packages("survminer");require(survminer)}
  univ_formulas <- sapply(covariates, function(x)as.formula(paste(paste0('Surv(',time,',', status,')~'), x)))

  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dat)})
  
  # Extract data 
  univ_results <- lapply(univ_models, function(x){ x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           CI <- paste0(HR.confint.lower, " - ", HR.confint.upper)
                           res<-c(beta, HR, CI, wald.test, p.value)
                           names(res)<-c("beta", "HR", "CI(95)", "wald.test", "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(as.data.frame(res))
}

