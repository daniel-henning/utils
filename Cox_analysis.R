
#### FUNCTION TO RUN MULTI-VARIABLES COXPH ANALYSIS.
mul_coxph <- function(dat, covariates, time = "time", status = "status"){

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
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(as.data.frame(res))
}

