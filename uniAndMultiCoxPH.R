# multi-variables cox-ph model
multivarCoxph <- function(covariates, data, time = "time", status = "status"){
  res.cox <- coxph(as.formula(paste(paste('Surv(',time,',',status,')~'),paste(covariates,collapse = " + "))), data = data)
  x <- summary(res.cox)
  p.value<- signif(x$coefficients[,5], digits=2)
  wald.test<-signif(x$wald["test"], digits=2)
  beta <- signif(x$coef[,1], digits=2);#coeficient beta
  HR   <- signif(x$coef[,2], digits=2);#exp(beta)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  HR <- paste0(HR, " (",
               HR.confint.lower, "-", HR.confint.upper, ")")
  res <- data.frame(beta, HR, p.value, wald.test)
  res$C.index <- x$concordance[1]
  res$C.index.se <- x$concordance[2]
  colnames(res) <- c("beta", "HR (95% CI for HR)","p.value", "wald.test","C.index","se(C)")
  return(res)
}

# uni-variables cox-ph model
univarCoxph <- function(covariates, data, time = "time", status = "status"){
  univ_formulas <- sapply(covariates, function(x) as.formula(paste(paste('Surv(',time,',',status,')~'), x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
  # Extract data
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  return(as.data.frame(res))
}

