# This scripts the analysis and data fit to candidate
# PDFs addresed in the research paper:
#
# 
# P. Pastoriza, I.G. Torre, F. Diéguez, I. Gomez, S. Gelado, J. Bello, 
# A. Ávila, J. Matías, V. Pytell, A. Hernandez-Fernandez (2020) Speech
# pause distribution as an early marker for Alzheimer’s disease. Submitted
#
#
# Maintainer <- Ivan G Torre
# Contact <- igonzalez@vicomtech.org
library(fitdistrplus)
library(truncdist)
options(warn=-1)


## Import data
MyRData <- read.csv2("listaPausasUmral50_FINAL.csv", header = TRUE, sep = ",")
MyRData$durPausa <-as.numeric(as.matrix(MyRData$durPausa))
unique(MyRData$Diagnóstico)

HC <- subset(MyRData, MyRData$Diagnóstico == "HC")$durPausa
amdMCIr <- subset(MyRData, MyRData$Diagnóstico == "a-mdMCI-R")$durPausa
amdMCIe <- subset(MyRData, MyRData$Diagnóstico == "a-mdMCI-E")$durPausa
AD <- subset(MyRData, MyRData$Diagnóstico == "AD")$durPausa

## Truncated distribution to fit
threshold <- 0.16

dtruncated_log_normal <- function(x, meanlog, sdlog) 
  dtrunc(x, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)
ptruncated_log_normal <- function(q, meanlog, sdlog) 
  ptrunc(q, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)
qtruncated_log_normal <- function(p, meanlog, sdlog) 
  qtrunc(p, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)

dtruncated_gamma <- function(x, shape, scale) 
  dtrunc(x, "gamma", a=threshold, shape=shape, scale=scale)
ptruncated_gamma <- function(q, shape, scale) 
  ptrunc(q, "gamma", a=threshold, shape=shape, scale=scale)

dtruncated_weibull <- function(x, shape, scale) 
  dtrunc(x, "weibull", a=threshold, shape=shape, scale=scale)
ptruncated_weibull <- function(q, shape, scale) 
  ptrunc(q, "weibull", a=threshold, shape=shape, scale=scale)



## Fit data to candidate distributions. Goodness of fit
goodnesoffit <- function(datos){
  fit1 <- fitdist(datos, "truncated_log_normal", start = list(meanlog=0, sdlog=1), method = "mle")
  fit2 <- fitdist(datos, "truncated_gamma", method = "mle",start = list(shape=1, scale=1))
  fit3 <- fitdist(datos, "truncated_weibull", method = "mle",start = list(shape=1, scale=1))
  return(round(c(fit1$aic, fit1$bic, fit2$aic, fit2$bic, fit3$aic, fit3$bic)))
}


## Fit data to LND
AdjustAndSave<-function(datos, name){
  fit <- fitdist(datos, "truncated_log_normal", start = list(meanlog=0, sdlog=1), method = "mle")
  gof<- gofstat(fit)
  print(gof$kstest)
  dir.create(file.path(getwd(), "adjusts"), showWarnings = FALSE)
  write.csv(data.frame(fit$estimate, fit$sd), paste("adjusts/fit_",name,".csv"))
  #print(fit$estimate[1])
  print(fit)
  ksdistance <- ks.test(datos, "ptruncated_log_normal", fit$estimate[1], fit$estimate[2], exact = FALSE)
  D <- ksdistance$statistic
  #print(D)
  #print(ksdistance)
  
  # Comprobar si datas sintenticos
  ntimes <- 100
  pvalue <- 0
  for (i in 1:ntimes) {
    stochas_data <- rtrunc(length(datos), spec="truncated_log_normal", meanlog=fit$estimate[1], sdlog=fit$estimate[2])
    stochast_ks <- ks.test(stochas_data, "ptruncated_log_normal", fit$estimate[1], fit$estimate[2], exact = FALSE)
    stochas_D <- stochast_ks$statistic
    # Count times stochastic data is greater than data
    if (stochas_D<D){
      pvalue <- pvalue +1
    }
  }
  #pvalue <- pvalue/ntimes
  #print("KS test:")
  #print(pvalue)
  #print(ksdistance$D)
  write.csv(data.frame(gof$aic, gof$bic, fit$loglik, gof$ks, gof$kstest), paste("adjusts/gof_",name,".csv"))
}


AIC_BIC <- t(data.frame("HC"=goodnesoffit(HC[HC>threshold]),
                    "amdMCIe" = goodnesoffit(amdMCIe[amdMCIe>threshold]),
                    "amdMCIr" = goodnesoffit(amdMCIr[amdMCIr>threshold]),
                    "AD" = goodnesoffit(AD[AD>threshold])))
colnames(AIC_BIC) <- c("LND_AIC", "LND_BIC","GAM_AIC", "GAM_BIC","WEIB_AIC", "WEIB_BIC")
write.csv(AIC_BIC,"adjusts/AIC_BIC.csv")

AdjustAndSave(HC[HC>threshold], "HC")
length(HC[HC>threshold])
AdjustAndSave(amdMCIr[amdMCIr>threshold], "amdMCIr")
length(amdMCIr[amdMCIr>threshold])
AdjustAndSave(amdMCIe[amdMCIe>threshold], "amdMCIe")
length(amdMCIe[amdMCIe>threshold])
AdjustAndSave(AD[AD>threshold], "AD")
length(AD[AD>threshold])
