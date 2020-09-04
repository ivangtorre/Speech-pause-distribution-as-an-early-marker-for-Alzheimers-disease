# This script replicates the analysis and data fit to candidate
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


## Import data
MyRData <- read.csv2("FullData.csv", header = TRUE, sep = ";")
unique(MyRData$Grupo)


## Subset data into groups
threshold <- 0.249

CS <- subset(MyRData, MyRData$Grupo == "CS")$duration
CS <- CS[CS>threshold]

DCLej <- subset(MyRData, MyRData$Grupo == "DCLej")$duration
DCLej <- DCLej[DCLej>threshold]

DCLtm <- subset(MyRData, MyRData$Grupo == "DCLtm")$duration
DCLtm <- DCLtm[DCLtm>threshold]

DTA <- subset(MyRData, MyRData$Grupo == "DTA")$duration
DTA <- DTA[DTA>threshold]


## Truncated distribution to fit
library(fitdistrplus)
library(truncdist)
dtruncated_log_normal <- function(x, meanlog, sdlog) 
  dtrunc(x, "lnorm", a=.25, meanlog=meanlog, sdlog=sdlog)
ptruncated_log_normal <- function(q, meanlog, sdlog) 
  ptrunc(q, "lnorm", a=.25, meanlog=meanlog, sdlog=sdlog)

dtruncated_gamma <- function(x, shape, scale) 
  dtrunc(x, "gamma", a=.25, shape=shape, scale=scale)
ptruncated_gamma <- function(q, shape, scale) 
  ptrunc(q, "gamma", a=.25, shape=shape, scale=scale)

dtruncated_gamma <- function(x, shape, scale) 
  dtrunc(x, "gamma", a=.25, shape=shape, scale=scale)
ptruncated_gamma <- function(q, shape, scale) 
  ptrunc(q, "gamma", a=.25, shape=shape, scale=scale)

dtruncated_weibull <- function(x, shape, scale) 
  dtrunc(x, "weibull", a=.25, shape=shape, scale=scale)
ptruncated_weibull <- function(q, shape, scale) 
  ptrunc(q, "weibull", a=.25, shape=shape, scale=scale)


## Fit data to candidate distributions. Goodness of fit
goodnesoffit <- function(datos){
  fit1 <- fitdist(datos, "truncated_log_normal", start = list(meanlog=0, sdlog=1), method = "mle")
  fit2 <- fitdist(datos, "truncated_gamma", method = "mle",start = list(shape=1, scale=1))
  fit3 <- fitdist(datos, "truncated_weibull", method = "mle",start = list(shape=1, scale=1))
  return(round(c(fit1$aic, fit1$bic, fit2$aic, fit2$bic, fit3$aic, fit3$bic)))
}
AIC_BIC <- t(data.frame("CS"=goodnesoffit(CS),
                    "DCLej" = goodnesoffit(DCLej),
                    "DCLtm" = goodnesoffit(DCLtm),
                    "DTA" = goodnesoffit(DTA)))
colnames(AIC_BIC) <- c("LND_AIC", "LND_BIC","GAM_AIC", "GAM_BIC","WEIB_AIC", "WEIB_BIC")
write.csv(AIC_BIC,"adjusts/AIC_BIC.csv")


## Fit data to LND
AdjustAndSave<-function(datos, name){
  fit <- fitdist(datos, "truncated_log_normal", start = list(meanlog=0, sdlog=1), method = "mle")
  gof<- gofstat(fit)
  dir.create(file.path(getwd(), "adjusts"), showWarnings = FALSE)
  write.csv(data.frame(fit$estimate, fit$sd), paste("adjusts/fit_",name,".csv"))
  write.csv(data.frame(gof$aic, gof$bic, fit$loglik, gof$ks, gof$kstest), paste("adjusts/gof_",name,".csv"))
}

AdjustAndSave(CS, "CS")
AdjustAndSave(DCLej, "DCLej")
AdjustAndSave(DCLtm, "DCLtm")
AdjustAndSave(DTA, "DTA")

