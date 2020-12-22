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
KS_HC<-c()
KS_amdMCIr <- c()
KS_amdMCIe <- c()
KS_AD <- c()
thresholds <- c()

for (threshold in seq(from = 0.05, to = 0.3, by = 0.01)){
  thresholds <- c(thresholds, threshold)
  dtruncated_log_normal <- function(x, meanlog, sdlog) 
    dtrunc(x, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)
  ptruncated_log_normal <- function(q, meanlog, sdlog) 
    ptrunc(q, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)
  qtruncated_log_normal <- function(p, meanlog, sdlog) 
    qtrunc(p, "lnorm", a=threshold, meanlog=meanlog, sdlog=sdlog)

  ## Fit data to LND
  AdjustAndSave<-function(datos, name){
    fit <- fitdist(datos, "truncated_log_normal", start = list(meanlog=0, sdlog=1), method = "mle")
    gof<- gofstat(fit)
    ksdistance <- ks.test(datos, "ptruncated_log_normal", fit$estimate[1], fit$estimate[2], exact = FALSE)
    D <- ksdistance$statistic

    # Comprobar si datas sintenticos
    ntimes <- 1000
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
    pvalue <- pvalue/ntimes
    return(pvalue)
  }


  KS_HC <- c(KS_HC, AdjustAndSave(HC[HC>threshold], "HC"))
  KS_amdMCIr <- c(KS_amdMCIr, AdjustAndSave(amdMCIr[amdMCIr>threshold], "amdMCIr"))
  KS_amdMCIe <- c(KS_amdMCIe, AdjustAndSave(amdMCIe[amdMCIe>threshold], "amdMCIe"))
  KS_AD <- c(KS_AD, AdjustAndSave(AD[AD>threshold], "AD"))
  mean_statistic <- (KS_AD+KS_HC+ KS_amdMCIe + KS_amdMCIr)/4
  #print(mean_statistic)
  plot(thresholds, mean_statistic)
}


adjustcutoff <- data.frame("thresholds"=thresholds, "mean_statistic" = mean_statistic)
colnames(adjustcutoff) <- c("thresholds", "mean_statistic")
write.csv(adjustcutoff,"adjusts/cutoffselection.csv")
