# Mariam Fofana
# October 2016
# Modeling analysis of drivers of TB drug resistance trajectories
# Data analysis code
# Calibration of overall TB prevalence and incidence

# Change working directory as necessary
setwd("C:/Users/mariam/Dropbox/TB DST Model")

# Load pacakges and source necessary files
if(!require(abind)){
  install.packages("abind")
  library(abind)
}
if(!require(fmsb)){
  install.packages("fmsb")
  library(fmsb)
}
if(!require(ggm)){
  install.packages("ggm")
  library(ggm)
}
if(!require(pROC)){
  install.packages("pROC")
  library(pROC)
}
if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}

source("R code/pcor.R")
source("R code/AnalysisFxnsDetRevis.R")


# I: LOAD & FORMAT DATA--------------------------------------------------------

# Load data from formatted files
# Baseline simulations
initALL <- as.matrix(read.csv("Most recent output/Output/InputAllBsl_110915.csv"))
init <- as.matrix(read.csv("Most recent output/Output/InputSelBsl_110915.csv"))
initALLb <- as.matrix(read.csv("Most recent output/Output/InputAllNP_110915.csv"))
initb <- as.matrix(read.csv("Most recent output/Output/InputSelNP_110915.csv"))

dataALL <- as.matrix(read.csv("Most recent output/Output/DataAllBsl_110915.csv"))
dataSA1 <- as.matrix(read.csv("Most recent output/Output/DataSelBsl_110915.csv")) 

# dataSA2 <- as.matrix(read.csv("Output/SA/dstData110915_AllPZA.csv"))
# dataSA3 <- as.matrix(read.csv("Output/SA/dstData110915_F-RF.csv"))
# dataSA4 <- as.matrix(read.csv("Output/SA/dstData110915_R-RF.csv"))
# dataSA5 <- as.matrix(read.csv("Output/SA/dstData110915_P-RP.csv"))
# dataSA6 <- as.matrix(read.csv("Output/SA/dstData110915_R-RP.csv"))
# dataSA7 <- as.matrix(read.csv("Output/SA/dstData110915_FP-RFP.csv"))
# dataSA8 <- as.matrix(read.csv("Output/SA/dstData110915_RP-RFP.csv"))
# dataSA9 <- as.matrix(read.csv("Output/SA/dstData110915_R-RP-RFP.csv"))
# dataSA10 <- as.matrix(read.csv("Output/SA/dstData110915_P-RP-RFP.csv"))
# 
# dataNP <- as.matrix(read.csv("Output/DataSelNP_110915.csv"))

# Total number of SA scenarios 
SAtot <- 1
# 
# Check dimensions (if csv conversion has added extra column, remove it)
if(ncol(initALL) > 74) {initALL <- initALL[, -1]}
if(ncol(init) > 74) {init <- init[, -1]}
if(ncol(initALLb) > 74) {initALLb <- initALLb[, -1]}
if(ncol(initb) > 74) {initb <- initb[, -1]}

if(ncol(dataALL) > 6292) {dataALL <- dataALL[, -1]}
if(ncol(dataNP) > 6292) {dataNP <- dataNP[, -1]}

for (i in 1:SAtot) {
  eval(parse(text=
               paste(
                 "if(ncol(dataSA", i, ")>6292){dataSA", i, "<-dataSA", i,"[,-1]}", sep=""
               )
  ))
}
# 
# # Separate contents of "Input" files
# # These files contain  (1) 34 randomly sampled input variables
# #                     (2) 12 summary outputs corresponding to epidemiologic values
# #                     (3) 12 "calibration scores" (0=no; 1=yes) 
# 
# # Remove unused parameters (sampled but not actually used in model)
# init <- init[, -c(24, 25, 28, 30, 35:49)] 
# initALL <- initALL[, -c(24, 25, 28, 30, 35:49)]
# initb <- initb[, -c(24, 25, 28, 30, 35:49)] 
# initALLb <- initALLb[, -c(24, 25, 28, 30, 35:49)]
# 
# inputsel <- init[, 1:31]  # input variables 
# outputsel <- init[, 32:43]  # summary outputs 
# scoresel <- init[, 44:55]  # calibration scores 
# inputselb <- initb[, 1:31]  # input variables 
# outputselb <- initb[, 32:43]  # summary outputs 
# scoreselb <- initb[, 44:55]  # calibration scores 
# 
# inputall <- initALL[, 1:31]  # input variables
# outputall <- initALL[, 32:43]  # summary outputs
# scoreall <- initALL[, 44:55]  # calibration scores
# inputallb <- initALLb[, 1:31]  # input variables
# outputallb <- initALLb[, 32:43]  # summary outputs
# scoreallb <- initALLb[, 44:55]  # calibration scores

# Reconstruct trajectory data into 3D arrays
# All trajectories
DSTlhsData <- Convert3D(dataALL, 52, 121)
dim(DSTlhsData)
# DSTlhsDataNP <- Convert3D(dataNP, 52, 121)
# dim(DSTlhsDataNP)

# Selected trajectories and SAs
for (i in 1:SAtot) {
  eval(parse(
    text=paste("DSTlhsDataSA", i, "<-Convert3D(dataSA", i, ",52,121)", sep="")))
}






# II: SIMULATION SELECTION CRITERIA-------------------------------------------

# Identify calibration criteria most influential in selecting simulations
# Produce a table tallying the number/% of trajectories that are
# too low, on target, or too high for each of the criteria

totsamp <- nrow(initALL)  # total number of simulations
selsamp <- nrow(init)  # number of selected simulations
totsampb <- nrow(initALLb)  # total number of simulations
selsampb <- nrow(initb)  # number of selected simulations

# Extract output
criteria <- outputall
criteriab <- outputallb
colnames(criteria) <- c("inc13", "inc10", "inc05", "inc00", "inc95", "inc90",  # Incidence 
                        "mdrrerx13", "mdrnew13",  # MDR in retreatment & new cases
                        "moxfinal13", "pzafinal13",   # PZA/FQ monoresistance, new cases
                        "rpfinal13", "rffinal13")  # PZA/FQ resistance among MDR

criteria2 <- matrix(NA, 3, 12)
rownames(criteria2) <- c("low", "on target", "high")
colnames(criteria2) <- colnames(criteria)
criteria2b <- criteria2
colnames(criteria2b) <- c("inc13", "inc10", "inc05", "inc00", "inc95", "inc90",  # Incidence 
                          "mdrrerx13", "mdrnew13",  # MDR in retreatment & new cases
                          "pzafinal13", "moxfinal13",  # PZA/FQ monoresistance, new cases
                          "rpfinal13", "rffinal13") 

# Calibration bounds
lo <- c(c(183, 194, 213, 220, 218, 218) * 0.75, c(16, 2.2) * 0.5, 0, 0, 40, 10)
hi <- c(c(183, 194, 213, 220, 218, 218) * 1.25, c(16, 2.2) * 1.5, 5, 5, 70, 40)







# Identify # of simulations with epi value above, below, within range
for (i in 1:12) {
  # Get index #s for simulations in each category
  tmp1 <- which(criteria[, i] < lo[i])  # too low
  tmp2 <- which(criteria[, i] <= hi[i] & criteria[, i] >= lo[i])  # in range
  tmp3 <- which(criteria[, i] > hi[i])  # too high
  tmp1b <- which(criteriab[, i] < lo[i])  # too low
  tmp2b <- which(criteriab[, i] <= hi[i] & criteriab[, i] >= lo[i])  # in range
  tmp3b <- which(criteriab[, i] > hi[i])  # too high
  
  # Get # of simulations in each category
  tmp1 <- length(tmp1)  # too low
  tmp2 <- length(tmp2)  # in range
  tmp3 <- length(tmp3)  # too high
  tmp1b <- length(tmp1b)  # too low
  tmp2b <- length(tmp2b)  # in range
  tmp3b <- length(tmp3b)  # too high
  
  # Fill in table
  criteria2[1,i] <- tmp1
  criteria2[2,i] <- tmp2
  criteria2[3,i] <- tmp3
  criteria2b[1,i] <- tmp1b
  criteria2b[2,i] <- tmp2b
  criteria2b[3,i] <- tmp3b
}

# Check that all columns add up to total number of simulations
if (sum(colSums(criteria2) == rep(totsamp, 12)) != 12) {
  print(colSums(criteria2))
}
if (sum(colSums(criteria2b) == rep(totsampb, 12)) != 12) {
  print(colSums(criteria2b))
}

criteria3 <- criteria2 / totsamp  # percentages
criteria3b <- criteria2b / totsampb  # percentages

# Write files to csv
write.csv(criteria2, "Tables/criterianum.csv")
write.csv(criteria2b, "Tables/criterianum_NP.csv")
write.csv(criteria3, "Tables/criteriapct.csv")
write.csv(criteria3b, "Tables/criteriapct_NP.csv")

# Identify simulations that meet all incidence criteria combined
allinc <- matrix(NA, 2, 3)
rownames(allinc) <- c("Number", "Percent")
colnames(allinc) <- c("Lo", "In range", "Hi")
allincb <- allinc

# Too low for >=1 incidence value
allinc[1, 1] <- length(which(criteria[, 1] < lo[1] |
                               criteria[, 2] < lo[2] |
                               criteria[, 3] < lo[3] |
                               criteria[, 4] < lo[4] |
                               criteria[, 5] < lo[5] |
                               criteria[, 6] < lo[6] ))
allincb[1, 1] <- length(which(criteriab[, 1] < lo[1] |
                                criteriab[, 2] < lo[2] |
                                criteriab[, 3] < lo[3] |
                                criteriab[, 4] < lo[4] |
                                criteriab[, 5] < lo[5] |
                                criteriab[, 6] < lo[6] ))

# On target for all incidence values
allinc[1, 2] <- length(which(criteria[, 1] <= hi[1] & criteria[, 1] >= lo[1] &
                               criteria[, 2] <= hi[2] & criteria[, 2] >= lo[2] &
                               criteria[, 3] <= hi[3] & criteria[, 3] >= lo[3] &
                               criteria[, 4] <= hi[4] & criteria[, 4] >= lo[4] &
                               criteria[, 5] <= hi[5] & criteria[, 5] >= lo[5] &
                               criteria[, 6] <= hi[6] & criteria[, 6] >= lo[6] ))
allincb[1, 2] <- length(which(criteriab[, 1] <= hi[1] & criteriab[, 1] >= lo[1] &
                                criteriab[, 2] <= hi[2] & criteriab[, 2] >= lo[2] &
                                criteriab[, 3] <= hi[3] & criteriab[, 3] >= lo[3] &
                                criteriab[, 4] <= hi[4] & criteriab[, 4] >= lo[4] &
                                criteriab[, 5] <= hi[5] & criteriab[, 5] >= lo[5] &
                                criteriab[, 6] <= hi[6] & criteriab[, 6] >= lo[6] ))

# Too high for >=1 incidence value
allinc[1, 3] <- length(which(criteria[, 1] > hi[1] |
                               criteria[, 2] > hi[2] |
                               criteria[, 3] > hi[3] |
                               criteria[, 4] > hi[4] |
                               criteria[, 5] > hi[5] |
                               criteria[, 6] > hi[6] ))
allincb[1, 3] <- length(which(criteriab[, 1] > hi[1] |
                                criteriab[, 2] > hi[2] |
                                criteriab[, 3] > hi[3] |
                                criteriab[, 4] > hi[4] |
                                criteriab[, 5] > hi[5] |
                                criteriab[, 6] > hi[6] ))

allinc[2, ] <- allinc[1, ] / totsamp  # percentages
allincb[2, ] <- allincb[1, ] / totsampb  # percentages

write.csv(allinc, "Tables/criteriainc.csv")
write.csv(allincb, "Tables/criteriainc_NP.csv")

remove(criteria, criteriab)

# Distribution of selection criteria values (with calibration range)
histlab <- c("Incidence 2013 (137-229)",
             "Incidence 2010 (145-242)",
             "Incidence 2005 (160-266)",
             "Incidence 2000 (165-275)",
             "Incidence 1995 (163-272)",
             "Incidence 1990 (163-272)",
             "% MDR, retreatment (8-24%)",
             "% MDR, new (1.1-3.3%)",
             "% FQ monoresistant, new (<RIF)",
             "% PZA monoresistant, new (<RIF)",
             "% PZAr among MDR, retreatment (40-70%)",
             "% FQr among MDR, retreatment (10-40%)")

histlo <- c(rep(50, 6), rep(0, 6))
histhi <- c(rep(500, 6), 80, 30, 5, 20, 100, 15)

par(mfrow=c(3, 4), mai=c(1, 1, 0.2, 0.2))
for (i in 1:12) {
  h <- hist(outputall[, i], breaks=50, plot=F)
  clr <- ifelse((h$breaks >= lo[i] & h$breaks <= hi[i]), 
                "grey", "white")[-length(h$breaks)]
  plot(h, main="", col=clr, xlim=c(histlo[i], histhi[i]), xlab=histlab[i])
}
dev.copy(png, "Plots/criteriahist.png", width=1100, height=450)
dev.off()

breaklist <- c(rep(50, 8), 200, rep(50, 3))
par(mfrow=c(3, 4), mai=c(1, 1, 0.2, 0.2))
for (i in 1:12) {
  h <- hist(outputallb[, i], breaks=breaklist[i], plot=F)
  clr <- ifelse((h$breaks >= lo[i] & h$breaks <= hi[i]), 
                "grey", "white")[-length(h$breaks)]
  plot(h, main="", col=clr, xlim=c(histlo[i], histhi[i]), xlab=histlab[i])
}
dev.copy(png, "Plots/criteriahist_NP.png", width=1100, height=450)
dev.off()


# III: DISTRIBUTION OF INPUT VALUES--------------------------------------------

# Compare distribution of input values before and after selection
# Parameter names
varnames <- c(expression(t[1]), expression(t[2]), expression(eta["DS,RIFr | HRZE"]),
              expression(eta["DS,FQr | HRZE"]), expression(eta["DS,PZAr | HRZE"]),
              expression(f[RIFr]), expression(f[FQr]), expression(f[PZAr]),
              expression(f["RIF/FQr"]), expression(f["RIF/PZAr"]),
              expression(f["FQ/PZAr"]), expression(f["RIF/FQ/PZAr"]),
              expression(c["DS | HRZE"]), expression(c["RIFr | HRZE"]),
              expression(c["PZAr | HRZE"]), expression(c["RIF/PZAr | HRZE"]),
              expression(c["RIFr | STR"]), expression(c["RIF/FQr | STR"]),
              expression(c["RIF/PZAr | STR"]), expression(c["RIF/FQ/PZAr | STR"]),
              expression(m^R), expression(eta["RIFr,RIF/FQr | HRZE"]),
              expression(eta["RIFr,RIF/PZAr | HRZE"]), 
              expression(eta["PZAr,RIF/PZAr | HRZE"]),
              expression(eta["PZAr,FQ/PZAr | HRZE"]),
              expression(eta["RIF/PZAr,RIF/FQ/PZAr | HRZE"]),
              expression(eta["RIFr,RIF/FQr | STR"]),
              expression(eta["RIFr,RIF/PZAr | STR"]),
              expression(eta["RIF/FQr,RIF/FQ/PZAr | STR"]),
              expression(eta["RIF/PZAr,RIF/FQ/PZAr | STR"]), 
              expression(beta[0]))

# Save plots as pdf
pdf("Plots/histpanel.pdf")
par(mfrow=c(8, 4), mai=rep(0.3, 4))
for(i in 1:31){
  p1 <- hist(inputall[, i], breaks=50, plot=F)
  p2 <- hist(inputsel[, i], breaks=50, plot=F)
  plot(p2, col=rgb(1, 0, 0, 1), border="darkred", main=varnames[i], 
       freq=F, xlim=c(min(inputall[, i]), max(inputall[, i])))
  plot(p1, col=rgb(1, 1, 1, 0), border="darkgray", add=T, main=NULL, freq=F)
}
dev.off()

pdf("Plots/histpanel_NP.pdf")
par(mfrow=c(8, 4), mai=rep(0.3, 4))
for(i in 1:31){
  p1 <- hist(inputallb[, i], breaks=50, plot=F)
  p2 <- hist(inputselb[, i], breaks=50, plot=F)
  plot(p2, col=rgb(1, 0, 0, 1), border="darkred", main=varnames[i], 
       freq=F, xlim=c(min(inputallb[, i]), max(inputallb[, i])))
  plot(p1, col=rgb(1, 1, 1, 0), border="darkgray", add=T, main=NULL, freq=F)
}
dev.off()

# Quantitative assessment of strength of conditioning for each sampled input
# Kolmogorov-Smirnov test for difference in ECDFs
# Make a table summarizing the sampled inputs and selected inputs
sampletab <- matrix(nrow=31, ncol=10)
colnames(sampletab) <- c("Lower sampling bound", "Upper sampling bound",
                         "Median (sampled)", "Median (selected",
                         "25 %ile (sampled)", "75 %ile (sampled)",
                         "25 %ile (selected)", "75 %ile (selected)",
                         "KS D statistic", "KS p-value")
rownames(sampletab) <- varnames
sampletabb <- sampletab

for (i in 1:31) {
  sampletab[i, 1:2] <- c(min(inputall[, i]), max(inputall[, i]))
  sampletab[i, 3:4] <- c(median(inputall[, i]), median(inputsel[, i]))
  sampletab[i, 5:6] <- quantile(inputall[, i], c(0.25, 0.75))
  sampletab[i, 7:8] <- quantile(inputsel[, i], c(0.25, 0.75))
  tmp <- ks.test(inputall[, i], inputsel[, i])
  sampletab[i, 9:10] <- c(tmp[[1]], tmp[[2]])
  remove(tmp)
}
for (i in 1:31) {
  sampletabb[i, 1:2] <- c(min(inputallb[, i]), max(inputallb[, i]))
  sampletabb[i, 3:4] <- c(median(inputallb[, i]), median(inputselb[, i]))
  sampletabb[i, 5:6] <- quantile(inputallb[, i], c(0.25, 0.75))
  sampletabb[i, 7:8] <- quantile(inputselb[, i], c(0.25, 0.75))
  tmp <- ks.test(inputall[, i], inputselb[, i])
  sampletabb[i, 9:10] <- c(tmp[[1]], tmp[[2]])
  remove(tmp)
}
head(sampletab)
head(sampletabb)

# Write to csv
write.csv(sampletab, "Tables/sampletab.csv")
write.csv(sampletabb, "Tables/sampletab_NP.csv")







# IV: OUTCOMES (PRE-XDR, MDR, % vs. RAW PREVALENCE)----------------------------

# Drug resistance outcomes 
# Create 3D array of drug resistance outcomes over time, by category
drout1 <- array(dim=c(selsamp, 50, 6))
dimnames(drout1)[[2]] <- seq(1, 50, 1)
dimnames(drout1)[[3]] <- c("prex all", "prex rerx", "prex new",
                           "mdr all", "mdr rerx", "mdr new")

for (i in 2:SAtot) {
  assign(paste("drout", i, sep=""), drout1)
}

for(i in 1:SAtot){
  for (yr in 1:50) {
    eval(parse(
      text=paste("tmp<-data.frame(DSTlhsDataSA", i,
                 "[1:selsamp,(62+yr),-1])", sep="")))
    
    # Total pop w/ active TB
    totpop <- rowSums(tmp[,  c(3, 9, 15, 21, 27, 33, 39, 45, 6, 12, 18, 24, 30, 36, 42, 48)])
    # Total pop retreatment cases
    totrerx <- rowSums(tmp[, c(6, 12, 18, 24, 30, 36, 42, 48)])
    # Total new cases
    totnew <- rowSums(tmp[, c(3, 9, 15, 21, 27, 33, 39, 45)])
    
    # Total pre-XDR
    prexdr <- rowSums(tmp[, c(27, 30, 45, 48)])
    # Retreatment pre-XDR cases
    prexdrrerx <- rowSums(tmp[, c(30, 48)])
    # New pre-XDR cases
    prexdrnew <- rowSums(tmp[, c(27, 45)])
    
    # Total MDR
    mdr <- rowSums(tmp[, c(9, 27, 33, 45, 12, 30, 36, 48)])
    # Retreatment MDR cases
    mdrrerx <- rowSums(tmp[, c(12, 30, 36, 48)])
    # New MDR cases
    mdrnew <- rowSums(tmp[, c(9, 27, 33, 45)])
    
    eval(parse(text=paste("drout", i, "[,yr,1]<-prexdr", sep="")))
    eval(parse(text=paste("drout", i, "[,yr,2]<-prexdrrerx", sep="")))
    eval(parse(text=paste("drout", i, "[,yr,3]<-prexdrnew", sep="")))
    eval(parse(text=paste("drout", i, "[,yr,4]<-mdr", sep="")))
    eval(parse(text=paste("drout", i, "[,yr,5]<-mdrrerx", sep="")))
    eval(parse(text=paste("drout", i, "[,yr,6]<-mdrnew", sep="")))
    
    remove(totpop, totrerx, totnew, prexdr, prexdrrerx, prexdrnew, mdr, mdrrerx, mdrnew)
  }  
}

# Write to csv
write.csv(drout1, "Tables/outprev.csv")
write.csv(drout2, "Tables/outprev_AllPZA.csv")
write.csv(drout3, "Tables/outprev_F-RF.csv")
write.csv(drout4, "Tables/outprev_R-RF.csv")
write.csv(drout5, "Tables/outprev_P-RP.csv")
write.csv(drout6, "Tables/outprev_R-RP.csv")
write.csv(drout7, "Tables/outprev_FP-RFP.csv")
write.csv(drout8, "Tables/outprev_RP-RFP.csv")
write.csv(drout9, "Tables/outprev_R-RP-RFP.csv")
write.csv(drout10, "Tables/outprev_P-RP-RFP.csv")

# Distribution of outcomes at 20 years
# Pre-XDR
for (i in 1:SAtot) {
  eval(parse(text=paste("prexdr20", i, "<-drout", i, "[,20,1:3]", sep="")))
}
# MDR
for (i in 1:SAtot) {
  eval(parse(text=paste("mdr20", i, "<-drout", i, "[,20,4:6]", sep="")))
}

# Summary of all SA results
drSA <- array(dim=c(SAtot, 6, 3))
dimnames(drSA)[[1]] <- c("Baseline", "All PZA", "F->RF", "R->RF", "P->RP", "R->RP",
                         "FP->RFP", "RP->RFP", "R->RP->RFP", "P->RP->RFP")
dimnames(drSA)[[2]] <- c("PreXDR all", "PreXDR rerx", "PreXDR new",
                         "MDR all", "MDR rerx", "MDR new")
dimnames(drSA)[[3]] <- c("25th %ile", "Median", "75th %ile")

tmpfunc <- function(x) {
  tmp <- quantile(x[, 1], c(0.25, 0.5, 0.75))
  for (i in 2:6) {
    tmp <- rbind(tmp, quantile(x[, i], c(0.25, 0.5, 0.75)))
  }
  return(tmp)
}

for (i in 1:SAtot) {
  eval(parse(
    text=paste("drSA[", i, ",,]<-tmpfunc(drout", i, "[,20,])", sep="")))
}

write.csv(drSA, "Tables/outsum.csv")

tmp<-data.frame(DSTlhsDataSA1[, 62, -1])
# Total pop w/ active TB
totpop <- rowSums(tmp[,  c(3, 9, 15, 21, 27, 33, 39, 45, 6, 12, 18, 24, 30, 36, 42, 48)])
# Total pre-XDR
prexdr <- rowSums(tmp[, c(27, 30, 45, 48)])
summary(prexdr/totpop)


write.csv(drSA, "Tables/outsum15.csv")

# Incidence over time
incsum <- matrix(NA, nrow=6, ncol=6)
rownames(incsum) <- c("All", "MDR", "PreXDR", "All, SA", "MDR, SA", "PreXDR, SA")
colnames(incsum) <- c("Min", "25", "Med", "Mean", "75", "Max")
incsum[1, ] <- summary(DSTlhsDataSA1[, 82, 50])
incsum[2, ] <- summary(DSTlhsDataSA1[, 82, 51])
incsum[3, ] <- summary(DSTlhsDataSA1[, 82, 52])
incsum[4, ] <- summary(DSTlhsDataSA2[, 82, 50])
incsum[5, ] <- summary(DSTlhsDataSA2[, 82, 51])
incsum[6, ] <- summary(DSTlhsDataSA2[, 82, 52])
write.csv(incsum, "Tables/incsum.csv")








# V: TRAJECTORY PLOTS----------------------------------------------------------

## May need to clear all plots before running this...
# Overall TB prevalence
# All trajectories up to 2015
PlotDet(1, 1, DSTlhsData, 32, 62, DSTlhsDataSA1, axisyr=2015, ymax=250,
        savepath="Plots/", filename="TBAll_prev.png", save=F)
# Selected trajectories up to 2015
PlotDet(1, 1, DSTlhsDataSA1, 32, 62, col1="gray", axisyr=2015, ymax=250,
        savepath="Plots/", filename="TBSel_prev.png", save=F)
# Selected trajectories to 20 years
PlotDet(1, 1, DSTlhsDataSA1, 32, 82, col1="gray", ymax=350,
        savepath="Plots/", filename="TBSel20_prev.png", save=F)


# Overall TB incidence
# All trajectories up to 2015
PlotDet(3, 1, DSTlhsData, 32, 62, DSTlhsDataSA1, axisyr=2015, ymax=300,
        savepath="Plots/", filename="TBAll_pct.png", save=F)
# Selected trajectories to 20 years
PlotDet(3, 1, DSTlhsDataSA1, 32, 82, col1="gray", ymax=500,
        savepath="Plots/", filename="TBSel_pct.png", save=F)

# MDR TB incidence
# All trajectories up to 2015
PlotDet(3, 2, DSTlhsData, 32, 62, DSTlhsDataSA1, axisyr=2015, ymax=300,
        savepath="Plots/", filename="TBAll_pct.png", save=F)
# Selected trajectories to 20 years
PlotDet(3, 2, DSTlhsDataSA1, 32, 82, col1="gray", ymax=500,
        savepath="Plots/", filename="TBSel_pct.png", save=F)

# XDR TB incidence
# All trajectories up to 2015
PlotDet(3, 3, DSTlhsData, 32, 62, DSTlhsDataSA1, axisyr=2015, ymax=20,
        savepath="Plots/", filename="TBAll_pct.png", save=F)
# Selected trajectories to 20 years
PlotDet(3, 3, DSTlhsDataSA1, 32, 82, col1="gray", ymax=10,
        savepath="Plots/", filename="TBSel_pct.png", save=F)


# VI: Correlation between overall, MDR, and XDR TB incidence/prevalence
tmp1 <- tmp2 <- tmp3 <- NULL
for (i in 1:1015) {
    tmp1 <- c(tmp1, (DSTlhsDataSA1[i, 62:82, 50] - DSTlhsDataSA1[i, 62:82, 51]))
    tmp2 <- c(tmp2, (DSTlhsDataSA1[i, 62:82, 51] - DSTlhsDataSA1[i, 62:82, 52]))
    tmp3 <- c(tmp3, DSTlhsDataSA1[i, 62:82, 52])
}

# Plot labels
lab1 <- "DS TB Incidence"
lab2 <- "MDR TB Incidence"
lab3 <- "Pre-XDR TB Incidence"

# Compute r
r1 <- round(cor(tmp1, tmp2), 2)
r2 <- round(cor(tmp1, tmp3), 2)
r3 <- round(cor(tmp2, tmp3), 2)

# Plot
plot(tmp1, tmp2, xlab=lab1, ylab=lab2, pch=23, cex=0.5)
mtext(paste("r=",r1),side=4,line=0.5, las=1)
plot(tmp1, tmp2, xlim=c(150,250), ylim=c(0,50), xlab=lab1, ylab=lab2, pch=23, cex=0.5)

plot(tmp1, tmp3, xlab=lab1, ylab=lab3, pch=23, cex=0.5)
mtext(paste("r=",r2),side=4,line=0.5, las=1)
plot(tmp1, tmp3, xlim=c(150,250), ylim=c(0,5), xlab=lab1, ylab=lab3, pch=23, cex=0.5)

plot(tmp2, tmp3, xlab=lab2, ylab=lab3, pch=23, cex=0.5)
mtext(paste("r=",r3),side=4,line=0.5, las=1)
plot(tmp2, tmp3, xlim=c(0,50), ylim=c(0,5), xlab=lab2, ylab=lab3, pch=23, cex=0.5)

# Cross-correlation of time series
tmp1 <- DSTlhsDataSA1[, 62:82, 50] - DSTlhsDataSA1[, 62:82, 51]
tmp2 <- DSTlhsDataSA1[, 62:82, 51] - DSTlhsDataSA1[, 62:82, 52]
tmp3 <- DSTlhsDataSA1[, 62:82, 52]

val1 <- val2 <- val3 <- NULL

for (i in 1: 1015) {
  t1 <- ts(tmp1 [i, ])
  t2 <- ts(tmp2 [i, ])
  t3 <- ts(tmp3 [i, ])
  
  c1 <- ccf(t1, t2, plot=F)
  c2 <- ccf(t1, t3, plot=F)
  c3 <- ccf(t2, t3, plot=F)
  
  val1 <- rbind(val1, c(c1$acf[c1$lag == 0], c1$acf[c1$lag == 5], c1$acf[c1$lag == 10])) 
  val2 <- rbind(val2, c(c2$acf[c2$lag == 0], c2$acf[c2$lag == 5], c2$acf[c2$lag == 10]))
  val3 <- rbind(val3, c(c3$acf[c3$lag == 0], c3$acf[c3$lag == 5], c3$acf[c3$lag == 10]))
}  

# Median [IQR] of cross-correlation values give a sense
# of whether theres is a *consistent* pattern of cross-correlation
summary(val1)
summary(val2)
summary(val3)





