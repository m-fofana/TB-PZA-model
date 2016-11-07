# Mariam Fofana
# January 2016
# Modeling analysis of drivers of TB drug resistance trajectories
# Data analysis code
# Impact of RIF DST w/ standardized treatment

# Change working directory as necessary
setwd("C:/Users/mariam/Dropbox/TB DST Model")

# Load pacakges and source necessary files
if(!require(abind)){
  install.packages("abind")
  library(abind)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(fmsb)){
  install.packages("fmsb")
  library(fmsb)
}
source("R code/AnalysisFxnsDet.R")


# I: LOAD & FORMAT DATA--------------------------------------------------------

# Load data from formatted files
# Control (DST levels unchanged)
dataSA0 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestNone.csv"))

# RIF DST for Rep vs. Fail
dataSA1 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRep.csv")) 
dataSA2 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestFail.csv")) 

# RIF DST for varying proportion of retreatment (Rep and Fail combined)
dataSA3 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRerx50.csv"))
dataSA4 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRerx75.csv"))
dataSA5 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRerx100.csv"))

# Add RIF DST for 25, 50 75, 100% new cases
dataSA6 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRN25.csv"))
dataSA7 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRN50.csv"))
dataSA8 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRN75.csv"))
dataSA9 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestRN100.csv"))

# Improve treatment outcomes associated w/ FQ resistance
dataSA10 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestFQrxR.csv"))
dataSA11 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestFQrxRN.csv"))

# Block acquisition of FQ resistance on STR
dataSA12 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestFQacqR.csv"))
dataSA13 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestFQacqRN.csv"))

# Lower sensitivity 98->93% (may not be enough of a drop to see a big effect...)
dataSA14 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestSensR.csv"))
dataSA15 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestSensRN.csv"))

# Lower specificity 98->93% (may not be enough of a drop to see a big effect...)
dataSA16 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestSpecR.csv"))
dataSA17 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestSpecRN.csv"))

# Increased treatment completion 77->94%
dataSA18 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestCompR.csv"))
dataSA19 <- as.matrix(read.csv("Most recent output/Output/Prelim/DataTestCompRN.csv"))

# Total number of SA scenarios 
SAtot <- 19

# Check dimensions (if csv conversion has added extra column, remove it)
for (i in 0:SAtot) {
  eval(parse(text=
    paste(
      "if(ncol(dataSA", i, ")>6655){dataSA", i, "<-dataSA", i,"[,-1]}", sep=""
    )
  ))
}

# Reconstruct trajectory data into 3D arrays
for (i in 0:SAtot) {
	eval(parse(
		text=paste("DSTlhsDataSA", i, "<-Convert3D(dataSA", i, ",55,121)", sep="")))
}

selsamp <- nrow(dataSA1)


# II: OUTCOMES (PRE-XDR, MDR, % vs. RAW PREVALENCE)----------------------------

# Drug resistance outcomes 
# Create 3D array of drug resistance outcomes over time, by category
drout0 <- array(dim=c(selsamp, 21, 9))
dimnames(drout0)[[2]] <- seq(1, 21, 1)
dimnames(drout0)[[3]] <- c("prex all", "prex rerx", "prex new",
													 "mdr all", "mdr rerx", "mdr new",
													 "rp all", "rprerx", "rp new") # MDR +/- PZA (no FQ) res
propout0 <- drout0
for (i in 1:SAtot) {
	assign(paste("drout", i, sep=""), drout0)
  assign(paste("propout", i, sep=""), drout0)
  
}

for(i in 0:SAtot){
	for (yr in 1:21) {
		 eval(parse(
      text=paste("tmp<-data.frame(DSTlhsDataSA", i,
                          "[1:selsamp,(61+yr),-1])", sep="")))
	
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

    eval(parse(text=paste("propout", i, "[,yr,1]<-prexdr/totpop", sep="")))
    eval(parse(text=paste("propout", i, "[,yr,2]<-prexdrrerx/totrerx", sep="")))
    eval(parse(text=paste("propout", i, "[,yr,3]<-prexdrnew/totnew", sep="")))
    eval(parse(text=paste("propout", i, "[,yr,4]<-mdr/totpop", sep="")))
    eval(parse(text=paste("propout", i, "[,yr,5]<-mdrrerx/totrerx", sep="")))
    eval(parse(text=paste("propout", i, "[,yr,6]<-mdrnew/totnew", sep="")))
		
		remove(totpop, totrerx, totnew, prexdr, prexdrrerx, prexdrnew, mdr, mdrrerx, mdrnew)
	}  
}

# Write to csv
# write.csv(drout2, "Tables/outprev_DSTRifNew.csv")
# write.csv(drout3, "Tables/outprev_DSTRifRep.csv")
# write.csv(drout4, "Tables/outprev_DSTRifFail.csv")

# Distribution of outcomes at 20 years
# Pre-XDR
for (i in 0:SAtot) {
  eval(parse(text=paste("prexdr20", i, "<-drout", i, "[,21,1:3]", sep="")))
  eval(parse(text=paste("prexdrb20", i, "<-propout", i, "[,21,1:3]", sep="")))
}
# MDR
for (i in 0:SAtot) {
  eval(parse(text=paste("mdr20", i, "<-drout", i, "[,21,4:6]", sep="")))
  eval(parse(text=paste("mdrb20", i, "<-propout", i, "[,21,4:6]", sep="")))
}

# Summary of all SA results
drSA <- array(dim=c((SAtot+1), 6, 3))
dimnames(drSA)[[1]] <- c("Ctrl", "Baseline", 
                         "New", "Rep", "Fail",
                         "New25", "New50", "New75",
                         "Rerx50", "Rerx75", "Rerx100",
                         "All50", "All75", "All100")
dimnames(drSA)[[2]] <- c("PreXDR all", "PreXDR rerx", "PreXDR new",
												 "MDR all", "MDR rerx", "MDR new")
dimnames(drSA)[[3]] <- c("25th %ile", "Median", "75th %ile")
drSA2 <- drSA

tmpfunc <- function(x) {
  tmp <- quantile(x[, 1], c(0.25, 0.5, 0.75))
  for (i in 2:6) {
    tmp <- rbind(tmp, quantile(x[, i], c(0.25, 0.5, 0.75)))
  }
  return(tmp)
}

for (i in 0:SAtot) {
	eval(parse(
		text=paste("drSA[", (i+1), ",,]<-tmpfunc(drout", i, "[,21,])", sep="")))
  eval(parse(
    text=paste("drSA2[", (i+1), ",,]<-tmpfunc(propout", i, "[,21,])", sep="")))  
}

# write.csv(drSA, "Most recent output/Tables/DSTsum.csv")
# write.csv(drSA2, "Most recent output/Tables/DSTpropsum.csv")

## Results I---------------
# Prioritizing failure cases is better than prioritizing other retreatment cases
# Plot cumulative mortality
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,53])
tmp[,2]<-rowSums(DSTlhsDataSA2[,62:82,53])
tmp[,3]<-rowSums(DSTlhsDataSA1[,62:82,53])
tmp[,4]<-rowSums(DSTlhsDataSA5[,62:82,53])
boxplot(tmp, xlab="DST coverage", ylab="TB mortality (per 100K)", 
        main="Cumulative mortality at 20 years", xaxt="n")
mtext(at=c(1,2,3,4), side=1, c("Baseline","Failure","Other retreatment", "All retreatment"))

# Plot cumulative # treated
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA2[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA1[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA5[,62:82,54])
boxplot(tmp, xlab="DST coverage", ylab="# treated (per 100K)", 
        main="Cumulative # treated for MDR at 20 years", xaxt="n", log="y")
mtext(at=c(1,2,3,4), side=1, c("Baseline","Failure","Other retreatment", "All retreatment"))

# Plot MDR prevalence at 20 years
tmp<-cbind(drout0[,21,4], drout2[,21,4], drout1[,21,4])
boxplot(tmp, xlab="DST coverage", ylab="MDR prevalence (per 100K)", 
        main="MDR at 20 years",xaxt="n", log="y")
mtext(at=c(1,2,3), side=1, c("Baseline","Failure","Other retreatment"))

# Plot PreXDR prevalence at 20 years
tmp<-cbind(drout0[,21,1], drout2[,21,1], drout1[,21,1])
boxplot(tmp, xlab="DST coverage", ylab="Pre-XDR prevalence (per 100K)", 
        main="Pre-XDR at 20 years",xaxt="n", log="y")
mtext(at=c(1,2,3), side=1, c("Baseline","Failure","Other retreatment"))

# Plot Number needed to treat to avert 1 MDR case
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-rowSums(DSTlhsDataSA2[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA1[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA5[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])

tmp[,1]<-tmp[,1]/(drout2[,21,4]-drout0[,21,4])
tmp[,2]<-tmp[,2]/(drout1[,21,4]-drout0[,21,4])
tmp[,3]<-tmp[,3]/(drout5[,21,4]-drout0[,21,4])

boxplot(-tmp, xlab="DST coverage", ylab="Addl treatment per MDR case averted", 
        main="Number needed to treat",xaxt="n")
mtext(at=c(1,2,3), side=1, c("Failure","Other retreatment","All retreatment"))

# Plot # averted MDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=2)
tmp[,1]<-(drout2[,21,4]-drout0[,21,4])
tmp[,2]<-(drout1[,21,4]-drout0[,21,4])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted MDR cases",xaxt="n", log="y")
mtext(at=c(1,2), side=1, c("Failure","Other retreatment"))

# Plot percent reduction in MDR at 20 years
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout2[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,2]<-(drout1[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,3]<-(drout5[,21,4]-drout0[,21,4])/drout0[,21,4]*100

boxplot(-tmp, xlab="DST coverage", main="Reduction in MDR prevalence at 20 yrs", 
        ylab="% reduction",xaxt="n")
mtext(at=c(1,2,3), side=1, c("Failure","Other retreatment","All retreatment"))

## Results II---------------------
# Expanding MDR treatment access among retreatment cases
# Plot cumulative mortality
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,53])
tmp[,2]<-rowSums(DSTlhsDataSA3[,62:82,53])
tmp[,3]<-rowSums(DSTlhsDataSA4[,62:82,53])
tmp[,4]<-rowSums(DSTlhsDataSA5[,62:82,53])
boxplot(tmp, xlab="DST coverage, retreatment", ylab="TB mortality (per 100K)", 
        main="Cumulative mortality at 20 years", xaxt="n")
mtext(at=c(1,2,3,4), side=1, c("Baseline","50%","75%","100%"))

# Plot cumulative # treated
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA3[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA4[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA5[,62:82,54])
boxplot(tmp, xlab="DST coverage, retreatment", ylab="# treated (per 100K)", 
        main="Cumulative # treated for MDR at 20 years", xaxt="n", log="y")
mtext(at=c(1,2,3,4), side=1, c("Baseline","50%","75%","100%"))

# Plot MDR prevalence at 20 years
tmp<-cbind(drout0[,21,4], drout3[,21,4], drout4[,21,4], drout5[,21,4])
boxplot(tmp, xlab="DST coverage, retreatment", ylab="MDR prevalence (per 100K)", 
        main="MDR at 20 years",xaxt="n", log="y")
mtext(at=c(1,2,3,4), side=1, c("Baseline","50%","75%","100%"))

# Plot PreXDR prevalence at 20 years
tmp<-cbind(drout0[,21,1], drout3[,21,1], drout4[,21,1], drout5[,21,1])
boxplot(tmp, xlab="DST coverage, retreatment", ylab="Pre-XDR prevalence (per 100K)", 
        main="Pre-XDR at 20 years",xaxt="n", log="y")
mtext(at=c(1,2,3,4), side=1, c("Baseline","50%","75%","100%"))

# Plot Number needed to treat to avert 1 MDR case
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-rowSums(DSTlhsDataSA3[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA4[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA5[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,1]<-tmp[,1]/(drout3[,21,4]-drout0[,21,4])
tmp[,2]<-tmp[,2]/(drout4[,21,4]-drout0[,21,4])
tmp[,3]<-tmp[,3]/(drout5[,21,4]-drout0[,21,4])
boxplot(-tmp, xlab="DST coverage, retreatment", ylab="Addl treatment per MDR case averted", 
        main="Number needed to treat",xaxt="n")
mtext(at=c(1,2,3), side=1, c("50%","75%","100%"))

# Plot # averted MDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout3[,21,4]-drout0[,21,4])
tmp[,2]<-(drout4[,21,4]-drout0[,21,4])
tmp[,3]<-(drout5[,21,4]-drout0[,21,4])
boxplot(-tmp, xlab="DST coverage, retreatment", ylab="Averted cases", 
        main="Averted MDR cases",xaxt="n", log="y")
mtext(at=c(1,2,3), side=1, c("50%","75%","100%"))

# Plot percent reduction in MDR at 20 years
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout3[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,2]<-(drout4[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,3]<-(drout5[,21,4]-drout0[,21,4])/drout0[,21,4]*100
boxplot(-tmp, xlab="DST coverage, retreatment", main="Reduction in MDR prevalence at 20 yrs", 
        ylab="% reduction",xaxt="n")
mtext(at=c(1,2,3), side=1, c("50%","75%","100%"))

## Results III----------------
# Expansion of STR coverage to new cases might have perverse effects
# Plot cumulative mortality
tmp<-matrix(nrow=selsamp, ncol=6)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,53])
tmp[,2]<-rowSums(DSTlhsDataSA5[,62:82,53])
tmp[,3]<-rowSums(DSTlhsDataSA6[,62:82,53])
tmp[,4]<-rowSums(DSTlhsDataSA7[,62:82,53])
tmp[,5]<-rowSums(DSTlhsDataSA8[,62:82,53])
tmp[,6]<-rowSums(DSTlhsDataSA9[,62:82,53])

boxplot(tmp, xlab="DST coverage", ylab="TB mortality (per 100K)", 
        main="Cumulative mortality at 20 years", xaxt="n")
mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
                        "+ 50% New", "+75% New", "+100% New"))

# Plot cumulative # treated
tmp<-matrix(nrow=selsamp, ncol=6)
tmp[,1]<-rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA5[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA6[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA7[,62:82,54])
tmp[,5]<-rowSums(DSTlhsDataSA8[,62:82,54])
tmp[,6]<-rowSums(DSTlhsDataSA9[,62:82,54])
boxplot(tmp, xlab="DST coverage", ylab="# treated (per 100K)", 
        main="Cumulative # treated for MDR at 20 years", xaxt="n", log="y")
mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Plot MDR prevalence at 20 years
tmp<-cbind(drout0[,21,4], drout5[,21,4], drout6[,21,4], drout7[,21,4],
           drout8[,21,4], drout9[,21,4])
boxplot(tmp, xlab="DST coverage", ylab="MDR prevalence (per 100K)", 
        main="MDR at 20 years",xaxt="n", log="y")
mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Plot PreXDR prevalence at 20 years
tmp<-cbind(drout0[,21,1], drout5[,21,1], drout6[,21,1], drout7[,21,1],
           drout8[,21,1], drout9[,21,1])
boxplot(tmp, xlab="DST coverage", ylab="Pre-XDR prevalence (per 100K)", 
        main="Pre-XDR at 20 years",xaxt="n", log="y")
mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Plot Number needed to treat to avert 1 MDR case
tmp<-matrix(nrow=selsamp, ncol=5)
tmp[,1]<-rowSums(DSTlhsDataSA5[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA6[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA7[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA8[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,5]<-rowSums(DSTlhsDataSA9[,62:82,54]) - rowSums(DSTlhsDataSA0[,62:82,54])
tmp[,1]<-tmp[,1]/(drout5[,21,4]-drout0[,21,4])
tmp[,2]<-tmp[,2]/(drout6[,21,4]-drout0[,21,4])
tmp[,3]<-tmp[,3]/(drout7[,21,4]-drout0[,21,4])
tmp[,4]<-tmp[,4]/(drout8[,21,4]-drout0[,21,4])
tmp[,5]<-tmp[,5]/(drout9[,21,4]-drout0[,21,4])
boxplot(-tmp, xlab="DST coverage", ylab="Addl treatment per MDR case averted", 
        main="Number needed to treat",xaxt="n", ylim=c(0,300))
mtext(at=1:5, cex=0.8,side=1, c("Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Plot # averted MDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=5)
tmp[,1]<-(drout5[,21,4]-drout0[,21,4])
tmp[,2]<-(drout6[,21,4]-drout0[,21,4])
tmp[,3]<-(drout7[,21,4]-drout0[,21,4])
tmp[,4]<-(drout8[,21,4]-drout0[,21,4])
tmp[,5]<-(drout9[,21,4]-drout0[,21,4])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted MDR cases vs. baseline",xaxt="n")
mtext(at=1:5, cex=0.8,side=1, c("Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout6[,21,4]-drout5[,21,4])
tmp[,2]<-(drout7[,21,4]-drout5[,21,4])
tmp[,3]<-(drout8[,21,4]-drout5[,21,4])
tmp[,4]<-(drout9[,21,4]-drout5[,21,4])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted MDR cases vs. Retreatment only",xaxt="n")
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=4)
tmp[1]<-sum(drout6[,21,4]>drout5[,21,4])/selsamp*100
tmp[2]<-sum(drout7[,21,4]>drout5[,21,4])/selsamp*100
tmp[3]<-sum(drout8[,21,4]>drout5[,21,4])/selsamp*100
tmp[4]<-sum(drout9[,21,4]>drout5[,21,4])/selsamp*100
barplot(tmp, xlab="DST coverage", ylab="% simuations with increased MDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion sims with MDR above thresholds
hitraj <- matrix(NA, 5, 3)
colnames(hitraj) <- thresh <- seq(10, 20, 5)
rownames(hitraj) <- c("Retreatment", "+25% New", 
                      "+50% New", "+75% New", "+100% new")

tmp <- drout5[, 21, 4]
for (i in 6:9) {
  eval(parse(text=paste("tmp<-cbind(tmp,drout", i, "[,21,4])", sep="")))
}

for (i in 1:5) {
  for (j in 1:3) {
    hitraj[i, j] <- length(which(tmp[, i] > thresh[j]))
  }
}
hitraj <- (hitraj / selsamp) * 100

hitraj2 <- t(hitraj)
colnames(hitraj2) <- c("Baseline", "PZA replaced")
rownames(hitraj2) <- c("1/100,000", "1.5/100,000", "2/100,000")
par(mai=rep(1.2, 4))
bp <- barplot(hitraj2, beside=T, ylab="% simulations above threshold", 
              ylim=c(0, 30), axis.lty=0, axisnames=F, las=1,
              col=c(rgb(0, 0, 0, 0.2), rgb(0, 0, 0, 0.5),	rgb(0, 0, 0, 0.7)))
text(bp, (hitraj2), round(hitraj2, 1), cex=0.8, pos=3) 
mtext(c("Retreatment", "+25% New", "+50% New", "+75% New", "+100% new"),
      side=1, at=c(2,6,10,14,18), line=0.2, cex=0.8)
segments(0.3, 0, 16, 0, lty=1, col=1)		
legend("topright", c("1/100,000", "1.5/100,000", "2/100,000"), 
       title="MDR threshold in 2035", bty="n", cex=0.8, 
       fill=c(rgb(0, 0, 0, 0.2), rgb(0, 0, 0, 0.5), rgb(0, 0, 0, 0.7)))

# Plot percent reduction in MDR at 20 years
tmp<-matrix(nrow=selsamp, ncol=5)
tmp[,1]<-(drout5[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,2]<-(drout6[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,3]<-(drout7[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,4]<-(drout8[,21,4]-drout0[,21,4])/drout0[,21,4]*100
tmp[,5]<-(drout9[,21,4]-drout0[,21,4])/drout0[,21,4]*100
boxplot(-tmp, xlab="DST coverage", main="Reduction in MDR prevalence at 20 yrs", 
        ylab="% reduction v. baseline",xaxt="n")
mtext(at=1:5, cex=0.8,side=1, c("Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Change in MDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout6[,21,4]-drout5[,21,4])/drout5[,21,4]*100
tmp[,2]<-(drout7[,21,4]-drout5[,21,4])/drout5[,21,4]*100
tmp[,3]<-(drout8[,21,4]-drout5[,21,4])/drout5[,21,4]*100
tmp[,4]<-(drout9[,21,4]-drout5[,21,4])/drout5[,21,4]*100
boxplot(tmp, xlab="DST coverage", main="Change in MDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-20,60))
abline(h=0, col=2)
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))




# Plot # averted PreXDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=5)
tmp[,1]<-(drout5[,21,1]-drout0[,21,1])
tmp[,2]<-(drout6[,21,1]-drout0[,21,1])
tmp[,3]<-(drout7[,21,1]-drout0[,21,1])
tmp[,4]<-(drout8[,21,1]-drout0[,21,1])
tmp[,5]<-(drout9[,21,1]-drout0[,21,1])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted PreXDR cases vs. baseline",xaxt="n")
mtext(at=1:5, cex=0.8,side=1, c("Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout6[,21,1]-drout5[,21,1])
tmp[,2]<-(drout7[,21,1]-drout5[,21,1])
tmp[,3]<-(drout8[,21,1]-drout5[,21,1])
tmp[,4]<-(drout9[,21,1]-drout5[,21,1])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted PreXDR cases vs. Retreatment only",xaxt="n")
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=4)
tmp[1]<-sum(drout6[,21,1]>drout5[,21,1])/selsamp*100
tmp[2]<-sum(drout7[,21,1]>drout5[,21,1])/selsamp*100
tmp[3]<-sum(drout8[,21,1]>drout5[,21,1])/selsamp*100
tmp[4]<-sum(drout9[,21,1]>drout5[,21,1])/selsamp*100
barplot(tmp, xlab="DST coverage", ylab="% simuations with increased PreXDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion sims with PreXDR above thresholds
hitraj <- matrix(NA, 5, 3)
colnames(hitraj) <- thresh <- seq(1, 2, 0.5)
rownames(hitraj) <- c("Retreatment", "+25% New", 
                      "+50% New", "+75% New", "+100% new")

tmp <- drout5[, 21, 1]
for (i in 6:9) {
  eval(parse(text=paste("tmp<-cbind(tmp,drout", i, "[,21,1])", sep="")))
}

for (i in 1:5) {
  for (j in 1:3) {
    hitraj[i, j] <- length(which(tmp[, i] > thresh[j]))
  }
}
hitraj <- (hitraj / selsamp) * 100

hitraj2 <- t(hitraj)
colnames(hitraj2) <- c("Baseline", "PZA replaced")
rownames(hitraj2) <- c("1/100,000", "1.5/100,000", "2/100,000")
par(mai=rep(1.2, 4))
bp <- barplot(hitraj2, beside=T, ylab="% simulations above threshold", 
              ylim=c(0, 100), axis.lty=0, axisnames=F, las=1,
              col=c(rgb(0, 0, 0, 0.2), rgb(0, 0, 0, 0.5),	rgb(0, 0, 0, 0.7)))
text(bp, (hitraj2), round(hitraj2, 1), cex=0.8, pos=3) 
mtext(c("Retreatment", "+25% New", "+50% New", "+75% New", "+100% new"),
      side=1, at=c(2,6,10,14,18), line=0.2, cex=0.8)
segments(0.3, 0, 16, 0, lty=1, col=1)		
legend("topright", c("1/100,000", "1.5/100,000", "2/100,000"), 
       title="PreXDR threshold in 2035", bty="n", cex=0.8, 
       fill=c(rgb(0, 0, 0, 0.2), rgb(0, 0, 0, 0.5), rgb(0, 0, 0, 0.7)))

# Plot percent reduction in PreXDR at 20 years
tmp<-matrix(nrow=selsamp, ncol=5)
tmp[,1]<-(drout5[,21,1]-drout0[,21,1])/drout0[,21,1]*100
tmp[,2]<-(drout6[,21,1]-drout0[,21,1])/drout0[,21,1]*100
tmp[,3]<-(drout7[,21,1]-drout0[,21,1])/drout0[,21,1]*100
tmp[,4]<-(drout8[,21,1]-drout0[,21,1])/drout0[,21,1]*100
tmp[,5]<-(drout9[,21,1]-drout0[,21,1])/drout0[,21,1]*100
boxplot(-tmp, xlab="DST coverage", main="Reduction in preXDR prevalence at 20 yrs", 
        ylab="% reduction v. baseline",xaxt="n")
mtext(at=1:5, cex=0.8,side=1, c("Retreatment","+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Change in PreXDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout6[,21,1]-drout5[,21,1])/drout5[,21,1]*100
tmp[,2]<-(drout7[,21,1]-drout5[,21,1])/drout5[,21,1]*100
tmp[,3]<-(drout8[,21,1]-drout5[,21,1])/drout5[,21,1]*100
tmp[,4]<-(drout9[,21,1]-drout5[,21,1])/drout5[,21,1]*100
boxplot(tmp, xlab="DST coverage", main="Change in PreXDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-50,100))
abline(h=0, col=2)
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))






# Change in each strain vs. baseline
sum0 <- array(dim=c(selsamp, 6))  # For distribution of values at year 20
lab <- c("DS", "RIFr",  "FQr", "PZAr", "RIF/FQr","RIF/PZAr", "FQ/PZAr", "RIF/FQ/PZAr")
colnames(sum1)<-lab
sum9<-sum5<-sum0
for (i in c(0,5,9)) {
  eval(parse(text=paste("tmp <- DSTlhsDataSA",i,sep="")))
  
  res0 <- rowSums(tmp[, 82, c(4, 7)])  # DS
  res1 <- rowSums(tmp[, 82, c(10, 13)])  # RIF
  res2 <- rowSums(tmp[, 82, c(16, 19)])  # FQ
  res3 <- rowSums(tmp[, 82, c(22, 25)])  # PZA
  res4 <- rowSums(tmp[, 82, c(28, 31)])  # RIF/FQ
  res5 <- rowSums(tmp[, 82, c(34, 37)])  # RIF/PZA
  res6 <- rowSums(tmp[, 82, c(40, 43)])  # FQ/PZA
  res7 <- rowSums(tmp[, 82, c(46, 49)])  # RIF/FQ/PZA

  eval(parse(text=paste("sum",i," <- cbind(res0,res1 ,res2, res3, res4, res5, res6,res7)",sep="")))
} 

boxplot((sum5-sum0)/sum0*100, main="% change in prevalence, retreatment",
        xaxt="n", ylim=c(-100,150))
mtext(side=1, at=1:8, lab, cex=0.6)
abline(h=0, col=2)

boxplot((sum9-sum0)/sum0*100, main="% change in prevalence, all",
        xaxt="n", ylim=c(-100,150))
mtext(side=1, at=1:8, lab, cex=0.6)
abline(h=0, col=2)


boxplot((sum5-sum0), main="change in prevalence, retreatment",
        xaxt="n", ylim=c(-10,10))
mtext(side=1, at=1:8, lab, cex=0.6)
abline(h=0, col=2)

boxplot((sum9-sum0), main="change in prevalence, all",
        xaxt="n", ylim=c(-10,10))
mtext(side=1, at=1:8, lab, cex=0.6)
abline(h=0, col=2)


## Results IV--------------------
# Effect of FQ

# Plot cumulative # treated
tmp<-matrix(nrow=selsamp, ncol=6)
tmp[,1]<-rowSums(DSTlhsDataSA5[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA9[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA10[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA11[,62:82,54])
tmp[,5]<-rowSums(DSTlhsDataSA12[,62:82,54])
tmp[,6]<-rowSums(DSTlhsDataSA13[,62:82,54])
boxplot(tmp, xlab="DST coverage", ylab="# treated (per 100K)", 
        main="Cumulative # treated for MDR at 20 years", xaxt="n", log="y")
# mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
#                                 "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=3)
tmp[1]<-sum(drout9[,21,4]>drout5[,21,4])/selsamp*100
tmp[2]<-sum(drout11[,21,4]>drout10[,21,4])/selsamp*100
tmp[3]<-sum(drout13[,21,4]>drout12[,21,4])/selsamp*100
barplot(tmp, xlab="DST coverage", ylab="% simulations with increased MDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:3, cex=0.8,side=1, c("Baseline", "Decreased resistance", "Improved outcomes"))


# Change in MDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout9[,21,4]-drout5[,21,4])/drout5[,21,4]*100
tmp[,2]<-(drout11[,21,4]-drout10[,21,4])/drout5[,21,4]*100
tmp[,3]<-(drout13[,21,4]-drout12[,21,4])/drout5[,21,4]*100
boxplot(tmp, xlab="DST coverage", main="Change in MDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-20,60))
abline(h=0, col=2)
mtext(at=1:3, cex=0.8,side=1, c("Baseline", "Decreased resistance", "Improved outcomes"))




# Plot # averted PreXDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout9[,21,1]-drout5[,21,1])
tmp[,2]<-(drout11[,21,1]-drout10[,21,1])
tmp[,3]<-(drout13[,21,1]-drout12[,21,1])
boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted PreXDR cases vs. Retreatment only",xaxt="n")
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=3)
tmp[1]<-sum(drout9[,21,1]>drout5[,21,1])/selsamp*100
tmp[2]<-sum(drout11[,21,1]>drout10[,21,1])/selsamp*100
tmp[3]<-sum(drout13[,21,1]>drout12[,21,1])/selsamp*100
barplot(tmp, xlab="DST coverage", ylab="% simuations with increased PreXDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:3, cex=0.8,side=1, c("Baseline", "Decreased resistance", "Improved outcomes"))



# Change in PreXDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=3)
tmp[,1]<-(drout9[,21,1]-drout5[,21,1])/drout5[,21,1]*100
tmp[,2]<-(drout11[,21,1]-drout10[,21,1])/drout5[,21,1]*100
tmp[,3]<-(drout13[,21,1]-drout12[,21,1])/drout5[,21,1]*100
boxplot(tmp, xlab="DST coverage", main="Change in PreXDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-50,100))
abline(h=0, col=2)
mtext(at=1:3, cex=0.8,side=1, c("Baseline", "Decreased resistance", "Improved outcomes"))


# Results V----------------------------
# Plot cumulative # treated
tmp<-matrix(nrow=selsamp, ncol=8)
tmp[,1]<-rowSums(DSTlhsDataSA5[,62:82,54])
tmp[,2]<-rowSums(DSTlhsDataSA9[,62:82,54])
tmp[,3]<-rowSums(DSTlhsDataSA14[,62:82,54])
tmp[,4]<-rowSums(DSTlhsDataSA15[,62:82,54])
tmp[,5]<-rowSums(DSTlhsDataSA16[,62:82,54])
tmp[,6]<-rowSums(DSTlhsDataSA17[,62:82,54])
tmp[,7]<-rowSums(DSTlhsDataSA18[,62:82,54])
tmp[,8]<-rowSums(DSTlhsDataSA19[,62:82,54])
boxplot(tmp, xlab="DST coverage", ylab="# treated (per 100K)", 
        main="Cumulative # treated for MDR at 20 years", xaxt="n", log="y")
# mtext(at=1:6, cex=0.8,side=1, c("Baseline","Retreatment","+ 25% new",
#                                 "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=4)
tmp[1]<-sum(drout9[,21,4]>drout5[,21,4])/selsamp*100
tmp[2]<-sum(drout15[,21,4]>drout14[,21,4])/selsamp*100
tmp[3]<-sum(drout17[,21,4]>drout16[,21,4])/selsamp*100
tmp[3]<-sum(drout19[,21,4]>drout18[,21,4])/selsamp*100

barplot(tmp, xlab="DST coverage", ylab="% simulations with increased MDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:3, cex=0.8,side=1, c("Baseline", "Decreased resistance", "Improved outcomes"))





# Change in MDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout9[,21,4]-drout5[,21,4])/drout5[,21,4]*100
tmp[,2]<-(drout15[,21,4]-drout14[,21,4])/drout5[,21,4]*100
tmp[,3]<-(drout17[,21,4]-drout16[,21,4])/drout5[,21,4]*100
tmp[,4]<-(drout19[,21,4]-drout18[,21,4])/drout5[,21,4]*100

boxplot(tmp, xlab="DST coverage", main="Change in MDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-20,60))
abline(h=0, col=2)
mtext(at=1:4, cex=0.8,side=1, c("Baseline",
                                "Lower sens", "Lower spec", "Improved completion"))



# Plot # averted PreXDR cases at 20 years
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout9[,21,1]-drout5[,21,1])
tmp[,2]<-(drout15[,21,1]-drout14[,21,1])
tmp[,3]<-(drout17[,21,1]-drout16[,21,1])
tmp[,4]<-(drout19[,21,1]-drout18[,21,1])

boxplot(-tmp, xlab="DST coverage", ylab="Averted cases", 
        main="Averted PreXDR cases vs. Retreatment only",xaxt="n")
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Proportion of sims in which treating new cases is worse than retreatment only
tmp<-vector(length=4)
tmp[1]<-sum(drout9[,21,1]>drout5[,21,1])/selsamp*100
tmp[2]<-sum(drout15[,21,1]>drout14[,21,1])/selsamp*100
tmp[3]<-sum(drout17[,21,1]>drout16[,21,1])/selsamp*100
tmp[4]<-sum(drout19[,21,1]>drout18[,21,1])/selsamp*100

barplot(tmp, xlab="DST coverage", ylab="% simuations with increased PreXDR", 
        main="Treating new patients vs. retreatment only",xaxt="n", ylim=c(0,100))
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))


# Change in PreXDR vs. retreatment only
tmp<-matrix(nrow=selsamp, ncol=4)
tmp[,1]<-(drout9[,21,1]-drout5[,21,1])/drout5[,21,1]*100
tmp[,2]<-(drout15[,21,1]-drout14[,21,1])/drout5[,21,1]*100
tmp[,3]<-(drout17[,21,1]-drout16[,21,1])/drout5[,21,1]*100
tmp[,4]<-(drout19[,21,1]-drout18[,21,1])/drout5[,21,1]*100
boxplot(tmp, xlab="DST coverage", main="Change in PreXDR prevalence at 20 yrs", 
        ylab="% change v. retreatment",xaxt="n", ylim=c(-50,100))
abline(h=0, col=2)
mtext(at=1:4, cex=0.8,side=1, c("+ 25% new",
                                "+ 50% New", "+75% New", "+100% New"))

# Sampled vars------------------

init <- as.matrix(read.csv("Most recent output/Output/InputSelBsl_110915.csv"))

# Check dimensions (if csv conversion has added extra column, remove it)
if(ncol(init) > 74) {init <- init[, -1]}

init <- init[, -c(24, 25, 28, 30, 35:49)] 

inputsel <- init[, 1:31]  # input variables 

  himdr <- rep(0, selsamp)
  himdr[which(drout9[,21,4]>drout5[,21,4])] <- 1
  mdrdif<-drout9[,,4]-drout5[,,4]
  
  varnames <- c("time1", "time2",
                "acqr/M->RM", "acqf", "acqp/M->MZ",
                "fit_r", "fit_f", "fit_p", "fit_rf", "fit_rp", "fit_fp", "fit_rfp",
                "rx_ds", "rx_r", "rx_p", "rx_rp", "rx_r2", "rx_rf2", "rx_rp2", "rx_rfp2",
                "resmult",
                "R->RM", "R->RZ/RM->RMZ", "Z->RZ","Z->MZ",
                "RZ->RMZ",
                "R->RM STR", "R->RZ STR", "RM->RMZ STR", "RZ->RMZ STR",
                "beta")
  
  for(i in 1:31){
    plot(inputsel[,i],mdrdif, xlab=varnames[i])
  }
  
  # Univariate regression
  # Create scaled values
  inputs2 <- matrix(ncol=31, nrow=selsamp)
  insd <- apply(inputsel, 2, sd)
  inave <- apply(inputsel, 2, mean)
  for (j in 1:31) {
    inputs2[, j] <- ((inputsel[, j] - inave[j]) / insd[j]) * 10
  }
  summary(inputs2)  
  
  # Compile input and output variables into dataframe
  regdata <- inputs2  # These are the sampled parameters (i.e., covariates)
  colnames(inputs2) <- varnames
  out1 <- himdr  

  tmpdata1 <- cbind(regdata, out1)
  tmpdata1 <- data.frame(tmpdata1)

  # Construct results table to compile regression results
  regprex <- matrix(nrow=31, ncol=9)  # Pre-XDR
  colnames(regprex) <- c("Uni coeff", "Uni SE", "Uni p",
                         "Multi coeff", "Multi SE", "Multi p",
                         "VIF coeff", "VIF SE", "VIF p")
  rownames(regprex) <- varnames

  # Univariate regression using scaled values
  for (i in 1:31) {
    tmp1 <- summary(glm(out1 ~ tmpdata1[, i],
                        family=binomial("logit"), data=tmpdata1))

    
    # Pre-XDR results
    regprex[i, 1] <- tmp1$coef[2, 1]  # estimate
    regprex[i, 2] <- tmp1$coef[2, 2]  # SE
    regprex[i, 3] <- tmp1$coef[2, 4]  # p-value

  }
  
  # Multivariate regression using scaled values, all covariates
  tmp1 <- summary(glm(out1 ~ .,
                      family=binomial("logit"), data=tmpdata1))

  
  regprex[, 4] <- tmp1$coef[-1, 1]  # estimate
  regprex[, 5] <- tmp1$coef[-1, 2]  # SE
  regprex[, 6] <- tmp1$coef[-1, 4]  # p-value
  

  # Multivar regression excluding collinear variables based on VIF, threshold=10
  tmpdata <- tmpdata1[, 1:31]
  vif10 <- vif_func(tmpdata, thresh=10, trace=T) 
  
  # Re-run multiple regression model with restricted set of variables
  form <- vif10[1]
  for (i in 2:length(vif10)) {
    form <- paste(form, "+", vif10[i], sep="")
  }
  
  form1 <- formula(paste("out1 ~ ", form))

  
  tmp1 <- summary(glm(form1,
                      family=binomial("logit"), data=tmpdata1))

  # Compile regression results into table
  regprex[, 7] <- tmp1$coef[-1, 1]  # estimate
  regprex[, 8] <- tmp1$coef[-1, 2]  # SE
  regprex[, 9] <- tmp1$coef[-1, 4]  # p-value

  
  # Write results to csv
#   write.csv(regprex, paste("Tables/Reg/Regprexdr", prev, "prev.csv", sep=""))
#   write.csv(regmdr, paste("Tables/Reg/Regmdr", prev2, "prev.csv", sep=""))

  varnames4 <- c("Time emerged, RIF resistance", "Time emerged, FQ resistance", 
                 "RIF resistance acquisition", "FQ resistance acquisition", 
                 "PZA resistance acquisition", "Fitness, RIFr", "Fitness, FQr", 
                 "Fitness, PZAr", "Fitness, RIF/FQr", "Fitness, RIF/PZAr", 
                 "Fitness, FQ/PZAr", "Fitness, RIF/FQ/PZAr", "1st-line cure, DS",
                 "1st-line cure, RIFr", "1st-line cure, PZAr", 
                 "1st-line cure, RIF/PZAr", "2nd-line cure, RIFr", 
                 "2nd-line cure, RIF/FQr", "2nd-line cure, RIF/PZAr", 
                 "2nd-line cure, RIF/FQ/PZAr", 
                 "RR amplification, retreatment", 
                 expression(paste(RIFr%->%RIF/FQr, " amplification")),
                 expression(paste(RIFr%->%RIF/PZAr, " amplification")),
                 expression(paste(PZAr%->%RIF/PZAr, " amplification")),
                 expression(paste(PZAr%->%FQ/PZAr, " amplification")), 
                 expression(paste(RIF/PZAr%->%RIF/FQ/PZAr, " amplification")),
                 expression(paste(RIF/FQr%->%RIF/FQ/PZAr, " amplification")),
                 expression(paste(RIFr%->%RIF/PZAr, " amplification")),
                 expression(paste(RIF/FQr%->%RIF/FQ/PZAr, " amplification")),
                 expression(paste(RIF/PZAr%->%RIF/FQ/PZAr, " amplification")),
                 "Transmission rate") 

# Ranking of regression coefficients
SARegp <- matrix(nrow=12, ncol=3)
SARegp[2, ] <- c("Param", "Coeff", "SE")

  # Read in values from saved csv
  tmp <- regprex

  # Remove non-sig p-values
  tmp2 <- which(tmp[, 9] > 0.05)
  if (sum(tmp2) > 0) {
    tmp <- tmp[-tmp2, ]
    tmpnames<-varnames4[-tmp2]
    }
  
  # Remove NA values
  tmp2 <- which(is.na(tmp[, 7]))
  if (sum(tmp2) > 0) {
    tmp <- tmp[-tmp2, ]
    tmpnames<-tmpnames[-tmp2]
    }
  
  # Order by absolute value
  tmp2 <- order(abs(tmp[, 7]), decreasing=T)
  tmp2 <- tmp2[1:10]
  tmp <- tmp[tmp2, ]
  tmpnames<-tmpnames[tmp2]
  
  # Fill in master table
  SARegp[3:12, 1] <- rownames(tmp)
  SARegp[3:12,2:3]<-as.numeric(tmp[,7:8])

  # Plot
  tmppt <- as.numeric(tmp[, 7])
  tmpse <- as.numeric(tmp[, 8])
  or <- data.frame(matrix(nrow=10, ncol=4))
  or[, 1] <- as.character(tmp[, 1])
  or[, 2] <- exp(tmppt)
  or[, 3] <- exp(tmppt - 1.96 * tmpse)
  or[, 4] <- exp(tmppt + 1.96 * tmpse)
  
 

  par(mai=c(1, 3, 1, 1))
  plot(or[, 2], seq(10, 1), pch=15, xaxp=c(0.5, 3, -1), log="x",
       ylab="", xlab="", bty="n", yaxt="n", tck=0, xlim=c(0.5, 2.5), cex=1.3,
       col=1, cex.axis=0.8)
  arrows(or[, 3], seq(10, 1), or[, 4], seq(10, 1),
         col=1, angle=0, length=0.1, code=3)
  mtext(c(1, 2), side=1, at=c(1, 2), line=1, cex=0.8)
  mtext(expression(paste("Odds ratio ", "("%+-%95, "% confidence interval)")),
        side=1, line=2.2, cex=0.8)
  segments(1, 0, 1, 12, lty=1, col=1, lwd=2)
  mtext(tmpnames, side=2, line=0, at=seq(10, 1), las=1, cex=0.7, padj=-0.1)
  

# PRCCs
  prexall <- matrix(NA, 31, 20)
  rownames(prexall) <- varnames
  colnames(prexall) <- seq(1, 20, 1)
  for (yr in 2:21) {
    tmpprex1 <- cbind(inputsel, mdrdif[, yr])

    for (i in 1:31) {
      tmp1 <- pcor.test(tmpprex1[, i], tmpprex1[, 32], tmpprex1[, -c(i, 32)], 
                        method="spearman")
      prexall[i, (yr-1)] <- tmp1$estimate
      remove(tmp1)
    }
    print(yr)  # Counter
  }  

  prcc <- t(prexall)
  cols <- c("deeppink4", "dodgerblue1")
  cex.lab <- 2
  sp <- .3
  
  lab.txt <- c(expression(t[1]), expression(t[2]), expression(eta["DS,RIFr | HRZE"]),
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
  prcc.order <- order(abs(prcc[20, ]), decreasing=T)
  
  # pdf(height=11, width=8, file="Plots/prcc_mdr.pdf")
  par(mar=c(5, 12, 8, 2))
  plot(range(0, 20), range(0, 10 * (2 + sp)), type="n", axes=F, xlab="", ylab="")
  gp <- 1
  for (jp in 1:10) {
    jp2 <- 11 - jp
    rect(0:19, gp + rep(0, times=20), 1:20, gp + prcc[, prcc.order[jp2]],
         col=cols[ceiling(prcc[, prcc.order[jp2]] + 1)])
    axis(side=2, at=c(gp - 1, gp, gp + 1), c(-1, 0, 1), tcl=-0.25, padj=0.5, 
         las=1, cex.axis=0.7, hadj=0.5)
    axis(side=4, at=c(gp - 1, gp, gp + 1), c(-1, 0, 1), tcl=-0.25, padj=0.5, 
         las=1, cex.axis=0.7, hadj=0.5)
    mtext(side=2, at=gp + 0, line=2, adj=1, las=1, lab.txt[prcc.order[jp2]], cex=1.3)
    text(10 - 2, gp + 1, signif(prcc[10, prcc.order[jp2]], 3))
    text(20 - 2, gp + 1, signif(prcc[20, prcc.order[jp2]], 3))
    gp <- gp + 2 + sp
  }
  axis(side=1, at=seq(0, 20, by=5), labels=seq(2015, 2035, 5))
  axis(side=3, at=seq(0, 20, by=5), labels=seq(2015, 2035, 5))
  
  mtext(side=1, "Year", padj=3, cex=1.3)
  mtext(side=3, "Year", padj=-3, cex=1.3)
  mtext(side=3, "Partial rank correlation coefficient, RIF resistance prevalence", 
        padj=-4.5, cex=1.5, at=-7, adj=0)
  mtext(side=3, at=-2, "Parameter", padj=-1, cex=1.3, adj=1)
  abline(v=c(10, 20), lty=2, col="grey65")
  # dev.off()
  

# IV: TRAJECTORY MEDIAN + ENVELOPE---------------------------------------------			
# Reformat data to show plot of median + IQR/median + range, overlaid
dat1 <- data.frame(t(apply(drout0[, , 4], 2, 
                     function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
dat2 <- data.frame(t(apply(drout5[, , 4], 2, 
                     function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
dat3 <- data.frame(t(apply(drout6[, , 4], 2, 
                     function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
dat4 <- data.frame(t(apply(drout9[, , 4], 2, 
                           function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))


datb1 <- data.frame(t(apply(propout0[, , 4]*100, 2, 
                           function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
datb2 <- data.frame(t(apply(propout5[, , 4]*100, 2, 
                           function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
datb3 <- data.frame(t(apply(propout6[, , 4]*100, 2, 
                           function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))
datb4 <- data.frame(t(apply(propout9[, , 4]*100, 2, 
                            function(x) {quantile(x, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))})))

dat1$yr <- dat2$yr <- dat3$yr <- dat4$yr <- 2015:2035
yr2 <- dat1$yr[c(0:21, 21:0)] 


# Plot with base graphics
plot(dat1$yr, dat1$X50., xlim=c(2015, 2035), xaxp=c(2015, 2035, 4),
     type="l", lwd=2, col=1, xlab="Year", ylab="MDR prevalence (per 100K)",
     ylim=c(0, 15), las=1)
polygon(yr2, c(dat1$X5., rev(dat1$X95.)), col=rgb(0.5, 0.5, 0.5, 0.3), border=F)
lines(dat2$yr, dat2$X50., col="blue", lwd=2)
polygon(yr2, c(dat2$X5., rev(dat2$X95.)), col=rgb(0, 0, 0.75, 0.3), border=F)
lines(dat4$yr, dat4$X50., col="red", lwd=2)
polygon(yr2, c(dat4$X5., rev(dat4$X95.)), col=rgb(0.75, 0, 0, 0.3), border=F)

# 
# 
# plot(dat1$yr, dat1$X50., xlim=c(2015, 2035), xaxp=c(2015, 2035, 4),
#      type="l", lwd=2, col=1, xlab="Year", ylab="MDR prevalence (per 100K)",
#      ylim=c(0, 30), las=1)
# polygon(yr2, c(dat1$X5., rev(dat1$X95.)), col=rgb(0.5, 0.5, 0.5, 0.3), border=F)
# lines(dat3$yr, dat3$X50., col="blue", lwd=2)
# polygon(yr2, c(dat3$X5., rev(dat3$X95.)), col=rgb(0, 0, 0.75, 0.3), border=F)
# legend("topleft", c("None", "RIF DST, all failure cases", "5th-95th percentile"), 
#        bty="n", col=c(1, "blue", 0), lwd=2, cex=0.8, title="DST scale-up")
# polygon(c(2014.75, 2015.3, 2015.3, 2014.75), c(23, 23, 25, 25), 
#         col=rgb(0, 0, 0.75, 0.3), border=F)
# polygon(c(2015.4, 2016, 2016, 2015.4), c(23, 23, 25, 25), col="gray", border=F)
# 
# Proportion MDR
plot(dat1$yr, datb1$X50., xlim=c(2015, 2035), xaxp=c(2015, 2035, 4),
     type="l", lwd=2, col=1, xlab="Year", ylab="% MDR",
     ylim=c(0, 15), las=1)
polygon(yr2, c(datb1$X25., rev(datb1$X75.)), col=rgb(0.5, 0.5, 0.5, 0.3), border=F)
lines(dat3$yr, datb2$X50., col="blue", lwd=2)
polygon(yr2, c(datb2$X25., rev(datb2$X75.)), col=rgb(0, 0, 0.5, 0.3), border=F)
lines(dat3$yr, datb4$X50., col="red", lwd=2)
polygon(yr2, c(datb4$X25., rev(datb4$X75.)), col=rgb(0.75, 0, 0, 0.3), border=F)



# plot(dat1$yr, datb1$X50., xlim=c(2015, 2035), xaxp=c(2015, 2035, 4),
#      type="l", lwd=2.5, col=1, xlab="Year", ylab="% MDR",
#      ylim=c(0, 12), las=1)
# polygon(yr2, c(datb1$X5., rev(datb1$X95.)), col=rgb(0.5, 0.5, 0.5, 0.25), border=F)
# lines(dat3$yr, datb3$X50., col="red", lwd=2.5)
# polygon(yr2, c(datb3$X5., rev(datb3$X95.)), col=rgb(1, 0, 0, 0.25), border=F)
# lines(dat4$yr, datb4$X50., col="blue", lwd=2.5)
# polygon(yr2, c(datb4$X5., rev(datb4$X95.)), col=rgb(0, 0, 1, 0.25), border=F)
# legend("topleft", c("None", "All failure cases", "All previously treated"), 
#        bty="n", col=c(1, "red", "blue"), lwd=2, cex=0.8, title="RIF DST scale-up")
# #polygon(c(2014.75, 2015.3, 2015.3, 2014.75), c(10.5, 10.5, 11.5, 11.5), 
# #        col=rgb(0, 0, 0.75, 0.3), border=F)
# #polygon(c(2015.4, 2016, 2016, 2015.4), c(10.5, 10.5, 11.5, 11.5), col="gray", border=F)
# mtext(side=2, line=1.9, las=3, cex=0.8,
#       text=expression(paste("(Median, ", 5^th-95^th, " percentiles)")))

