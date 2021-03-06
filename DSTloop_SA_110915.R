taskid <- Sys.getenv("SGE_TASK_ID")
taskid <- as.numeric(taskid)
print(paste("Running task # ", taskid, sep=""))

#Mariam Fofana
#November 2015
#Model of DST for new FQ-based first-line TB regimens 
#SA: Reduce/null out resistance amplification pathways

#Change working directory as necessary
setwd("~/DST2")
source("Rcode/DSTfxns_110915.R")

#Load required packages
library(abind)
library(deSolve)

#Load data from formatted files
init<-read.csv("Baseline/Init/InputSelBsl_110915.csv")
data<-read.csv("Baseline/Data/DataSelBsl_110915.csv")

#Separate contents of "Input" files
if(ncol(init)>74){init<-init[,-1]}
inputsel<-init[,1:50] #input variables

#Check dimensions (if csv conversion has added extra column, remove it)
if(ncol(data)>6292){data<-data[,-1]}

#Reconstruct trajectory data into 3D array
#Selected trajectories
DSTlhsDatasel<-data[,1:121]
for(j in 2:52){
  print(paste("j=",j))
  min<-(j-1)*121+1
  max<-j*121
  temp2<-data[,min:max]
  DSTlhsDatasel<-abind(DSTlhsDatasel,temp2,along=3)
  remove(temp2)
}

dim(DSTlhsDatasel)

stateinit<-DSTlhsDatasel[,62,2:49] #Set to [,1,2:49] if starting from beginning, [,62,2:49] if starting from 2015

#SET PARAMETERS

#DST availability and test characteristics
#CURRENTLY NOT ACCOUNTING FOR IMPERFECT SPECIFICITY
#Create empty 3x8 matrix of proportion getting each type of DST
#Rows are target pop (new cases, non-failure rerx, failure)
#Columns are type of DST received
#This is combo of lab-confirmed & presumptive dx of DR (hence why #s so high for retreatment cases)

dst_access<-array(0,dim=c(3,8,2))
dimnames(dst_access)[[2]]<-c("None","R","F","P","RF","RP","FP","RFP")
dimnames(dst_access)[[1]]<-c("New","Rep","Fail")
dimnames(dst_access)[[3]]<-c("Pre","Post") #pre vs post intro of FQ regimen

dst_access[,2,1]<-c(0.054,0.261,0.261) #At baseline assume that 15% of new cases and 26% rerx are getting tested for Rif

#Change access to DST after intro of FQ regimen by modifying lines below as necessary
dst_access[,2,2]<-dst_access[,2,1]

for(i in 1:3){
  dst_access[i,1,1]<-1-sum(dst_access[i,-1,1])
  dst_access[i,1,2]<-1-sum(dst_access[i,-1,2])
}

#Sensitivity
DSTr_sens<-0.98 #Sensitivity of Rif resistance diagnosis
DSTf_sens<-0.93 #Sensitivity of FQ resistance diagnosis
DSTp_sens<-0.80 #Sensitivity of PZA resistance diagnosis

DSTrf_sens<-DSTr_sens*DSTf_sens #Sensitivity of combined Rif+FQ resistance diagnosis
DSTrp_sens<-DSTr_sens*DSTp_sens #Sensitivity of combined Rif+PZA resistance diagnosis
DSTfp_sens<-DSTf_sens*DSTp_sens #Sensitivity of combined FQ+PZA resistance diagnosis
DSTx_sens<-DSTr_sens*DSTf_sens*DSTp_sens #Sensitivity of combined DST for 3 drugs

#Choice of treatment regimen
#Create 8x6x2 array
#Rows are dx'd DR state
#Columns are choice of treatment regimen
#3rd dimension is pre- vs post-intro of new regimen, rx-naive vs. retreatment
#Each cell = proportion getting a particular regimen
#FQrx = FQ-based 1st-line regimen
#STR = standardized MDR regimen (used w/o comprehensive DST)
#DST = DST-based individualized regimen

dxrx<-array(0,dim=c(8,4,2))
dimnames(dxrx)[[1]]<-c("DS","R","F","P","RF","RP","FP","RFP")
dimnames(dxrx)[[2]]<-c("HRZE","FQrx","STR","DST") 
dimnames(dxrx)[[3]]<-c("Pre","Post")

#At baseline:
#In "pre" era, only those dx'd w/ Rif resistance get STR regimen, all others get HRZE
#In "post" era, only those dx'd as DS get FQrx, all others get DST-based rx
dxrx[1,1,1]<-dxrx[1,2,2]<-1 #100% of DS-TB get HRZE in pre era, 100% get FQrx in post era
dxrx[2,3,1]<-dxrx[2,4,2]<-1 #100% of Rif-R TB get STR in pre era, 100% get DST-based Rx in post era
dxrx[3,1,1]<-dxrx[3,4,2]<-1 #100% of FQ-R TB get HRZE in pre era, DST-based rx in post
dxrx[4,1,1]<-dxrx[4,4,2]<-1 #100% of PZA-R get HRZE in pre era (no PZA dx), 100% get DST-based Rx in post era
dxrx[5,3,1]<-dxrx[5,4,2]<-1 #100% of Rif+FQ-R get STR in pre era, 100% get DST-based rx in post era
dxrx[6,3,1]<-dxrx[6,4,2]<-1 #100% of Rif+PZA-R get STR in pre era, 100% get DST-based rx in post era
dxrx[7,1,1]<-dxrx[7,4,2]<-1 #100% of FQ+PZA-R get HRZE in pre era, 100% get DST-based rx in post era
dxrx[8,3,1]<-dxrx[8,4,2]<-1 #100% of triple-res get STR in pre era, 100% get DST-based rx in post era


#To change timesteps
totaltime<-60 #Total # of years simulated (120 if starting from beginning, set to 60 if starting in 2015)
increment<-1 #Size of time steps, in years
intro<-121 #Time of introduction of FQ-based regimen (If this is >totaltime, then FQ-based regimen is never implemented)

#Treatment completion
comp1<-1-0.06 # % completing rx 1st-line (Source: WHO)
comp2<-1-0.17 # % completing empiric retreatment (Source: WHO)
comp3<-1-0.23 # % completing 2nd-line rx and STR (Source: Ahuja et al)
comp<-c(comp1,comp2,comp3)


#READ IN SAMPLED PARAMETERS FROM CSV FILE
#***Comment out this section if generating a new set of sampled parameters instead***
#To read in sampled parameters from csv file, uncomment this section 
#and modify file path & name as necessary

dstlhs1<-as.matrix(inputsel)
sampleno<-nrow(dstlhs1)


#GENERATE MATRIX TO COLLECT RELEVANT INPUT & OUTPUT DATA
dstlhs2<-matrix(ncol=(50+12+12), nrow=sampleno) #Matrix to hold 50 sampled values + 12 output var + 12 selection var
DSTlhsData<-array(0,dim=c(sampleno,121,52)) #3D array to hold # in each of 48 compartments + time (total 49 values)
#Plus 3 incidence values (overall, MDR, Pre-XDR)
#for each year simulated (totaltime + 1 for year 0) 
#and each individual simulation (sampleno)


#START SIMULATIONS!

#Create dummy variable to track # of simulations that meet selection crieria, reset to 0 
row<-0 

#BEGIN SIMULATION LOOP
# If no effect of PZA resistance on further amplification, change parameters
for(loop in 1:sampleno){
  if (taskid == 1) {
		dstlhs1[loop,26]<-dstlhs1[loop,3] #Set PZA->RIF/PZA same as DS->RIF
		dstlhs1[loop,27]<-dstlhs1[loop,4]# PZA->FQ/PZA, HRZE
		dstlhs1[loop,29]<-dstlhs1[loop,22]# RIF/PZA->3x, HRZE
		dstlhs1[loop,34]<-dstlhs1[loop,31]# RIF/PZA->3x, STR	
	}

  ##################
  #DEFINE PARAMETERS
  ##################
  
  #Times of introduction of resistance
  time1<-dstlhs1[loop,1]
  time2<-dstlhs1[loop,2]
  
  
  #Prob of monoresistance acquisition, for DS-TB
  acqr<-dstlhs1[loop,3]
  acqf<-dstlhs1[loop,4]
  acqp<-dstlhs1[loop,5]
  
  #"resistance multiplier"
  resmult<-dstlhs1[loop,21]
  
  #Acquisition of multiple resistance
  acqr_rf<-dstlhs1[loop,22] #Rif-->Rif+ FQ
  acqr_rp<-dstlhs1[loop,23] #Rif-->Rif+PZA
  acqf_rf<-dstlhs1[loop,24] #FQ-->Rif+FQ
  acqf_fp<-dstlhs1[loop,25] #FQ-->FQ+PZA
  acqp_rp<-dstlhs1[loop,26] #PZA-->Rif+PZA
  acqp_fp<-dstlhs1[loop,27] #PZA-->FQ+PZA
  acqrf_rfp<-dstlhs1[loop,28] #Rif+FQ-->triple resistance
  acqrp_rfp<-dstlhs1[loop,29] #Rif+PZA-->triple resistance
  acqfp_rfp<-dstlhs1[loop,26] #FQ+PZA-->triple resistance
  
  #Acquisition of multiple resistance on standardized retreatment regimen (STR)
  acqr_rf2<-dstlhs1[loop,31] #Rif-->Rif+FQ
  acqr_rp2<-dstlhs1[loop,32] #Rif-->Rif+PZA
  acqrf_rfp2<-dstlhs1[loop,33] #Rif+FQ-->triple resistance
  acqrp_rfp2<-dstlhs1[loop,34] #Rif+PZA-->triple resistance
  
	# Modify resistance acquisition for SAs
	# Set specific resistance acquisition values to 0
	if (taskid == 2) {
			acqf_rf<-0 #FQ-->Rif+FQ
	} else if (taskid == 3) {
		acqr_rf<-0 #Rif-->Rif+ FQ	  
	} else if (taskid == 4) {
		acqp_rp<-0 #PZA-->Rif+PZA
	} else if (taskid == 5) {
	  acqr_rp<-0 #Rif-->Rif+PZA
		acqr_rp2<-0 #Rif-->Rif+PZA
	} else if (taskid == 6) {
	  acqfp_rfp<-0 #FQ+PZA-->triple resistance
	} else if (taskid == 7) {
		acqrp_rfp<-0 #Rif+PZA-->triple resistance
		acqrp_rfp2<-0 #Rif+PZA-->triple resistance
	} else if (taskid == 8) {
	  acqr_rp<-0 #Rif-->Rif+PZA
		acqrp_rfp<-0 #Rif+PZA-->triple resistance
		acqrp_rfp2<-0 #Rif+PZA-->triple resistance
	} else if (taskid == 9) {
	  acqp_rp<-0 #PZA-->Rif+PZA
		acqrp_rfp<-0 #Rif+PZA-->triple resistance
		acqrp_rfp2<-0 #Rif+PZA-->triple resistance
	}
	
	
  #Acquisition of resistance on FQ-based regimen
  acqr3<-dstlhs1[loop,3] #DS->Rif dstlhs1[loop,30] #
  acqf3<-dstlhs1[loop,35] #DS->FQ dstlhs1[loop,36] #
  acqp3<-dstlhs1[loop,5] #DS->PZA dstlhs1[loop,37] #
  acqr_rf3<-dstlhs1[loop,38] #Rif->Rif+FQ dstlhs1[loop,36] #
  acqr_rp3<-dstlhs1[loop,23] #Rif->Rif+PZA dstlhs1[loop,37] #
  acqf_rf3<-dstlhs1[loop,39] #FQ->Rif+FQ dstlhs1[loop,30] #
  acqf_fp3<-dstlhs1[loop,40] #FQ->FQ+PZA dstlhs1[loop,41] #
  acqp_rp3<-dstlhs1[loop,26] #PZA->Rif+PZA dstlhs1[loop,30] #
  acqp_fp3<-dstlhs1[loop,42] #PZA->FQ+PZA dstlhs1[loop,43] #
  acqrf_rfp3<-dstlhs1[loop,44] #Rif+FQ->triple dstlhs1[loop,41] #
  acqrp_rfp3<-dstlhs1[loop,45] #Rif+PZA->triple dstlhs1[loop,42] #
  acqfp_rfp3<-dstlhs1[loop,46] #FQ+PZA->triple dstlhs1[loop,30] #
  
  #Acquisition of resistance on DST-based regimen
  acqrf_rfp4<-dstlhs1[loop,47] #Rif+FQ->triple
  acqrp_rfp4<-dstlhs1[loop,48] #Rif+PZA->triple
  acqfp_rfp4<-dstlhs1[loop,49] #FQ+PZA->triple
  
  #General TB Nat Hx parameters
  tr <- dstlhs1[loop,50]/100000 # transmissions per person-yr
  lp <- 0.5 # reduction in reinfection rate if latently infected (as % susceptible to reinfection)
  fp <- 0.15 # proportion of infections progressing rapidly to active TB
  m1 <- 1/70 # baseline mortality rate, per yr (70-year LE)
  m2 <- 1/6 #active TB mortality rate
  sc <- 1/6 #rate of spontaneous cure without treatment (= rate of mortality, based on Vynnycky et al.)
  er <- 0.0007 # endogenous reactivation/relapse, per year (~5% lifetime)
  dx <- 1/(271/187) # rate of successful diagnosis of TB (but not drug resistance!) among rx-naive, per yr 
  # = Prev/Inc 
  relinf <- 0.2 #relative infectiousness when failing rx compared to untreated active TB
  
  #Rate of diagnosis for retreatment (per yr)
  dx_fail <- 2 #Assuming dx of failure after completion of 6-month regimen on average
  dx_rep <- dx #Dx rate for relapsed, reactivated, reinfected
  #Set similar to rate of initial dx to reflect delay to return to healthcare system
  #Incorporates time to seeking care, sensitivity of dx, repeated care-seeking
  
  #Rates of relapse by resistance state
  #Inputs based on Menzies PLoS Med 2009 papers
  rl <- 0.04 #DS-TB (Gillespie et al; WHO TB Report; Leung et al)
  rlr <- rl*4 #All strains resistant to Rif and Rif+other (Menzies et al 1-2 mo Rif vs. 6 mo Rif)
  rlf <- rl*3 #Mox and Mox+PZA resistant strains (Menzies et al: ~3-fold increase in relapse w/ 3 vs. 4 active drugs in IP, DS vs. INHr)
  rlp <- rl*2 #PZA-monoresistant strain (Menzies et al: 0.5 RR of relapse w/ PZA vs. no PZA)
  
  rlmat <- array(dim=c(8,4)) #Relapse rates for 8 strains x 4 regimens
  rlmat[,1] <- c(rl, rlr, rl, rlp, rlr, rlr, rlp, rlr) #HRZE (res to Mox does not matter)
  rlmat[,2] <- c(rl, rl, rlf, rlp, rlf, rlp, rlf, rlf) #STR
  rlmat[,3] <- c(rl, rlr, rlf, rlp, rlr, rlr, rlf, rlr) #REMox
  # For PaMZ: c(rl, rl, rlf, rlp, rlf, rlp, rlf, rlf)
  rlmat[,4] <- c(rl, rl, rl, rl, rlf, rlp, rl, rlf) #DST
  
  #Drug resistance Nat Hx and testing parameters
  #Matrix of probabilities of resistance acquisition: 8x8 DR states
  #Resistance states in order: 1=DS, 2=R, 3=F, 4=P, 5=RF, 6=RP, 7=FP, 8=RFP
  #rows=from, columns=to. E.g., [1,3] = prob for transition from state 1 to state 3
  #Simulating 3 eras: pre-HRZE (no acquisition of resistance to any drugs)
  #                   pre-FQ (acquisition of Rif and PZA resistance but not FQ)
  #                   current (acquisition of resistance to all drugs, availability of STR)
  
  #Create empty matrices 
  acq1<-matrix(rep(0,64),ncol=8) #For new cases receiving 1st-line rx (prior to intro of FQ regimen)
  colnames(acq1)<-rownames(acq1)<-c("DS","R","F","P","RF","RP","FP","RFP")
  
  acq3<-acq1 #STR, current era
  acq4<-acq1 #FQ-based 1st-line
  acq5<-acq1 #2nd-line pre-FQ regimen
  acq6<-acq1 #pre-FQ
  acq7<-acq1 #pre-HRZE
  
  #Fill in matrices based on DR state and choice of regimen
  #DS-TB, HRZE
  acq1[1,c(2:4)]<-c(acqr,acqf,acqp) #acquisition of monoresistance
  acq1[2,5:6]<-c(acqr_rf,acqr_rp) #Rif-> Rif+FQ and Rif+PZA
  acq1[3,c(5,7)]<-c(acqf_rf,acqf_fp) #FQ-> Rif+FQ and FQ+PZA
  acq1[4,6:7]<-c(acqp_rp,acqp_fp) #PZA->Rif+PZA and FQ+PZA
  acq1[5:7,8]<-c(acqrf_rfp,acqrp_rfp,acqfp_rfp) #Rif+FQ, Rif+PZA, and FQ+PZA->triple resistance
  
  #HRZE, Retreatment
  #Assuming prob of FQ resistance acquisition is higher 
  #Reflects higher chance of prior exposure to FQ
  #Modify values for acquisition of FQ resistance
  acq2<-acq1 #For retreatment cases receiving 1st-line rx (prior to intro of FQ regimen)
  acq2[1,3]<-acq2[1,3]*resmult
  acq2[2,5]<-acq2[2,5]*resmult
  acq2[4,7]<-acq2[4,7]*resmult
  acq2[6,8]<-acq2[6,8]*resmult
  
  #STR (No rif so 0 chance of Rif res acquisition)
  acq3[1,c(2:4)]<-c(0,dstlhs1[loop,35],dstlhs1[loop,5]) #acquisition of monoresistance
  acq3[2,5:6]<-c(acqr_rf2,acqr_rp2) #Rif-> Rif+FQ and Rif+PZA
  acq3[3,c(5,7)]<-c(0,dstlhs1[loop,40]) #FQ-> Rif+FQ and FQ+PZA; same as FQrx
  acq3[4,6:7]<-c(0,dstlhs1[loop,42]) #PZA->Rif+PZA and FQ+PZA; same as FQrx
  acq3[5:7,8]<-c(acqrf_rfp2,acqrp_rfp2,0) #Rif+FQ, Rif+PZA, and FQ+PZA->triple resistance
  
  #FQ-based 1st-line regimen
  acq4[1,c(2:4)]<-c(acqr3,acqf3,acqp3) #acquisition of monoresistance
  acq4[2,5:6]<-c(acqr_rf3,acqr_rp3) #Rif-> Rif+FQ and Rif+PZA
  acq4[3,c(5,7)]<-c(acqf_rf3,acqf_fp3) #FQ-> Rif+FQ and FQ+PZA
  acq4[4,6:7]<-c(acqp_rp3,acqp_fp3) #PZA->Rif+PZA and FQ+PZA
  acq4[5:7,8]<-c(acqrf_rfp3,acqrp_rfp3,acqfp_rfp3) #Rif+FQ, Rif+PZA, and FQ+PZA->triple resistance
  
  #DST-based regimen, pre-REMox
  acq5[1,c(2:4)]<-acq1[1,c(2:4)] #acquisition of monoresistance same as HRZE
  acq5[2,5:6]<-acq3[2,5:6] #Rif-> Rif+FQ and Rif+PZA same as STR
  acq5[3,c(5,7)]<-acq1[3,c(5,7)] #FQ-> Rif+FQ and FQ+PZA same as on HRZE
  acq5[4,6:7]<-c(dstlhs1[loop,3],dstlhs1[loop,35]) #PZA->Rif+PZA and FQ+PZA same as DS->Rif & DS->FQ on MRZE
  acq5[5:7,8]<-c(acqrf_rfp4,acqrp_rfp4,acqfp_rfp4) #Rif+FQ, Rif+PZA, and FQ+PZA->triple resistance
  
  acq8<-acq5 #2nd-line after intro of FQ regimen
  acq8[1,c(2:4)]<-acq4[1,c(2:4)] #acquisition of monoresistance same as MRZE
  
  
  #Pre-FQ; prior to availability of FQ, 0 prob of acquiring FQ res
  acq6[1,c(2:4)]<-c(acqr,0,acqp) #acquisition of monoresistance
  acq6[2,5:6]<-c(0,acqr_rp) #Rif-> Rif+FQ and Rif+PZA
  acq6[3,c(5,7)]<-c(0,0) #FQ-> Rif+FQ and FQ+PZA
  acq6[4,6:7]<-c(acqp_rp,0) #PZA->Rif+PZA and FQ+PZA
  acq6[5:7,8]<-c(0,0,0) #Rif+FQ, Rif+PZA, and FQ+PZA->triple resistance
  
  
  #sequential acquisition of resistance
  #Probabilities of acquiring resistance to multiple drugs in a single treatment episode
  #Assuming sequential acquisition of resistance
  #E.g., P(DS-->Rif+FQ)= P((DS-->Rif)*(Rif-->Rif+FQ)) + P((DS-->FQ)*(FQ-->Rif+FQ))
  #Acquisition of resistance to 1st drug increases prob of acquiring further resistance
  #Note: acq5 has 0 prob of FQ resistance, so 0 prob of DS-->Rif+FQ or DS-->FQ+PZA
  #Note: acq6 not included here as pre-treatment era, no acquisition of resistance
  
  #Double resistance acquisition from DS
  acq1[1,5]<-acq1[1,2]*acq1[2,5]+acq1[1,3]*acq1[3,5] #DS->Rif+FQ
  acq1[1,6]<-acq1[1,2]*acq1[2,6]+acq1[1,4]*acq1[4,6] #DS->Rif+PZA
  acq1[1,7]<-acq1[1,3]*acq1[3,7]+acq1[1,4]*acq1[4,7] #DS->FQ+PZA
  
  acq2[1,5]<-acq2[1,2]*acq2[2,5]+acq2[1,3]*acq2[3,5] #DS->Rif+FQ
  acq2[1,6]<-acq2[1,2]*acq2[2,6]+acq2[1,4]*acq2[4,6] #DS->Rif+PZA
  acq2[1,7]<-acq2[1,3]*acq2[3,7]+acq2[1,4]*acq2[4,7] #DS->FQ+PZA
  
  acq3[1,5]<-acq3[1,2]*acq3[2,5]+acq3[1,3]*acq3[3,5] #DS->Rif+FQ
  acq3[1,6]<-acq3[1,2]*acq3[2,6]+acq3[1,4]*acq3[4,6] #DS->Rif+PZA
  acq3[1,7]<-acq3[1,3]*acq3[3,7]+acq3[1,4]*acq3[4,7] #DS->FQ+PZA
  
  acq4[1,5]<-acq4[1,2]*acq4[2,5]+acq4[1,3]*acq4[3,5] #DS->Rif+FQ
  acq4[1,6]<-acq4[1,2]*acq4[2,6]+acq4[1,4]*acq4[4,6] #DS->Rif+PZA
  acq4[1,7]<-acq4[1,3]*acq4[3,7]+acq4[1,4]*acq4[4,7] #DS->FQ+PZA
  
  acq5[1,5]<-acq5[1,2]*acq5[2,5]+acq5[1,3]*acq5[3,5] #DS->Rif+FQ
  acq5[1,6]<-acq5[1,2]*acq5[2,6]+acq5[1,4]*acq5[4,6] #DS->Rif+PZA
  acq5[1,7]<-acq5[1,3]*acq5[3,7]+acq5[1,4]*acq5[4,7] #DS->FQ+PZA
  
  acq8[1,5]<-acq8[1,2]*acq8[2,5]+acq8[1,3]*acq8[3,5] #DS->Rif+FQ
  acq8[1,6]<-acq8[1,2]*acq8[2,6]+acq8[1,4]*acq8[4,6] #DS->Rif+PZA
  acq8[1,7]<-acq8[1,3]*acq8[3,7]+acq8[1,4]*acq8[4,7] #DS->FQ+PZA
  
  acq6[1,6]<-acq6[1,2]*acq6[2,6]+acq6[1,4]*acq6[4,6] #DS->Rif+PZA 
  
  #Single to triple resistance
  acq1[2,8]<-acq1[2,5]*acq1[5,8]+acq1[2,6]*acq1[6,8] #Rif->triple
  acq1[3,8]<-acq1[3,5]*acq1[5,8]+acq1[3,7]*acq1[7,8] #FQ->triple
  acq1[4,8]<-acq1[4,6]*acq1[6,8]+acq1[4,7]*acq1[7,8] #PZA->triple
  
  acq2[2,8]<-acq2[2,5]*acq2[5,8]+acq2[2,6]*acq2[6,8] #Rif->triple
  acq2[3,8]<-acq2[3,5]*acq2[5,8]+acq2[3,7]*acq2[7,8] #FQ->triple
  acq2[4,8]<-acq2[4,6]*acq2[6,8]+acq2[4,7]*acq2[7,8] #PZA->triple
  
  acq3[2,8]<-acq3[2,5]*acq3[5,8]+acq3[2,6]*acq3[6,8] #Rif->triple
  acq3[3,8]<-acq3[3,5]*acq3[5,8]+acq3[3,7]*acq3[7,8] #FQ->triple
  acq3[4,8]<-acq3[4,6]*acq3[6,8]+acq3[4,7]*acq3[7,8] #PZA->triple
  
  acq4[2,8]<-acq4[2,5]*acq4[5,8]+acq4[2,6]*acq4[6,8] #Rif->triple
  acq4[3,8]<-acq4[3,5]*acq4[5,8]+acq4[3,7]*acq4[7,8] #FQ->triple
  acq4[4,8]<-acq4[4,6]*acq4[6,8]+acq4[4,7]*acq4[7,8] #PZA->triple
  
  acq5[2,8]<-acq5[2,5]*acq5[5,8]+acq5[2,6]*acq5[6,8] #Rif->triple
  acq5[3,8]<-acq5[3,5]*acq5[5,8]+acq5[3,7]*acq5[7,8] #FQ->triple
  acq5[4,8]<-acq5[4,6]*acq5[6,8]+acq5[4,7]*acq5[7,8] #PZA->triple
  
  acq8[2,8]<-acq8[2,5]*acq8[5,8]+acq8[2,6]*acq8[6,8] #Rif->triple
  acq8[3,8]<-acq8[3,5]*acq8[5,8]+acq8[3,7]*acq8[7,8] #FQ->triple
  acq8[4,8]<-acq8[4,6]*acq8[6,8]+acq8[4,7]*acq8[7,8] #PZA->triple
  
  #Triple resistance acquisition from DS
  acq1[1,8]<-acq1[1,2]*acq1[2,8]+acq1[1,3]*acq1[3,8]+acq1[1,4]*acq1[4,8] #DS-> triple resistance
  acq2[1,8]<-acq2[1,2]*acq2[2,8]+acq2[1,3]*acq2[3,8]+acq2[1,4]*acq2[4,8] #DS-> triple resistance
  acq3[1,8]<-acq3[1,2]*acq3[2,8]+acq3[1,3]*acq3[3,8]+acq3[1,4]*acq3[4,8] #DS-> triple resistance
  acq4[1,8]<-acq4[1,2]*acq4[2,8]+acq4[1,3]*acq4[3,8]+acq4[1,4]*acq4[4,8] #DS-> triple resistance
  acq5[1,8]<-acq5[1,2]*acq5[2,8]+acq5[1,3]*acq5[3,8]+acq5[1,4]*acq5[4,8] #DS-> triple resistance
  acq8[1,8]<-acq8[1,2]*acq8[2,8]+acq8[1,3]*acq8[3,8]+acq8[1,4]*acq8[4,8] #DS-> triple resistance
  
  
  #Fill in probability of remaining in same resistance state
  # = 1 - (sum of probabilities of moving to any other resistance state)
  acq1[1,1]<-1-sum(acq1[1,-1])
  acq1[2,2]<-1-sum(acq1[2,-2])
  acq1[3,3]<-1-sum(acq1[3,-3])
  acq1[4,4]<-1-sum(acq1[4,-4])
  acq1[5,5]<-1-sum(acq1[5,-5])
  acq1[6,6]<-1-sum(acq1[6,-6])
  acq1[7,7]<-1-sum(acq1[7,-7])
  acq1[8,8]<-1-sum(acq1[8,-8])
  
  acq2[1,1]<-1-sum(acq2[1,-1])
  acq2[2,2]<-1-sum(acq2[2,-2])
  acq2[3,3]<-1-sum(acq2[3,-3])
  acq2[4,4]<-1-sum(acq2[4,-4])
  acq2[5,5]<-1-sum(acq2[5,-5])
  acq2[6,6]<-1-sum(acq2[6,-6])
  acq2[7,7]<-1-sum(acq2[7,-7])
  acq2[8,8]<-1-sum(acq2[8,-8])
  
  acq3[1,1]<-1-sum(acq3[1,-1])
  acq3[2,2]<-1-sum(acq3[2,-2])
  acq3[3,3]<-1-sum(acq3[3,-3])
  acq3[4,4]<-1-sum(acq3[4,-4])
  acq3[5,5]<-1-sum(acq3[5,-5])
  acq3[6,6]<-1-sum(acq3[6,-6])
  acq3[7,7]<-1-sum(acq3[7,-7])
  acq3[8,8]<-1-sum(acq3[8,-8])
  
  acq4[1,1]<-1-sum(acq4[1,-1])
  acq4[2,2]<-1-sum(acq4[2,-2])
  acq4[3,3]<-1-sum(acq4[3,-3])
  acq4[4,4]<-1-sum(acq4[4,-4])
  acq4[5,5]<-1-sum(acq4[5,-5])
  acq4[6,6]<-1-sum(acq4[6,-6])
  acq4[7,7]<-1-sum(acq4[7,-7])
  acq4[8,8]<-1-sum(acq4[8,-8])
  
  acq5[1,1]<-1-sum(acq5[1,-1])
  acq5[2,2]<-1-sum(acq5[2,-2])
  acq5[3,3]<-1-sum(acq5[3,-3])
  acq5[4,4]<-1-sum(acq5[4,-4])
  acq5[5,5]<-1-sum(acq5[5,-5])
  acq5[6,6]<-1-sum(acq5[6,-6])
  acq5[7,7]<-1-sum(acq5[7,-7])
  acq5[8,8]<-1-sum(acq5[8,-8])
  
  acq6[1,1]<-1-sum(acq6[1,-1])
  acq6[2,2]<-1-sum(acq6[2,-2])
  acq6[3,3]<-1-sum(acq6[3,-3])
  acq6[4,4]<-1-sum(acq6[4,-4])
  acq6[5,5]<-1-sum(acq6[5,-5])
  acq6[6,6]<-1-sum(acq6[6,-6])
  acq6[7,7]<-1-sum(acq6[7,-7])
  acq6[8,8]<-1-sum(acq6[8,-8])
  
  acq7[1,1]<-1-sum(acq7[1,-1])
  acq7[2,2]<-1-sum(acq7[2,-2])
  acq7[3,3]<-1-sum(acq7[3,-3])
  acq7[4,4]<-1-sum(acq7[4,-4])
  acq7[5,5]<-1-sum(acq7[5,-5])
  acq7[6,6]<-1-sum(acq7[6,-6])
  acq7[7,7]<-1-sum(acq7[7,-7])
  acq7[8,8]<-1-sum(acq7[8,-8])
  
  acq8[1,1]<-1-sum(acq8[1,-1])
  acq8[2,2]<-1-sum(acq8[2,-2])
  acq8[3,3]<-1-sum(acq8[3,-3])
  acq8[4,4]<-1-sum(acq8[4,-4])
  acq8[5,5]<-1-sum(acq8[5,-5])
  acq8[6,6]<-1-sum(acq8[6,-6])
  acq8[7,7]<-1-sum(acq8[7,-7])
  acq8[8,8]<-1-sum(acq8[8,-8])
  
  
  #Check that prob values add up to 1; if not, stop program
  for(i in 1:8){
    if(isTRUE(all.equal(rowSums(acq1)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 1")
    if(isTRUE(all.equal(rowSums(acq2)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 2")
    if(isTRUE(all.equal(rowSums(acq3)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 3")
    if(isTRUE(all.equal(rowSums(acq4)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 4")
    if(isTRUE(all.equal(rowSums(acq5)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 5")
    if(isTRUE(all.equal(rowSums(acq6)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 6")
    if(isTRUE(all.equal(rowSums(acq7)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 7")
    if(isTRUE(all.equal(rowSums(acq8)[i],1,check.attributes = F)) !=1) stop("error in resistance acquisition probs 8")
  }
  
  #Transmission fitness (as % transmissibility vs. DS-TB)
  fit_r<-dstlhs1[loop,6] #Relative fitness of Rif-resistant strain
  fit_f<-dstlhs1[loop,7] #Relative fitness of FQ-resistant srain compared to DS-TB
  fit_p<-dstlhs1[loop,8] #Relative fitness of PZA-resistant strain
  fit_rf<-dstlhs1[loop,9] #Relative fitness of Rif+FQ strain
  fit_rp<-dstlhs1[loop,10] #Relative fitness of Rif+PZA strain
  fit_fp<-dstlhs1[loop,11] #Relative fitness of FQ+PZA strain
  fit_r2<-dstlhs1[loop,12] #Relative fitness of Rif+2 strain
  
  #Cross-table of relative transmissibility for each pair of strains
  #Current strain in rows; superinfecting strain in columns
  #Matrix gives probability of superinfecting strain taking over (i.e., transition to new DR state)
  #Computed as (fitness of superinfecting strain)/(sum of fitness of current and superinfecting strains)
  
  cross<-matrix(nrow=8,ncol=8) #empty matrix
  colnames(cross)<-rownames(cross)<-c("ds","r","f","p","rf","rp","fp","r2") #assign row & column names for each DR state
  
  fit<-c(1,fit_r,fit_f,fit_p,fit_rf,fit_rp,fit_fp,fit_r2)
  for(i in 1:8){
    for (j in 1:8){
      cross[i,j]<-fit[j]/(fit[j]+fit[i])  #Prob of takeover by new= fit_new/(fit_new+fit_old)
    }
    cross[i,i]<-1 #Make the diagonal = 1 (otherwise with above formula, will all be = 0.5)
  }
  
    #INITIAL STATE PARAMETERS (steady-state conditions in pre-treatment era; no resistance)
    # Set initial matrices and population size
    # Reminder of the Matrix Structure:
    # matrix = 6 x 8
    # rows: 6 TB states (1=S, 2=L1, 3=A1, 4=F, 5=L2, 6=A2)
    # cols: 8 drug resistance states (1=DS, 2=Rif, 3=FQ, 4=PZA, 5=Rif+FQ, 6=Rif+PZA, 7=FQ+PZA, 8=Rif+2) 
    statevect<-stateinit[loop,]
    statename<-paste("s",1:48,sep="")  # 6x8=48 parameters
    names(statevect)<-statename #assign parameter names
  
  
  #Create master matrices for each era
  #Using master function, create transition matrix for each era, 
  #given resistance acquisition probabilities
  #Should use "master" function if REMox, "master2" function if PaMZ    
  
  #Prior to introduction of FQ-based regimen
  mastermat0<-(master(0,0,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp))
  
  mastermat1<-(master(1,0,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp))
  
  mastermat2<-(master(2,0,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp))
  
  #After introduction of FQ-based regimen
  mastermat2b<-(master(2,1,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp))
  
  
  # Put mastermat into a useable (vector) format
  mastervect0 <-c(mastermat0)
  mastervect1 <-c(mastermat1)
  mastervect2 <-c(mastermat2)
  mastervect2b <-c(mastermat2b)
  
  #Combine all into a single vector, with time of emergence of resistance
  mastervect<-c(mastervect0,mastervect1,mastervect2,mastervect2b,
                time1,time2)
  
  ##############################################################
  # SOLVE THE EQUATIONS!!
  
  # Times at which the model will be evaluated:
    times <- seq(61, 120, by = increment) #Complete duration of simulation; set first term to 0 if starting from beginning, 61 if starting in 2015
  
  out <- ode(statevect, times, TBdx, mastervect) #ODE solver function   
  
  Countermat<-array(dim=c(121,25)) #Create empty counter matrix
  for(i in 62:121){ #set to 1:121 if starting from beginning, 62:121 if starting in 2015
    Countermat[i,]<-Counter(out[(i-61),], #Set to i if starting from beginning, i-61 if from 2015
			    mastervect,increment) #Apply Counter function using output from ODE solver
  }
  
  #SELECTION OF TRAJECTORIES
  #Reconstruct population matrix at the end of simulation timeframe
  t <- out[,1] #time
  s <- out[,2] #Susceptible
  l1 <- out[,c(3,9,15,21,27,33,39,45)] #Latent rx-naive
  a1 <- out[,c(4,10,16,22,28,34,40,46)] #Active rx-naive
  f <- out[,c(5,11,17,23,29,35,41,47)] #Failing rx
  l2<- out[,c(6,12,18,24,30,36,42,48)] #Latent, rx-experienced
  a2<- out[,c(7,13,19,25,31,37,43,49)] #Active, rx-experienced
  
  total <-array(c(rowSums(out[,2:49]))) #Total population
  
  #Output values for selection of trajectories consistent with epi data
  atot<-rowSums(out[,c(4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49)]) #Total with active TB
  newtot<-rowSums(out[,c(4,10,16,22,28,34,40,46)]) #Total with active TB rx-naive
  rerxtot<-rowSums(out[,c(7,13,19,25,31,37,43,49)]) #Total with active TB, rx-experienced
  mdrnew<-rowSums(out[,c(10,28,34,46)])/newtot*100 #% with Rif-R active TB, rx-naive
  mdr<-rowSums(out[,c(13,31,37,49)])/rerxtot*100 #% with Rif-R active TB, rx-experienced
  pzanew<-out[,22]/newtot*100 #% with PZA monoresistance, rx-naive
  moxnew<-out[,16]/newtot*100 #% with FQ monoresistance, rx-naive
  rprerx<-rowSums(out[,c(37,49)])/rowSums(out[,c(13,31,37,49)])*100 #% with Rif+FQ resistance, rx-experienced
  rfrerx<-rowSums(out[,c(31,49)])/rowSums(out[,c(13,31,37,49)])*100 #% with Rif+PZA resistance, rx-experienced
  
  #Yearly incidence: overall, MDR, pre-XDR
  #For all below first dimension should be left blank if starting from beginning, 62:121 if starting from 2015
  incAll<-rowSums(Countermat[62:121,-25]) #Incidence overall 
  incMDR<-rowSums(Countermat[62:121,c(2,5,6,8,10,13,14,16,18,21,22,24)]) #Incidence of all MDR
  incXDR<-rowSums(Countermat[62:121,c(8,16,24)]) #Incidence of pre-XDR
  
 
  
  #Store population data into 3D array
  DSTlhsData[loop,62:121,]<-cbind(out, incAll, incMDR, incXDR)
  
#Overall incidence for years '90, '95, '00, '05, '10, '13
    incidence<-array(Countermat[60,-25],dim=c(3,8)) #Incidence for each strain
    inc13<-(sum(incidence)) #Total incidence
    
    incidence<-array(Countermat[57,-25],dim=c(3,8)) #Incidence for each strain
    inc10<-(sum(incidence)) #Total incidence

    incidence<-array(Countermat[52,-25],dim=c(3,8)) #Incidence for each strain
    inc05<-(sum(incidence)) #Total incidence
    
    incidence<-array(Countermat[47,-25],dim=c(3,8)) #Incidence for each strain
    inc00<-(sum(incidence)) #Total incidence
    
    incidence<-array(Countermat[42,-25],dim=c(3,8)) #Incidence for each strain
    inc95<-(sum(incidence)) #Total incidence
    
    incidence<-array(Countermat[37,-25],dim=c(3,8)) #Incidence for each strain
    inc90<-(sum(incidence)) #Total incidence
    
    #prevalence of MDR, PZA/FQ monoresistance, and MDR+PZA/FQ resistance 
    #for years '90, '95, '00, '05, '13
    mdrlo13<-mdrnew[60]
    mdrlo05<-mdrnew[52]
    mdrlo00<-mdrnew[47]
    mdrlo95<-mdrnew[42]
    mdrlo90<-mdrnew[37]
    
    mdrfinal13<-mdr[60]
    mdrfinal05<-mdr[52]
    mdrfinal00<-mdr[47]
    mdrfinal95<-mdr[42]
    mdrfinal90<-mdr[37]
    
    pzafinal13<-pzanew[60]
    pzafinal05<-pzanew[52]
    pzafinal00<-pzanew[47]
    pzafinal95<-pzanew[42]
    pzafinal90<-pzanew[37]
    
    moxfinal13<-moxnew[60]
    moxfinal05<-moxnew[52]
    moxfinal00<-moxnew[47]
    moxfinal95<-moxnew[42]
    moxfinal90<-moxnew[37]
    
    rpfinal13<-rprerx[60]
    rpfinal05<-rprerx[52]
    rpfinal00<-rprerx[47]
    rpfinal95<-rprerx[42]
    rpfinal90<-rprerx[37]
    
    rffinal13<-rfrerx[60]
    rffinal05<-rfrerx[52]
    rffinal00<-rfrerx[47]
    rffinal95<-rfrerx[42]
    rffinal90<-rfrerx[37]
    
    #Criteria for selection of simulations (may modify later)
    #incidence +/- 25%
    #MDR prevalence in retreatment +/- 25%
    #PZA and Mox monoresistance in new cases <5%
    #MDR in retreatment > in new
    #% resistant to mox or pza among mdr
    
    #Create empty vector to keep track of criteria met
    score<-rep(0,12)
    
    #Store input values, incidence and resistance prevalence, and selection criteria for each simulation
    dstlhs2[loop,]<-c(dstlhs1[loop,],inc13, inc10, inc05, inc00, inc95, inc90,
                      mdrfinal13, mdrlo13,pzafinal13,moxfinal13,rpfinal13,rffinal13,
                      score)

 
  #Keep track of: 
  print(loop) #Total simulations completed
  print(row) #Total simulations meeting selection criteria
  
} #END OF SIMULATION LOOP


# #WRITE THE DATA!!!
# #Uncomment this section and 
# #Modify file path/names as necessary
salist <- c("AllPZA", 
						"F-RF", "R-RF", "P-RP", "R-RP", 
						"FP-RFP", "RP-RFP", 
						"R-RP-RFP", "P-RP-RFP")
write.csv(DSTlhsData,file=paste("SA/Data/dstData110915_", salist[taskid], ".csv", sep=""))



