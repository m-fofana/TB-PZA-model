#Mariam Fofana
#Last updated: October 2015
#Model of DST for new fluoroquinolone-based first-line TB regimens

#This script is a compilation of the functions used in the model
#Run this prior to running the simulation loop

#############################################
#FUNCTION I: CREATE MASTER TRANSITION MATRIX#
#############################################

#This function determines the appropriate fixed transition rates
#(i.e., not dependent on force of infection)
#(e.g., active DS-TB --> active MDR TB)
#Inputs are: (1) era; era 0:pre-treatment; era 1: resistance to HRZE; era 2: resistance to HRZE and FQ
#            (2) opt; this is a 0/1 indicator for use of FQ-based first-line regimen
#            (3) acq1-acq7; these are 7 sets of values for probability of resistance acquisition,
#                 conditional on drug resistance state and choice of treatment regimen

master<-function(era,opt,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp){ #Define function inputs
    
  era<-era  #Assign values to variables
  opt<-opt
  acq1<-acq1
  acq2<-acq2
  acq3<-acq3
  acq4<-acq4
  acq5<-acq5
  acq6<-acq6
  acq7<-acq7
  acq8<-acq8
  
  comp<-comp
    
  #Create 4-dimensional array of transition rates for master matrix: 6x8x6x8
  mastermat <- array(0, dim=c(6,8,6,8)) 
  
  #Dimensions 1 & 2 indicate origin compartment; 6 nat hx states x 8 drug resistance states
  #Dimensions 3 & 4 indicate destination compartment; 6 nat hx states x 8 drug resistance states
  i <- rep(1:6, 8) # Rows; 1,2,3,4,5,6,1,2,3,4,5,6,etc...
  j <- rep(1:8, each=6) # Columns; 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,etc...
  ##i <- rep(1:6,each=8) #Rows; This creates 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4, etc...
  ##j <- rep(1:8,6) #Columns; This creates 1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6, etc...
  k <- array(c(i,j),dim=c(48,2)) #This creates combo of row/column for every possible state
  m <- array(c(k,k),dim=c(48,4)) #This creates combo for every  transition of state x-->state x (for mortality)
  
  ####################################################################################################################
  ##### MORTALITY #####
  
  mortmat <- array(0, dim=c(6,8))  # create empty matrix
  mortmat[,]<-m1  #Fill in with background mortality
  mortmat[3,]<-mortmat[3,]+m2  #Add active TB mortality for naive active TB compartment
  mortmat[6,]<-mortmat[6,]+m2  #Add active TB mortality for retreatment active TB compartment
  
  #No added mortality for failure compartment as assumed to have baseline mortality while receiving partially effective rx
  #Can change this by modifying and uncommenting line below
  #mortmat[4,]<-mortmat[4,]+0
  
  mortvect <-c(mortmat) #Convert mortality matrix to vector
  mastermat[m]<-mastermat[m]+mortvect #Add mortality rates into master matrix
  
  ####################################################################################################################
  ###### RESISTANCE ACQUISITION  ######
  
  #Switch between eras
  #If modeling prior to intro of HRZE (era 0), model will run with acq7 values
  #If modeling prior to intro of FQ (era 1), model will run with acq6 values
  #If modeling post intro of resistance to all drugs incl. FQ (era 2), model will use acq1-acq5
  acq<-array(,dim=c(8,8,5))
  dimnames(acq)[[1]]<-dimnames(acq)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
  dimnames(acq)[[3]]<-c("HRZE new","HRZE rerx","STR","MRZE","DST")
  
  if(era==0 & opt==0){
    acq[,,1]<-acq7 #No rx, so same probs for all
    acq[,,2]<-acq7 
    acq[,,3]<-acq7 
    acq[,,4]<-acq7
    acq[,,5]<-acq7
  }else if(era==1 & opt==0){
    acq[,,1]<-acq6 #No resistance to FQ, so same set of probs for all
    acq[,,2]<-acq6
    acq[,,3]<-acq6
    acq[,,4]<-acq6
    acq[,,5]<-acq6
  }else if(era==2 & opt==0){
    acq[,,1]<-acq1 #HRZE, new
    acq[,,2]<-acq2 #HRZE, rerx
    acq[,,3]<-acq3 #STR
    acq[,,4]<-acq4 #FQrx
    acq[,,5]<-acq5 #DST
  }else if(era==2 & opt==1){
    acq[,,1]<-acq1 #HRZE, new
    acq[,,2]<-acq2 #HRZE, rerx
    acq[,,3]<-acq3 #STR
    acq[,,4]<-acq4 #FQrx
    acq[,,5]<-acq8 #DST
  }

   
  #Array of true vs. diagnosed drug resistance (conditional on availability of DST)
  #True DR in rows, dx'd DR in columns, 3rd dimension is new vs retreatment vs. failure
  dxmat<-array(0,dim=c(8,8,3))
  dimnames(dxmat)[[1]]<-dimnames(dxmat)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
  dimnames(dxmat)[[3]]<-c("New","Rep","Fail")
                    
  #DS is always dx'd as DS (perfect specificity)
  dxmat[1,1,]<-1
  
  for(i in 1:3){ #Loop to fill in values for New, Rep, Fail
      #Rif-R always dx'd as either DS or Rif-R
      dxmat[2,2,i]<-DSTr_sens*(dst_access[i,2,(opt+1)]+ #Rif only DST
                                 dst_access[i,5,(opt+1)]+ #Rif+FQ DST
                                 dst_access[i,6,(opt+1)]+ #Rif+PZA DST
                                 dst_access[i,8,(opt+1)]) #comprehensive DST
      dxmat[2,1,i]<-1-dxmat[2,2,i]
      
      #FQ-R
      dxmat[3,3,i]<-DSTf_sens*(dst_access[i,3,(opt+1)]
                               +dst_access[i,5,(opt+1)]
                               +dst_access[i,7,(opt+1)]
                               +dst_access[i,8,(opt+1)]) #Correctly dx'd
      dxmat[3,1,i]<-1-dxmat[3,3,i]
        
      #PZA-R
      dxmat[4,4,i]<-DSTp_sens*(dst_access[i,4,(opt+1)]
                               +dst_access[i,6,(opt+1)]
                               +dst_access[i,7,(opt+1)]
                               +dst_access[i,8,(opt+1)]) #Correctly dx'd
      dxmat[4,1,i]<-1-dxmat[4,4,i]
      
      #Rif+FQ
      dxmat[5,5,i]<-DSTrf_sens*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
      dxmat[5,2,i]<-(DSTr_sens*(1-DSTf_sens))*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as Rif-R only after testing for both Rif and FQ
                  DSTr_sens*(dst_access[i,2,(opt+1)]+dst_access[i,6,(opt+1)]) #dx'd as Rif-R only after testing for Rif only or RP
      dxmat[5,3,i]<-(DSTf_sens*(1-DSTr_sens))*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as FQ-R only after testing for both Rif and FQ
                  DSTf_sens*(dst_access[i,3,(opt+1)]+dst_access[i,7,(opt+1)]) #dx'd as FQ-R only after testing for FQ only or FP
      dxmat[5,1,i]<-(1-DSTr_sens)*(1-DSTf_sens)*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both Rif and FQ
                  (1-DSTr_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS after testing for RP
                  (1-DSTf_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FP
                  (1-DSTr_sens)*dst_access[i,2,(opt+1)]+ #dx'd as DS only after testing for Rif only
                  (1-DSTf_sens)*dst_access[i,3,(opt+1)]+ #dx'd as DS only after testing for FQ only
                  dst_access[i,4,(opt+1)]+ #Received DST for P only, assumed to be DS
                  dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
         
      #Rif+PZA
      dxmat[6,6,i]<-DSTrp_sens*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
      dxmat[6,2,i]<-(DSTr_sens*(1-DSTp_sens))*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as Rif-R only after testing for both Rif and PZA
                  DSTr_sens*(dst_access[i,2,(opt+1)]+dst_access[i,5,(opt+1)]) #dx'd as Rif-R only after testing for Rif only or RF
      dxmat[6,4,i]<-(DSTp_sens*(1-DSTr_sens))*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as PZA-R only after testing for both Rif and PZA
                  DSTp_sens*(dst_access[i,4,(opt+1)]+dst_access[i,7,(opt+1)]) #dx'd as PZA-R only after testing for PZA only or FP
      dxmat[6,1,i]<-(1-DSTr_sens)*(1-DSTp_sens)*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both Rif and PZA
                  (1-DSTr_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS after testing for RF
                  (1-DSTp_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FP
                  (1-DSTr_sens)*dst_access[i,2,(opt+1)]+ #dx'd as DS only after testing for Rif only
                  (1-DSTp_sens)*dst_access[i,4,(opt+1)]+ #dx'd as DS only after testing for PZA only
                  dst_access[i,3,(opt+1)]+ #Received DST for F only, assumed to be DS
                  dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
        
      #FQ+PZA
      dxmat[7,7,i]<-DSTfp_sens*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
      dxmat[7,3,i]<-(DSTf_sens*(1-DSTp_sens))*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as FQ-R only after testing for both FQ and PZA
                  DSTf_sens*(dst_access[i,3,(opt+1)]+dst_access[i,5,(opt+1)]) #dx'd as FQ-R only after testing for FQ only or RF
      dxmat[7,4,i]<-(DSTp_sens*(1-DSTf_sens))*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as PZA-R only after testing for both FQ and PZA
                  DSTp_sens*(dst_access[i,4,(opt+1)]+dst_access[i,6,(opt+1)]) #dx'd as PZA-R only after testing for PZA only or RP
      dxmat[7,1,i]<-(1-DSTf_sens)*(1-DSTp_sens)*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both FQ and PZA
                  (1-DSTf_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS after testing for RF
                  (1-DSTp_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS after testing for RP
                  (1-DSTf_sens)*dst_access[i,3,(opt+1)]+ #dx'd as DS only after testing for FQ only
                  (1-DSTp_sens)*dst_access[i,4,(opt+1)]+ #dx'd as DS only after testing for PZA only
                  dst_access[i,2,(opt+1)]+ #Received DST for R only, assumed to be DS
                  dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
        
      #Triple-resistant
      dxmat[8,8,i]<-DSTx_sens*dst_access[i,8,(opt+1)] #Correcty dx'd
      dxmat[8,7,i]<-DSTfp_sens*(1-DSTr_sens)*dst_access[i,8,(opt+1)]+ #dx'd as FQ+PZA after testing for all
                  DSTfp_sens*dst_access[i,7,(opt+1)] #Dx'd as FQ+PZA after testing for FQ and PZA
      dxmat[8,6,i]<-DSTrp_sens*(1-DSTf_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif+PZA after testing for all
                  DSTrp_sens*dst_access[i,6,(opt+1)] #Dx'd as Rif+PZA after testing for Rif and PZA  
      dxmat[8,5,i]<-DSTrf_sens*(1-DSTp_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif+FQ after testing for all
                  DSTrf_sens*dst_access[i,5,(opt+1)] #Dx'd as Rif+FQ after testing for FQ and Rif
      dxmat[8,4,i]<-(DSTp_sens*(1-DSTr_sens)*(1-DSTf_sens))*dst_access[i,8,(opt+1)]+ #dx'd as PZA-R only after testing for all
                  (DSTp_sens*(1-DSTf_sens))*dst_access[i,7,(opt+1)]+ #dx'd as PZA-R only after testing for PZA and FQ
                  (DSTp_sens*(1-DSTr_sens))*dst_access[i,6,(opt+1)]+ #dx'd as PZA-R only after testing for PZA and Rif
                  DSTp_sens*dst_access[i,4,(opt+1)]#dx'd as PZA-R only after testing for PZA only
      dxmat[8,3,i]<-(DSTf_sens*(1-DSTr_sens)*(1-DSTp_sens))*dst_access[i,8,(opt+1)]+ #dx'd as FQ-R only after testing for all
                  (DSTf_sens*(1-DSTp_sens))*dst_access[i,7,(opt+1)]+ #dx'd as FQ-R only after testing for PZA and FQ
                  (DSTf_sens*(1-DSTr_sens))*dst_access[i,5,(opt+1)]+ #dx'd as FQ-R only after testing for FQ and Rif
                  DSTf_sens*dst_access[i,3,(opt+1)]#dx'd as FQ-R only after testing for FQ only
      dxmat[8,2,i]<-(DSTr_sens*(1-DSTf_sens)*(1-DSTp_sens))*dst_access[i,8,(opt+1)]+ #dx'd as Rif-R only after testing for all
                  (DSTr_sens*(1-DSTp_sens))*dst_access[i,6,(opt+1)]+ #dx'd as Rif-R only after testing for PZA and Rif
                  (DSTr_sens*(1-DSTf_sens))*dst_access[i,5,(opt+1)]+ #dx'd as Rif-R only after testing for FQ and Rif
                  DSTr_sens*dst_access[i,2,(opt+1)]#dx'd as Rif-R only after testing for Rif only
      dxmat[8,1,i]<-(1-DSTr_sens)*(1-DSTf_sens)*(1-DSTp_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif-R only after testing for all
                  (1-DSTf_sens)*(1-DSTp_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FQ and PZA
                  (1-DSTr_sens)*(1-DSTp_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS only after testing for PZA and Rif
                  (1-DSTr_sens)*(1-DSTf_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS only after testing for FQ and Rif
                  (1-DSTp_sens)*dst_access[i,4,(opt+1)]+#dx'd as DS after testing for PZA only
                  (1-DSTf_sens)*dst_access[i,3,(opt+1)]+#dx'd as DS  after testing for FQ only
                  (1-DSTr_sens)*dst_access[i,2,(opt+1)]+#dx'd as DS after testing for Rif only
                  dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
  }
      
 
  #Check rowSums to make sure that they add up to 1
  for(i in 1:8){
    if(isTRUE(all.equal(rowSums(dxmat)[i],3,check.attributes = F)) !=1) stop("error in dxmat")
  }
  
  ####################################################################################################################
  ##### PROBABILITIES OF CURE #####
  
  #Probabilities of cure by DR state and choice of treatment regimen
  rx<-array(0,dim=c(8,6)) #Create empty array
  rownames(rx)<-c("DS","R","F","P","RF","RP","FP","RFP") #Rows: DR states
  colnames(rx)<-c("HRZE","MRZE","DST","ReRx","STR","Def") #Columns: Regimen
  
    #Each row = DR state, each column= (1) first-line HRZE; (2) FQ-based first-line (MRZE), 
    #(3) individualized DST-based treatment, (4) Empiric regimen for retreatment cases ("category 2"),
    #(5) standardized MDR regimen, (6) special values for default cases
  
    #Assume that empiric retreatment regimen is no more efficacious than first-line rx
  
  #First-line
  ##DS-TB
  #Prob of cure for DS-TB assumed to be the same regardless of treatment type
  rx[1,1]<-dstlhs1[loop,13]
  ##Monoresistance to Rif, Rif/FQ
  rx[c(2,5),1]<-dstlhs1[loop,14]
  ##Monoresistance to FQ
  rx[3,1]<-rx[1,1] #Same as DS
  ##Monoresistance to PZA, FQ/PZA
  rx[c(4,7),1]<-dstlhs1[loop,15]
  ##Rif/PZA and triple resistance
  rx[c(6,8),1]<-dstlhs1[loop,16]

  
  #2nd-line treatment
  rx[c(1,3,4,7),5]<-rx[c(1,3,4,7),1] #Same as HRZE
  ##Rif monoresistance
  rx[2,5]<-dstlhs1[loop,17]
  ##Rif/FQ
  rx[5,5]<-dstlhs1[loop,18]
  ##Rif/PZA
  rx[6,5]<-dstlhs1[loop,19]
  ##Triple resitance
  rx[8,5]<-dstlhs1[loop,20]

  #MRZE [keep same as HRZE for now]
  rx[,2]<-rx[,1]
  
  #DST [keep same as STR for now]
  rx[,3]<-rx[,1]
  
  ##Retreatment regimens
  #Set empiric retreatment regimen as similar to HRZE (if opt=0) or MRZE (if opt=1)
  rx[,4]<-rx[,opt+1] 
  
  #Defaulters
  rx[,6]<-c(0.5,0,0.25,0.25,0,0,0,0)
  
  ####################################################################################################################
  ##### TREATMENT COMPLETION #####

  compmat<-array(dim=c(8,5)) #8 DR states (as diagnosed) x 5 regimens
  rownames(compmat)<-c("DS","R","F","P","RF","RP","FP","RFP") #Rows: DR states
  colnames(compmat)<-c("HRZE","MRZE","DST","ReRx","STR") #Columns: Regimen
  compmat[,1:2]<-rep(comp[1],8) #first-line
  compmat[,4]<-rep(comp[2],8) #ReRx
  compmat[,5]<-rep(comp[3],8) #STR
  
  #Values for DST-based regimen               
  compmat[,3]<-c(comp[1], #DS would get HRZE/MRZE
                 comp[3], #Rif-R would get a 2nd-line rx
                 comp[1], #FQ-R would get HRZE
                 comp[2], #PZA-R would get HRE/9 mo or similar
                 comp[3], #Rif+FQ-R would get 2nd-line rx
                 comp[3], #Rif+PZA-R would get 2nd-line rx
                 comp[2], #FQ+PZA-R would get HRE/9 mo or similar
                 comp[3]) #Triple-res would get 2nd-line rx
                 
  
  #Matrices of treatment completion and default based on dx'd DR, new or ReRx 
  #8 DR states (rows) by 4 % completion (columns; HRZE, MRZE, STR, DST)
  #Based on diagnosed DR, how likely are you to complete treatment?
  dxrx_comp_new <- dxrx[,,(1+opt)]*compmat[,c(1,2,5,3)] #Complete rx, new cases
  dxrx_comp_rerx <- dxrx[,,(1+opt)]*compmat[,c(4,4,5,3)] #Complete rx, retreatment
  dxrx_def_new <- dxrx[,,(1+opt)]*(1-compmat[,c(1,2,5,3)]) #Default, new cases
  dxrx_def_rerx <- dxrx[,,(1+opt)]*(1-compmat[,c(4,4,5,3)]) #Default, retreatment
    
     
  ####################################################################################################################
  ##### TREATMENT OUTCOMES #####
  
  #Treatment outcomes are conditional on DR state and choice of regimen
  #5 possible outcomes: dc (cure after default), df (continued active TB after default), rx (cure),
  #                     rl (relapse after initial cure), fl (treatment failure, continued active TB)
  
  #Variables below are named based on: (1) outcome; (2) DR state and new vs. retreatment; (3) choice of regimen
  
  #Note: In order to allow for uncertainty around prob of cure, each simulation includes a randomly sampled "failure multiplier"
  #Final prob of cure = 1 - (1 -  baseline prob)*failure multiplier
    
  #Matrix of proportion getting each regimen
  #Dot product of dxmat and dxrx
  #Based on initial DR, what is your prob of cure?
  #Will make matrix with (1) true initial DR state in rows, (2) rx in column
  #Add 3rd dimension for new, rep, fail
  
  #New cases
  rxmat<-array(,dim=c(8,4,3,2))
  dimnames(rxmat)[[1]]<-c("DS","R","F","P","RF","RP","FP","RFP")
  dimnames(rxmat)[[2]]<-c("HRZE","MRZE","STR","DST")
  dimnames(rxmat)[[3]]<-c("New","Rep","Fail")
  dimnames(rxmat)[[4]]<-c("Comp","Def") #Completers vs defaulters
  
  rxmat[,,1,1]<-dxmat[,,1] %*% dxrx_comp_new 
  rxmat[,,2,1]<-dxmat[,,2] %*% dxrx_comp_rerx
  rxmat[,,3,1]<-dxmat[,,3] %*% dxrx_comp_rerx
  rxmat[,,1,2]<-dxmat[,,1] %*% dxrx_def_new 
  rxmat[,,2,2]<-dxmat[,,2] %*% dxrx_def_rerx
  rxmat[,,3,2]<-dxmat[,,3] %*% dxrx_def_rerx
  
  #Check sums--should add up to 1
  for(i in 1:3){
    for(j in 1:8){
      if(isTRUE(all.equal(rowSums(rxmat[,,i,])[j],1,check.attributes = F)) !=1) stop("error in rxmat")
    }
  }
  
  
  #Adjusted prob of cure (taking into account failure multiplier)
  rx2<-rx

  #Treatment completion
  
  
  #   #Final prob of cure is dependent on choice of rx and final DR state
  #   #For new cases, first-line regimen
  #final outcomes array
  #8 rows (initial DR) x 8 columns (final DR) x 5 layers (each outcome)
  finalout<-array(0, dim=c(8,8,5,3))
  dimnames(finalout)[[1]]<-dimnames(finalout)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
  dimnames(finalout)[[3]]<-c("Cure","Relapse","Failure","Def/cure","Def/fail")
  dimnames(finalout)[[4]]<-c("New","Rep","Fail")  
  
  for(h in 1:3){ 
    tmp<-ceiling(0.5+h/3) #Will round to 1 for New, 2 for Rep and Fail
  
    for(i in 1:8){ #Each row of the rxmat
    
      #Cured
      finalout[i,,1,h]<-
        rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE 
          rx2[,1]*(1-rlmat[,1])+
          rxmat[i,3,h,1]*acq[i,,3]* #STR
          rx2[,5]*(1-rlmat[,2])+
          rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
          rx2[,2]*(1-rlmat[,1])+
          rxmat[i,4,h,1]*acq[i,,5]* #DST
          rx2[,3]*(1-rlmat[,2])
      
      #Relapse
      finalout[i,,2,h]<-
        rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE 
         rx2[,1]*rlmat[,1]+
         rxmat[i,3,h,1]*acq[i,,3]* #STR
         rx2[,5]*rlmat[,2]+
         rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
         rx2[,2]*rlmat[,1]+
         rxmat[i,4,h,1]*acq[i,,5]* #DST
         rx2[,3]*rlmat[,2]
        
      #Fail
      finalout[i,,3,h]<-
        rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE new
        (1-rx2[,1])+
        rxmat[i,3,h,1]*acq[i,,3]* #STR
        (1-rx2[,5])+
        rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
        (1-rx2[,2])+
        rxmat[i,4,h,1]*acq[i,,5]* #DST
        (1-rx2[,3])
      
      #Cure after default
      finalout[i,,4,h]<-
        (rxmat[i,1,h,2]*acq[i,,tmp]+ #HRZE new
        rxmat[i,3,h,2]*acq[i,,3]+ #STR
        rxmat[i,2,h,2]*acq[i,,4]+ #FQ 1st-line
        rxmat[i,4,h,2]*acq[i,,5])* #DST
        rx2[,6]
      
      #Continued active TB after default
      finalout[i,,5,h]<-
        (rxmat[i,1,h,2]*acq[i,,tmp]+ #HRZE new
         rxmat[i,3,h,2]*acq[i,,3]+ #STR
         rxmat[i,2,h,2]*acq[i,,4]+ #FQ 1st-line
         rxmat[i,4,h,2]*acq[i,,5])* #DST
         (1-rx2[,6])
      } #End of i loop
    } #End of h loop

  #Check that probabilities add up to 1
#   for(h in 1:3){
#     for(i in 1:8){
#   print(sum(finalout[i,,,h]))
#     }
#   }

#FILL IN MASTERMAT NOW!!!

  ####################################################################################################################
  ##### DISEASE PROGRESSION #####
  
  #Mastermat array dimensions 1 & 2 indicate ORIGIN compartment; 6 nat hx states x 8 drug resistance states
  #Dimensions 3 & 4 indicate DESTINATION compartment; 6 nat hx states x 8 drug resistance states
  #Nat hx states: (1) Susceptible (2) Latent rx-naive (3) Active TB x-naive 
  #               (4) Failing (5) Latent rx-experienced (6) Active TB rx-experienced
  #DR states: (1) DS (2) Rif only (3) FQ only (4) PZA only
  #           (5) Rif & FQ (6) Rif & PZA (7) PZA &FQ (8) RIF & FQ & PZA
  
  
  ##L->A progression from latent to active TB: endogenous reactivation
  ##No acquisition of resistance in this transition
  #New cases
  mastermat[2,1,3,1]<-er #From latent rx-naive DS-TB to active rx-naive DS-TB
  mastermat[2,2,3,2]<-er #From latent rx-naive Rif-R TB to active rx-naive Rif-R TB
  mastermat[2,3,3,3]<-er
  mastermat[2,4,3,4]<-er
  mastermat[2,5,3,5]<-er
  mastermat[2,6,3,6]<-er
  mastermat[2,7,3,7]<-er
  mastermat[2,8,3,8]<-er
  
  #Retreatment cases
  mastermat[5,1,6,1]<-er 
  mastermat[5,2,6,2]<-er
  mastermat[5,3,6,3]<-er
  mastermat[5,4,6,4]<-er
  mastermat[5,5,6,5]<-er
  mastermat[5,6,6,6]<-er
  mastermat[5,7,6,7]<-er
  mastermat[5,8,6,8]<-er 
  
  ##A->L spontaneous cure w/o treatment
  ##No acquisition of resistant in this transition
  #New cases
  mastermat[3,1,2,1]<-sc #From active rx-naive DS-TB to latent rx-naive DS-TB
  mastermat[3,2,2,2]<-sc #From active rx-naive Rif-R TB to latent rx-naive Rif-R TB
  mastermat[3,3,2,3]<-sc
  mastermat[3,4,2,4]<-sc
  mastermat[3,5,2,5]<-sc
  mastermat[3,6,2,6]<-sc
  mastermat[3,7,2,7]<-sc
  mastermat[3,8,2,8]<-sc
  
  #Retreatment cases
  mastermat[6,1,5,1]<-sc 
  mastermat[6,2,5,2]<-sc
  mastermat[6,3,5,3]<-sc
  mastermat[6,4,5,4]<-sc
  mastermat[6,5,5,5]<-sc
  mastermat[6,6,5,6]<-sc
  mastermat[6,7,5,7]<-sc
  mastermat[6,8,5,8]<-sc

  ###Initial treatment
  ##A1->L2: successful 1st treatment; resistance can be acquired during treatment, such that
  #         transition rate = (dx rate) * (prob of res acquisition) * (prop rx success + prop cure after default)
  
  #DS
  mastermat[3,1,5,]<-dx * (finalout[1,,1,1]+finalout[1,,4,1])
  
  #Rif monores
  # Rate is computed in the same manner as above but requires accounting for missed diagnosis of drug resistance
  # Prob of res acq and prob of cure are weighted averages for 
  #(1) those correclty dx'd as Rif-R and (2) those misdx'd as DS-TB
  mastermat[3,2,5,]<-dx * (finalout[2,,1,1]+finalout[2,,4,1])
    
  #FQ monores
  mastermat[3,3,5,]<-dx *(finalout[3,,1,1]+finalout[3,,4,1]) 
    
  #PZA
  mastermat[3,4,5,]<-dx * (finalout[4,,1,1]+finalout[4,,4,1])
    
  #Rif+FQ
  mastermat[3,5,5,]<-dx * (finalout[5,,1,1]+finalout[5,,4,1])
  
  #Rif+PZA
  mastermat[3,6,5,]<-dx * (finalout[6,,1,1]+finalout[6,,4,1])
  
  #FQ+PZA
  mastermat[3,7,5,]<-dx * (finalout[7,,1,1]+finalout[7,,4,1])
  
  #Rif+2
  mastermat[3,8,5,]<-dx * (finalout[8,,1,1]+finalout[8,,4,1])

    
  ##A1->F failure after initial treatment: rate of dx * prop fail
  #DS
  mastermat[3,1,4,]<-dx*finalout[1,,3,1]

  #Rif
  mastermat[3,2,4,]<-dx*finalout[2,,3,1]
    
  #FQ
  mastermat[3,3,4,]<-dx*finalout[3,,3,1]
  
  #PZA
  mastermat[3,4,4,]<-dx*finalout[4,,3,1]
  
  #Rif+FQ
  mastermat[3,5,4,]<-dx*finalout[5,,3,1]
    
  #Rif+PZA
  mastermat[3,6,4,]<-dx*finalout[6,,3,1]
    
  #PZA+FQ
  mastermat[3,7,4,]<-dx*finalout[7,,3,1]
  
  #Rif+2
  mastermat[3,8,4,]<-dx*finalout[8,,3,1]
  
  
  ##A1->A2 recurrence after default and Relapse
  #DS
  mastermat[3,1,6,]<-dx*(finalout[1,,5,1]+finalout[1,,2,1])
    
  #Rif
  mastermat[3,2,6,]<-dx*(finalout[2,,5,1]+finalout[2,,2,1])
  
  #FQ
  mastermat[3,3,6,]<-dx*(finalout[3,,5,1]+finalout[3,,2,1])
  
  #PZA
  mastermat[3,4,6,]<-dx*(finalout[4,,5,1]+finalout[4,,2,1])
      
  #Rif+FQ
  mastermat[3,5,6,]<-dx*(finalout[5,,5,1]+finalout[5,,2,1])
  
  #Rif+PZA
  mastermat[3,6,6,]<-dx*(finalout[6,,5,1]+finalout[6,,2,1])
  
  #PZA+FQ
  mastermat[3,7,6,]<-dx*(finalout[7,,5,1]+finalout[7,,2,1])
    
  #Rif+2
  mastermat[3,8,6,]<-dx*(finalout[8,,5,1]+finalout[8,,2,1])
  
  
  ##RETREATMENT 
  #A2->L2; F->L2 successful retreatment: rate of dx * (prop treatment success+prop cure after default)
  #DS
  mastermat[6,1,5,]<-dx_rep*(finalout[1,,1,2]+finalout[1,,4,2])
  mastermat[4,1,5,]<-dx_fail*(finalout[1,,1,3]+finalout[1,,4,3])
  
  #Rif
  mastermat[6,2,5,]<-dx_rep*(finalout[2,,1,2]+finalout[2,,4,2])
  mastermat[4,2,5,]<-dx_fail*(finalout[2,,1,3]+finalout[2,,4,3])
  
  #FQ
  mastermat[6,3,5,]<-dx_rep*(finalout[3,,1,2]+finalout[3,,4,2])
  mastermat[4,3,5,]<-dx_fail*(finalout[3,,1,3]+finalout[3,,4,3])
  
  #PZA
  mastermat[6,4,5,]<-dx_rep*(finalout[4,,1,2]+finalout[4,,4,2])  
  mastermat[4,4,5,]<-dx_fail*(finalout[4,,1,3]+finalout[4,,4,3])
  
  #Rif+FQ
  mastermat[6,5,5,]<-dx_rep*(finalout[5,,1,2]+finalout[5,,4,2]) 
  mastermat[4,5,5,]<-dx_fail*(finalout[5,,1,3]+finalout[5,,4,3])
  
  #Rif+PZA
  mastermat[6,6,5,]<-dx_rep*(finalout[6,,1,2]+finalout[6,,4,2])  
  mastermat[4,6,5,]<-dx_fail*(finalout[6,,1,3]+finalout[6,,4,3]) 
  
  #FQ+PZA
  mastermat[6,7,5,]<-dx_rep*(finalout[7,,1,2]+finalout[7,,4,2])  
  mastermat[4,7,5,]<-dx_fail*(finalout[7,,1,3]+finalout[7,,4,3]) 
  
  #Rif+2
  mastermat[6,8,5,]<-dx_rep*(finalout[8,,1,2]+finalout[8,,4,2])    
  mastermat[4,8,5,]<-dx_fail*(finalout[8,,1,3]+finalout[8,,4,3]) 
  
  
  #A2->F, F->F failure: rate of dx * prop fail
  #DS
  mastermat[6,1,4,]<-dx_rep*finalout[1,,3,2]
  
  #Rif
  mastermat[6,2,4,]<-dx_rep*finalout[2,,3,2]

  #FQ
  mastermat[6,3,4,]<-dx_rep*finalout[3,,3,2]
  
  #PZA
  mastermat[6,4,4,]<-dx_rep*finalout[4,,3,2]
  
  #Rif+FQ
  mastermat[6,5,4,]<-dx_rep*finalout[5,,3,2]
  
  #Rif+PZA
  mastermat[6,6,4,]<-dx_rep*finalout[6,,3,2]
  
  #PZA+FQ
  mastermat[6,7,4,]<-dx_rep*finalout[7,,3,2]
  
  #Rif+2
  mastermat[6,8,4,]<-dx_rep*finalout[8,,3,2]
  
  
  #F->A2, A2->A2 recurrence after default and relapse
  #DS
  mastermat[4,1,6,]<-dx_fail*(finalout[1,,5,3]+finalout[1,,2,3])
  
  #Rif
  mastermat[4,2,6,]<-dx_fail*(finalout[2,,5,3]+finalout[2,,2,3])
  
  #FQ
  mastermat[4,3,6,]<-dx_fail*(finalout[3,,5,3]+finalout[3,,2,3]) 
  
  #PZA
  mastermat[4,4,6,]<-dx_fail*(finalout[4,,5,3]+finalout[4,,2,3])
    
  #Rif+FQ
  mastermat[4,5,6,]<-dx_fail*(finalout[5,,5,3]+finalout[5,,2,3])
  
  #Rif+PZA
  mastermat[4,6,6,]<-dx_fail*(finalout[6,,5,3]+finalout[6,,2,3])
  
  #PZA+FQ
  mastermat[4,7,6,]<-dx_fail*(finalout[7,,5,3]+finalout[7,,2,3])
  
  #Rif+2
  mastermat[4,8,6,]<-dx_fail*(finalout[8,,5,3]+finalout[8,,2,3])
  
  return(mastermat)
}


# # Alternative version for PaMZ
# master2<-function(era,opt,acq1,acq2,acq3,acq4,acq5,acq6,acq7,acq8,comp){ #Define function inputs
#   
#   era<-era  #Assign values to variables
#   opt<-opt
#   acq1<-acq1
#   acq2<-acq2
#   acq3<-acq3
#   acq4<-acq4
#   acq5<-acq5
#   acq6<-acq6
#   acq7<-acq7
#   acq8<-acq8
#   
#   comp<-comp
#   
#   #Create 4-dimensional array of transition rates for master matrix: 6x8x6x8
#   mastermat <- array(0, dim=c(6,8,6,8)) 
#   
#   #Dimensions 1 & 2 indicate origin compartment; 6 nat hx states x 8 drug resistance states
#   #Dimensions 3 & 4 indicate destination compartment; 6 nat hx states x 8 drug resistance states
#   
#   i <- rep(1:6,each=8) #Rows; This creates 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4, etc...
#   j <- rep(1:8,6) #Columns; This creates 1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6, etc...
#   k <- array(c(i,j),dim=c(48,2)) #This creates combo of row/column for every possible state
#   m <- array(c(k,k),dim=c(48,4)) #This creates combo for every  transition of state x-->state x (for mortality)
#   
#   ####################################################################################################################
#   ##### MORTALITY #####
#   
#   mortmat <- array(0, dim=c(6,8))  # create empty matrix
#   mortmat[,]<-m1  #Fill in with background mortality
#   mortmat[3,]<-mortmat[3,]+m2  #Add active TB mortality for naive active TB compartment
#   mortmat[6,]<-mortmat[6,]+m2  #Add active TB mortality for retreatment active TB compartment
#   
#   #No added mortality for failure compartment as assumed to have baseline mortality while receiving partially effective rx
#   #Can change this by modifying and uncommenting line below
#   #mortmat[4,]<-mortmat[4,]+0
#   
#   mortvect <-c(mortmat) #Convert mortality matrix to vector
#   mastermat[m]<-mastermat[m]+mortvect #Add mortality rates into master matrix
#   
#   ####################################################################################################################
#   ###### RESISTANCE ACQUISITION  ######
#   
#   #Switch between eras
#   #If modeling prior to intro of HRZE (era 0), model will run with acq7 values
#   #If modeling prior to intro of FQ (era 1), model will run with acq6 values
#   #If modeling post intro of resistance to all drugs incl. FQ (era 2), model will use acq1-acq5
#   acq<-array(,dim=c(8,8,5))
#   dimnames(acq)[[1]]<-dimnames(acq)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
#   dimnames(acq)[[3]]<-c("HRZE new","HRZE rerx","STR","MRZE","DST")
#   
#   if(era==0 & opt==0){
#     acq[,,1]<-acq7 #No rx, so same probs for all
#     acq[,,2]<-acq7 
#     acq[,,3]<-acq7 
#     acq[,,4]<-acq7
#     acq[,,5]<-acq7
#   }else if(era==1 & opt==0){
#     acq[,,1]<-acq6 #No resistance to FQ, so same set of probs for all
#     acq[,,2]<-acq6
#     acq[,,3]<-acq6
#     acq[,,4]<-acq6
#     acq[,,5]<-acq6
#   }else if(era==2 & opt==0){
#     acq[,,1]<-acq1 #HRZE, new
#     acq[,,2]<-acq2 #HRZE, rerx
#     acq[,,3]<-acq3 #STR
#     acq[,,4]<-acq4 #FQrx
#     acq[,,5]<-acq5 #DST
#   }else if(era==2 & opt==1){
#     acq[,,1]<-acq1 #HRZE, new
#     acq[,,2]<-acq2 #HRZE, rerx
#     acq[,,3]<-acq3 #STR
#     acq[,,4]<-acq4 #FQrx
#     acq[,,5]<-acq8 #DST
#   }
#   
#   
#   #Array of true vs. diagnosed drug resistance (conditional on availability of DST)
#   #True DR in rows, dx'd DR in columns, 3rd dimension is new vs retreatment vs. failure
#   dxmat<-array(0,dim=c(8,8,3))
#   dimnames(dxmat)[[1]]<-dimnames(dxmat)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
#   dimnames(dxmat)[[3]]<-c("New","Rep","Fail")
#   
#   #DS is always dx'd as DS (perfect specificity)
#   dxmat[1,1,]<-1
#   
#   for(i in 1:3){ #Loop to fill in values for New, Rep, Fail
#     #Rif-R always dx'd as either DS or Rif-R
#     dxmat[2,2,i]<-DSTr_sens*(dst_access[i,2,(opt+1)]+ #Rif only DST
#                                dst_access[i,5,(opt+1)]+ #Rif+FQ DST
#                                dst_access[i,6,(opt+1)]+ #Rif+PZA DST
#                                dst_access[i,8,(opt+1)]) #comprehensive DST
#     dxmat[2,1,i]<-1-dxmat[2,2,i]
#     
#     #FQ-R
#     dxmat[3,3,i]<-DSTf_sens*(dst_access[i,3,(opt+1)]
#                              +dst_access[i,5,(opt+1)]
#                              +dst_access[i,7,(opt+1)]
#                              +dst_access[i,8,(opt+1)]) #Correctly dx'd
#     dxmat[3,1,i]<-1-dxmat[3,3,i]
#     
#     #PZA-R
#     dxmat[4,4,i]<-DSTp_sens*(dst_access[i,4,(opt+1)]
#                              +dst_access[i,6,(opt+1)]
#                              +dst_access[i,7,(opt+1)]
#                              +dst_access[i,8,(opt+1)]) #Correctly dx'd
#     dxmat[4,1,i]<-1-dxmat[4,4,i]
#     
#     #Rif+FQ
#     dxmat[5,5,i]<-DSTrf_sens*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
#     dxmat[5,2,i]<-(DSTr_sens*(1-DSTf_sens))*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as Rif-R only after testing for both Rif and FQ
#       DSTr_sens*(dst_access[i,2,(opt+1)]+dst_access[i,6,(opt+1)]) #dx'd as Rif-R only after testing for Rif only or RP
#     dxmat[5,3,i]<-(DSTf_sens*(1-DSTr_sens))*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as FQ-R only after testing for both Rif and FQ
#       DSTf_sens*(dst_access[i,3,(opt+1)]+dst_access[i,7,(opt+1)]) #dx'd as FQ-R only after testing for FQ only or FP
#     dxmat[5,1,i]<-(1-DSTr_sens)*(1-DSTf_sens)*(dst_access[i,5,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both Rif and FQ
#       (1-DSTr_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS after testing for RP
#       (1-DSTf_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FP
#       (1-DSTr_sens)*dst_access[i,2,(opt+1)]+ #dx'd as DS only after testing for Rif only
#       (1-DSTf_sens)*dst_access[i,3,(opt+1)]+ #dx'd as DS only after testing for FQ only
#       dst_access[i,4,(opt+1)]+ #Received DST for P only, assumed to be DS
#       dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
#     
#     #Rif+PZA
#     dxmat[6,6,i]<-DSTrp_sens*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
#     dxmat[6,2,i]<-(DSTr_sens*(1-DSTp_sens))*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as Rif-R only after testing for both Rif and PZA
#       DSTr_sens*(dst_access[i,2,(opt+1)]+dst_access[i,5,(opt+1)]) #dx'd as Rif-R only after testing for Rif only or RF
#     dxmat[6,4,i]<-(DSTp_sens*(1-DSTr_sens))*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as PZA-R only after testing for both Rif and PZA
#       DSTp_sens*(dst_access[i,4,(opt+1)]+dst_access[i,7,(opt+1)]) #dx'd as PZA-R only after testing for PZA only or FP
#     dxmat[6,1,i]<-(1-DSTr_sens)*(1-DSTp_sens)*(dst_access[i,6,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both Rif and PZA
#       (1-DSTr_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS after testing for RF
#       (1-DSTp_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FP
#       (1-DSTr_sens)*dst_access[i,2,(opt+1)]+ #dx'd as DS only after testing for Rif only
#       (1-DSTp_sens)*dst_access[i,4,(opt+1)]+ #dx'd as DS only after testing for PZA only
#       dst_access[i,3,(opt+1)]+ #Received DST for F only, assumed to be DS
#       dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
#     
#     #FQ+PZA
#     dxmat[7,7,i]<-DSTfp_sens*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)]) #Correcty dx'd
#     dxmat[7,3,i]<-(DSTf_sens*(1-DSTp_sens))*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as FQ-R only after testing for both FQ and PZA
#       DSTf_sens*(dst_access[i,3,(opt+1)]+dst_access[i,5,(opt+1)]) #dx'd as FQ-R only after testing for FQ only or RF
#     dxmat[7,4,i]<-(DSTp_sens*(1-DSTf_sens))*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as PZA-R only after testing for both FQ and PZA
#       DSTp_sens*(dst_access[i,4,(opt+1)]+dst_access[i,6,(opt+1)]) #dx'd as PZA-R only after testing for PZA only or RP
#     dxmat[7,1,i]<-(1-DSTf_sens)*(1-DSTp_sens)*(dst_access[i,7,(opt+1)]+dst_access[i,8,(opt+1)])+ #dx'd as DS after testing for both FQ and PZA
#       (1-DSTf_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS after testing for RF
#       (1-DSTp_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS after testing for RP
#       (1-DSTf_sens)*dst_access[i,3,(opt+1)]+ #dx'd as DS only after testing for FQ only
#       (1-DSTp_sens)*dst_access[i,4,(opt+1)]+ #dx'd as DS only after testing for PZA only
#       dst_access[i,2,(opt+1)]+ #Received DST for R only, assumed to be DS
#       dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
#     
#     #Triple-resistant
#     dxmat[8,8,i]<-DSTx_sens*dst_access[i,8,(opt+1)] #Correcty dx'd
#     dxmat[8,7,i]<-DSTfp_sens*(1-DSTr_sens)*dst_access[i,8,(opt+1)]+ #dx'd as FQ+PZA after testing for all
#       DSTfp_sens*dst_access[i,7,(opt+1)] #Dx'd as FQ+PZA after testing for FQ and PZA
#     dxmat[8,6,i]<-DSTrp_sens*(1-DSTf_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif+PZA after testing for all
#       DSTrp_sens*dst_access[i,6,(opt+1)] #Dx'd as Rif+PZA after testing for Rif and PZA  
#     dxmat[8,5,i]<-DSTrf_sens*(1-DSTp_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif+FQ after testing for all
#       DSTrf_sens*dst_access[i,5,(opt+1)] #Dx'd as Rif+FQ after testing for FQ and Rif
#     dxmat[8,4,i]<-(DSTp_sens*(1-DSTr_sens)*(1-DSTf_sens))*dst_access[i,8,(opt+1)]+ #dx'd as PZA-R only after testing for all
#       (DSTp_sens*(1-DSTf_sens))*dst_access[i,7,(opt+1)]+ #dx'd as PZA-R only after testing for PZA and FQ
#       (DSTp_sens*(1-DSTr_sens))*dst_access[i,6,(opt+1)]+ #dx'd as PZA-R only after testing for PZA and Rif
#       DSTp_sens*dst_access[i,4,(opt+1)]#dx'd as PZA-R only after testing for PZA only
#     dxmat[8,3,i]<-(DSTf_sens*(1-DSTr_sens)*(1-DSTp_sens))*dst_access[i,8,(opt+1)]+ #dx'd as FQ-R only after testing for all
#       (DSTf_sens*(1-DSTp_sens))*dst_access[i,7,(opt+1)]+ #dx'd as FQ-R only after testing for PZA and FQ
#       (DSTf_sens*(1-DSTr_sens))*dst_access[i,5,(opt+1)]+ #dx'd as FQ-R only after testing for FQ and Rif
#       DSTf_sens*dst_access[i,3,(opt+1)]#dx'd as FQ-R only after testing for FQ only
#     dxmat[8,2,i]<-(DSTr_sens*(1-DSTf_sens)*(1-DSTp_sens))*dst_access[i,8,(opt+1)]+ #dx'd as Rif-R only after testing for all
#       (DSTr_sens*(1-DSTp_sens))*dst_access[i,6,(opt+1)]+ #dx'd as Rif-R only after testing for PZA and Rif
#       (DSTr_sens*(1-DSTf_sens))*dst_access[i,5,(opt+1)]+ #dx'd as Rif-R only after testing for FQ and Rif
#       DSTr_sens*dst_access[i,2,(opt+1)]#dx'd as Rif-R only after testing for Rif only
#     dxmat[8,1,i]<-(1-DSTr_sens)*(1-DSTf_sens)*(1-DSTp_sens)*dst_access[i,8,(opt+1)]+ #dx'd as Rif-R only after testing for all
#       (1-DSTf_sens)*(1-DSTp_sens)*dst_access[i,7,(opt+1)]+ #dx'd as DS after testing for FQ and PZA
#       (1-DSTr_sens)*(1-DSTp_sens)*dst_access[i,6,(opt+1)]+ #dx'd as DS only after testing for PZA and Rif
#       (1-DSTr_sens)*(1-DSTf_sens)*dst_access[i,5,(opt+1)]+ #dx'd as DS only after testing for FQ and Rif
#       (1-DSTp_sens)*dst_access[i,4,(opt+1)]+#dx'd as DS after testing for PZA only
#       (1-DSTf_sens)*dst_access[i,3,(opt+1)]+#dx'd as DS  after testing for FQ only
#       (1-DSTr_sens)*dst_access[i,2,(opt+1)]+#dx'd as DS after testing for Rif only
#       dst_access[i,1,(opt+1)] #Did not receive DST, assumed to be DS
#   }
#   
#   
#   #Check rowSums to make sure that they add up to 1
#   for(i in 1:8){
#     if(isTRUE(all.equal(rowSums(dxmat)[i],3,check.attributes = F)) !=1) stop("error in dxmat")
#   }
#   
#   ####################################################################################################################
#   ##### PROBABILITIES OF CURE #####
#   
#   #Probabilities of cure by DR state and choice of treatment regimen
#   rx<-array(0,dim=c(8,6)) #Create empty array
#   rownames(rx)<-c("DS","R","F","P","RF","RP","FP","RFP") #Rows: DR states
#   colnames(rx)<-c("HRZE","PaMZ","DST","ReRx","STR","Def") #Columns: Regimen
#   
#   #Each row = DR state, each column= (1) first-line HRZE; (2) FQ-based first-line (MRZE), 
#   #(3) individualized DST-based treatment, (4) Empiric regimen for retreatment cases ("category 2"),
#   #(5) standardized MDR regimen, (6) special values for default cases
#   
#   #Assume that empiric retreatment regimen is no more efficacious than first-line rx
#   
#   ##DS-TB
#   #Prob of cure for DS-TB assumed to be the same regardless of treatment type
#   rx[1,1:5]<-0.85/0.87  
#   
#   
#   ##Monoresistance to Rif, FQ, PZA
#   #Prob of cure for monoresistant strains on HRZE = multiple <1 of prob of cure for DS-TB
#   rx[2:4,1]<-rx[1,1]*c(0.53,1,0.76/0.86) #Sources for RR of cure: Espinal, Yee
#   
#   #Set prob for FQ-R and PZA-R as multiple <1 of prob of DS-TB
#   rx[2:4,2]<-rx[1,2]*c(1,0.9,76/86) 
#   
#   #Assume FQ-monores and PZA-monores will get fully efficacious HRZE/MHRE
#   rx[2:4,3]<-c(64/70,rep(rx[1,1],2)) #Source for prob of cure with Rif-R: Orenstein
#   
#   #Set STR regimen as equivalent to DST-based for Rif-R (STR does not contain Rif or INH)
#   #Assume STR regimen for PZA-R and FQ-R has same outcomes as HRZE
#   rx[2:4,5]<-c(rx[2,3],rx[3,1],rx[4,1])
#   
#   
#   ##Polyresistant strains; Sources: Espinal, Falzon, Yee
#   rx[5,1:3]<-c(rx[2,1],rx[3,2],rx[2,3]-0.16) #Prob of cure for RF strain on HRZE, MRZE, DST (Missing H&R, missing M&R)
#   rx[5,5]<-rx[2,3]-0.26 #Prob of cure for RF strain on STR
#   
#   rx[6,1:3]<-c(rx[2,1]*76/86,rx[4,2],64/70) #Prob of cure for RP strain on HRZE, MRZE, DST
#   rx[6,5]<-(64/70*76/86) #Prob of cure for RP strain on STR
#   
#   rx[7,1:3]<-c(rx[3,1]*76/86,rx[6,1]+0.15,rx[3,3]) #Prob of cure for FP strain on HRZE, MRZE, DST; set as multiple of outcomes for FQ-R
#   rx[7,5]<-rx[3,5]*76/86 #Prob of cure for FP strain on STR
#   
#   rx[8,1:3]<-c(rx[5,1]*76/86,rx[7,2],rx[5,3]) #Prob of cure for RFP strain on HRZE, MRZE, DST; set as multiple of outcomes for RF
#   rx[8,5]<-rx[5,5]*76/86 #Prob of cure for RFP strain on STR
#   
#   ##Retreatment regimens
#   #Set empiric retreatment regimen as similar to HRZE (if opt=0) or MRZE (if opt=1)
#   rx[,4]<-rx[,opt+1] 
#   
#   #Defaulters
#   rx[,6]<-c(0.5,0,0.25,0.25,0,0,0,0)
#   
#   ####################################################################################################################
#   ##### TREATMENT COMPLETION #####
#   
#   compmat<-array(dim=c(8,5)) #8 DR states (as diagnosed) x 5 regimens
#   rownames(compmat)<-c("DS","R","F","P","RF","RP","FP","RFP") #Rows: DR states
#   colnames(compmat)<-c("HRZE","MRZE","DST","ReRx","STR") #Columns: Regimen
#   compmat[,1:2]<-rep(comp[1],8) #first-line
#   compmat[,4]<-rep(comp[2],8) #ReRx
#   compmat[,5]<-rep(comp[3],8) #STR
#   
#   #Values for DST-based regimen               
#   compmat[,3]<-c(comp[1], #DS would get HRZE/PaMZ
#                  comp[1], #Rif-R would get PaMZ
#                  comp[1], #FQ-R would get HRZE
#                  comp[2], #PZA-R would get HRE/9 mo or similar
#                  comp[3], #Rif+FQ-R would get 2nd-line rx
#                  comp[3], #Rif+PZA-R would get 2nd-line rx
#                  comp[2], #FQ+PZA-R would get HRE/9 mo or similar
#                  comp[3]) #Triple-res would get 2nd-line rx
#   
#   
#   #Matrices of treatment completion and default based on dx'd DR, new or ReRx 
#   #8 DR states (rows) by 4 % completion (columns; HRZE, MRZE, STR, DST)
#   #Based on diagnosed DR, how likely are you to complete treatment?
#   dxrx_comp_new <- dxrx[,,(1+opt)]*compmat[,c(1,2,5,3)] #Complete rx, new cases
#   dxrx_comp_rerx <- dxrx[,,(1+opt)]*compmat[,c(4,4,5,3)] #Complete rx, retreatment
#   dxrx_def_new <- dxrx[,,(1+opt)]*(1-compmat[,c(1,2,5,3)]) #Default, new cases
#   dxrx_def_rerx <- dxrx[,,(1+opt)]*(1-compmat[,c(4,4,5,3)]) #Default, retreatment
#   
#   
#   ####################################################################################################################
#   ##### TREATMENT OUTCOMES #####
#   
#   #Treatment outcomes are conditional on DR state and choice of regimen
#   #5 possible outcomes: dc (cure after default), df (continued active TB after default), rx (cure),
#   #                     rl (relapse after initial cure), fl (treatment failure, continued active TB)
#   
#   #Variables below are named based on: (1) outcome; (2) DR state and new vs. retreatment; (3) choice of regimen
#   
#   #Note: In order to allow for uncertainty around prob of cure, each simulation includes a randomly sampled "failure multiplier"
#   #Final prob of cure = 1 - (1 -  baseline prob)*failure multiplier
#   
#   #Matrix of proportion getting each regimen
#   #Dot product of dxmat and dxrx
#   #Based on initial DR, what is your prob of cure?
#   #Will make matrix with (1) true initial DR state in rows, (2) rx in column
#   #Add 3rd dimension for new, rep, fail
#   
#   #New cases
#   rxmat<-array(,dim=c(8,4,3,2))
#   dimnames(rxmat)[[1]]<-c("DS","R","F","P","RF","RP","FP","RFP")
#   dimnames(rxmat)[[2]]<-c("HRZE","MRZE","STR","DST")
#   dimnames(rxmat)[[3]]<-c("New","Rep","Fail")
#   dimnames(rxmat)[[4]]<-c("Comp","Def") #Completers vs defaulters
#   
#   rxmat[,,1,1]<-dxmat[,,1] %*% dxrx_comp_new 
#   rxmat[,,2,1]<-dxmat[,,2] %*% dxrx_comp_rerx
#   rxmat[,,3,1]<-dxmat[,,3] %*% dxrx_comp_rerx
#   rxmat[,,1,2]<-dxmat[,,1] %*% dxrx_def_new 
#   rxmat[,,2,2]<-dxmat[,,2] %*% dxrx_def_rerx
#   rxmat[,,3,2]<-dxmat[,,3] %*% dxrx_def_rerx
#   
#   #Check sums--should add up to 1
#   for(i in 1:3){
#     for(j in 1:8){
#       if(isTRUE(all.equal(rowSums(rxmat[,,i,])[j],1,check.attributes = F)) !=1) stop("error in rxmat")
#     }
#   }
#   
#   
#   #Adjusted prob of cure (taking into account failure multiplier)
#   rx2<-1-rx
#   for(i in 1:5){ #Do this for all reg but not for defaulters
#     rx2[,i]<-rx2[,i]*unlist(dstlhs1[loop,13:20])
#   }
#   rx2<-1-rx2
#   
#   #Treatment completion
#   
#   
#   #   #Final prob of cure is dependent on choice of rx and final DR state
#   #   #For new cases, first-line regimen
#   #final outcomes array
#   #8 rows (initial DR) x 8 columns (final DR) x 5 layers (each outcome)
#   finalout<-array(0, dim=c(8,8,5,3))
#   dimnames(finalout)[[1]]<-dimnames(finalout)[[2]]<-c("DS","R","F","P","RF","RP","FP","RFP")
#   dimnames(finalout)[[3]]<-c("Cure","Relapse","Failure","Def/cure","Def/fail")
#   dimnames(finalout)[[4]]<-c("New","Rep","Fail")  
#   
#   for(h in 1:3){ 
#     tmp<-ceiling(0.5+h/3) #Will round to 1 for New, 2 for Rep and Fail
#     
#     for(i in 1:8){ #Each row of the rxmat
#       
#       #Cured
#       finalout[i,,1,h]<-
#         rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE 
#         rx2[,1]*(1-rlmat[,1])+
#         rxmat[i,3,h,1]*acq[i,,3]* #STR
#         rx2[,5]*(1-rlmat[,2])+
#         rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
#         rx2[,2]*(1-rlmat[,1])+
#         rxmat[i,4,h,1]*acq[i,,5]* #DST
#         rx2[,3]*(1-rlmat[,2])
#       
#       #Relapse
#       finalout[i,,2,h]<-
#         rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE 
#         rx2[,1]*rlmat[,1]+
#         rxmat[i,3,h,1]*acq[i,,3]* #STR
#         rx2[,5]*rlmat[,2]+
#         rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
#         rx2[,2]*rlmat[,1]+
#         rxmat[i,4,h,1]*acq[i,,5]* #DST
#         rx2[,3]*rlmat[,2]
#       
#       #Fail
#       finalout[i,,3,h]<-
#         rxmat[i,1,h,1]*acq[i,,tmp]* #HRZE new
#         (1-rx2[,1])+
#         rxmat[i,3,h,1]*acq[i,,3]* #STR
#         (1-rx2[,5])+
#         rxmat[i,2,h,1]*acq[i,,4]* #FQ 1st-line
#         (1-rx2[,2])+
#         rxmat[i,4,h,1]*acq[i,,5]* #DST
#         (1-rx2[,3])
#       
#       #Cure after default
#       finalout[i,,4,h]<-
#         (rxmat[i,1,h,2]*acq[i,,tmp]+ #HRZE new
#            rxmat[i,3,h,2]*acq[i,,3]+ #STR
#            rxmat[i,2,h,2]*acq[i,,4]+ #FQ 1st-line
#            rxmat[i,4,h,2]*acq[i,,5])* #DST
#         rx2[,6]
#       
#       #Continued active TB after default
#       finalout[i,,5,h]<-
#         (rxmat[i,1,h,2]*acq[i,,tmp]+ #HRZE new
#            rxmat[i,3,h,2]*acq[i,,3]+ #STR
#            rxmat[i,2,h,2]*acq[i,,4]+ #FQ 1st-line
#            rxmat[i,4,h,2]*acq[i,,5])* #DST
#         (1-rx2[,6])
#     } #End of i loop
#   } #End of h loop
#   
#   #Check that probabilities add up to 1
#   #   for(h in 1:3){
#   #     for(i in 1:8){
#   #   print(sum(finalout[i,,,h]))
#   #     }
#   #   }
#   
#   #FILL IN MASTERMAT NOW!!!
#   
#   ####################################################################################################################
#   ##### DISEASE PROGRESSION #####
#   
#   #Mastermat array dimensions 1 & 2 indicate ORIGIN compartment; 6 nat hx states x 8 drug resistance states
#   #Dimensions 3 & 4 indicate DESTINATION compartment; 6 nat hx states x 8 drug resistance states
#   #Nat hx states: (1) Susceptible (2) Latent rx-naive (3) Active TB x-naive 
#   #               (4) Failing (5) Latent rx-experienced (6) Active TB rx-experienced
#   #DR states: (1) DS (2) Rif only (3) FQ only (4) PZA only
#   #           (5) Rif & FQ (6) Rif & PZA (7) PZA &FQ (8) RIF & FQ & PZA
#   
#   
#   ##L->A progression from latent to active TB: endogenous reactivation
#   ##No acquisition of resistance in this transition
#   #New cases
#   mastermat[2,1,3,1]<-er #From latent rx-naive DS-TB to active rx-naive DS-TB
#   mastermat[2,2,3,2]<-er #From latent rx-naive Rif-R TB to active rx-naive Rif-R TB
#   mastermat[2,3,3,3]<-er
#   mastermat[2,4,3,4]<-er
#   mastermat[2,5,3,5]<-er
#   mastermat[2,6,3,6]<-er
#   mastermat[2,7,3,7]<-er
#   mastermat[2,8,3,8]<-er
#   
#   #Retreatment cases
#   mastermat[5,1,6,1]<-er 
#   mastermat[5,2,6,2]<-er
#   mastermat[5,3,6,3]<-er
#   mastermat[5,4,6,4]<-er
#   mastermat[5,5,6,5]<-er
#   mastermat[5,6,6,6]<-er
#   mastermat[5,7,6,7]<-er
#   mastermat[5,8,6,8]<-er 
#   
#   ##A->L spontaneous cure w/o treatment
#   ##No acquisition of resistant in this transition
#   #New cases
#   mastermat[3,1,2,1]<-sc #From active rx-naive DS-TB to latent rx-naive DS-TB
#   mastermat[3,2,2,2]<-sc #From active rx-naive Rif-R TB to latent rx-naive Rif-R TB
#   mastermat[3,3,2,3]<-sc
#   mastermat[3,4,2,4]<-sc
#   mastermat[3,5,2,5]<-sc
#   mastermat[3,6,2,6]<-sc
#   mastermat[3,7,2,7]<-sc
#   mastermat[3,8,2,8]<-sc
#   
#   #Retreatment cases
#   mastermat[6,1,5,1]<-sc 
#   mastermat[6,2,5,2]<-sc
#   mastermat[6,3,5,3]<-sc
#   mastermat[6,4,5,4]<-sc
#   mastermat[6,5,5,5]<-sc
#   mastermat[6,6,5,6]<-sc
#   mastermat[6,7,5,7]<-sc
#   mastermat[6,8,5,8]<-sc
#   
#   ###Initial treatment
#   ##A1->L2: successful 1st treatment; resistance can be acquired during treatment, such that
#   #         transition rate = (dx rate) * (prob of res acquisition) * (prop rx success + prop cure after default)
#   
#   #DS
#   mastermat[3,1,5,]<-dx * (finalout[1,,1,1]+finalout[1,,4,1])
#   
#   #Rif monores
#   # Rate is computed in the same manner as above but requires accounting for missed diagnosis of drug resistance
#   # Prob of res acq and prob of cure are weighted averages for 
#   #(1) those correclty dx'd as Rif-R and (2) those misdx'd as DS-TB
#   mastermat[3,2,5,]<-dx * (finalout[2,,1,1]+finalout[2,,4,1])
#   
#   #FQ monores
#   mastermat[3,3,5,]<-dx *(finalout[3,,1,1]+finalout[3,,4,1]) 
#   
#   #PZA
#   mastermat[3,4,5,]<-dx * (finalout[4,,1,1]+finalout[4,,4,1])
#   
#   #Rif+FQ
#   mastermat[3,5,5,]<-dx * (finalout[5,,1,1]+finalout[5,,4,1])
#   
#   #Rif+PZA
#   mastermat[3,6,5,]<-dx * (finalout[6,,1,1]+finalout[6,,4,1])
#   
#   #FQ+PZA
#   mastermat[3,7,5,]<-dx * (finalout[7,,1,1]+finalout[7,,4,1])
#   
#   #Rif+2
#   mastermat[3,8,5,]<-dx * (finalout[8,,1,1]+finalout[8,,4,1])
#   
#   
#   ##A1->F failure after initial treatment: rate of dx * prop fail
#   #DS
#   mastermat[3,1,4,]<-dx*finalout[1,,3,1]
#   
#   #Rif
#   mastermat[3,2,4,]<-dx*finalout[2,,3,1]
#   
#   #FQ
#   mastermat[3,3,4,]<-dx*finalout[3,,3,1]
#   
#   #PZA
#   mastermat[3,4,4,]<-dx*finalout[4,,3,1]
#   
#   #Rif+FQ
#   mastermat[3,5,4,]<-dx*finalout[5,,3,1]
#   
#   #Rif+PZA
#   mastermat[3,6,4,]<-dx*finalout[6,,3,1]
#   
#   #PZA+FQ
#   mastermat[3,7,4,]<-dx*finalout[7,,3,1]
#   
#   #Rif+2
#   mastermat[3,8,4,]<-dx*finalout[8,,3,1]
#   
#   
#   ##A1->A2 recurrence after default and Relapse
#   #DS
#   mastermat[3,1,6,]<-dx*(finalout[1,,5,1]+finalout[1,,2,1])
#   
#   #Rif
#   mastermat[3,2,6,]<-dx*(finalout[2,,5,1]+finalout[2,,2,1])
#   
#   #FQ
#   mastermat[3,3,6,]<-dx*(finalout[3,,5,1]+finalout[3,,2,1])
#   
#   #PZA
#   mastermat[3,4,6,]<-dx*(finalout[4,,5,1]+finalout[4,,2,1])
#   
#   #Rif+FQ
#   mastermat[3,5,6,]<-dx*(finalout[5,,5,1]+finalout[5,,2,1])
#   
#   #Rif+PZA
#   mastermat[3,6,6,]<-dx*(finalout[6,,5,1]+finalout[6,,2,1])
#   
#   #PZA+FQ
#   mastermat[3,7,6,]<-dx*(finalout[7,,5,1]+finalout[7,,2,1])
#   
#   #Rif+2
#   mastermat[3,8,6,]<-dx*(finalout[8,,5,1]+finalout[8,,2,1])
#   
#   
#   ##RETREATMENT 
#   #A2->L2; F->L2 successful retreatment: rate of dx * (prop treatment success+prop cure after default)
#   #DS
#   mastermat[6,1,5,]<-dx_rep*(finalout[1,,1,2]+finalout[1,,4,2])
#   mastermat[4,1,5,]<-dx_fail*(finalout[1,,1,3]+finalout[1,,4,3])
#   
#   #Rif
#   mastermat[6,2,5,]<-dx_rep*(finalout[2,,1,2]+finalout[2,,4,2])
#   mastermat[4,2,5,]<-dx_fail*(finalout[2,,1,3]+finalout[2,,4,3])
#   
#   #FQ
#   mastermat[6,3,5,]<-dx_rep*(finalout[3,,1,2]+finalout[3,,4,2])
#   mastermat[4,3,5,]<-dx_fail*(finalout[3,,1,3]+finalout[3,,4,3])
#   
#   #PZA
#   mastermat[6,4,5,]<-dx_rep*(finalout[4,,1,2]+finalout[4,,4,2])  
#   mastermat[4,4,5,]<-dx_fail*(finalout[4,,1,3]+finalout[4,,4,3])
#   
#   #Rif+FQ
#   mastermat[6,5,5,]<-dx_rep*(finalout[5,,1,2]+finalout[5,,4,2]) 
#   mastermat[4,5,5,]<-dx_fail*(finalout[5,,1,3]+finalout[5,,4,3])
#   
#   #Rif+PZA
#   mastermat[6,6,5,]<-dx_rep*(finalout[6,,1,2]+finalout[6,,4,2])  
#   mastermat[4,6,5,]<-dx_fail*(finalout[6,,1,3]+finalout[6,,4,3]) 
#   
#   #FQ+PZA
#   mastermat[6,7,5,]<-dx_rep*(finalout[7,,1,2]+finalout[7,,4,2])  
#   mastermat[4,7,5,]<-dx_fail*(finalout[7,,1,3]+finalout[7,,4,3]) 
#   
#   #Rif+2
#   mastermat[6,8,5,]<-dx_rep*(finalout[8,,1,2]+finalout[8,,4,2])    
#   mastermat[4,8,5,]<-dx_fail*(finalout[8,,1,3]+finalout[8,,4,3]) 
#   
#   
#   #A2->F, F->F failure: rate of dx * prop fail
#   #DS
#   mastermat[6,1,4,]<-dx_rep*finalout[1,,3,2]
#   
#   #Rif
#   mastermat[6,2,4,]<-dx_rep*finalout[2,,3,2]
#   
#   #FQ
#   mastermat[6,3,4,]<-dx_rep*finalout[3,,3,2]
#   
#   #PZA
#   mastermat[6,4,4,]<-dx_rep*finalout[4,,3,2]
#   
#   #Rif+FQ
#   mastermat[6,5,4,]<-dx_rep*finalout[5,,3,2]
#   
#   #Rif+PZA
#   mastermat[6,6,4,]<-dx_rep*finalout[6,,3,2]
#   
#   #PZA+FQ
#   mastermat[6,7,4,]<-dx_rep*finalout[7,,3,2]
#   
#   #Rif+2
#   mastermat[6,8,4,]<-dx_rep*finalout[8,,3,2]
#   
#   
#   #F->A2, A2->A2 recurrence after default and relapse
#   #DS
#   mastermat[4,1,6,]<-dx_fail*(finalout[1,,5,3]+finalout[1,,2,3])
#   
#   #Rif
#   mastermat[4,2,6,]<-dx_fail*(finalout[2,,5,3]+finalout[2,,2,3])
#   
#   #FQ
#   mastermat[4,3,6,]<-dx_fail*(finalout[3,,5,3]+finalout[3,,2,3]) 
#   
#   #PZA
#   mastermat[4,4,6,]<-dx_fail*(finalout[4,,5,3]+finalout[4,,2,3])
#   
#   #Rif+FQ
#   mastermat[4,5,6,]<-dx_fail*(finalout[5,,5,3]+finalout[5,,2,3])
#   
#   #Rif+PZA
#   mastermat[4,6,6,]<-dx_fail*(finalout[6,,5,3]+finalout[6,,2,3])
#   
#   #PZA+FQ
#   mastermat[4,7,6,]<-dx_fail*(finalout[7,,5,3]+finalout[7,,2,3])
#   
#   #Rif+2
#   mastermat[4,8,6,]<-dx_fail*(finalout[8,,5,3]+finalout[8,,2,3])
#   
#   return(mastermat)
# }
# 


#######################################################
#FUNCTION II: COMPUTE RATES OF CHANGE AT EACH TIMESTEP#
#######################################################

#Output from this function is rates of change, to feed into ODE solver

TBdx<-function(t, state, parameters) { #function arguments are (1) timesteps, (2) initial state parameters, (3) master transition matrix
  x<-as.list(c(state,parameters))
  with(x,{
    statev <- unlist(x[1:48]) #vector of state parameters
    
    #vectors for state transitions
    masterv0 <- unlist(x[49:(49+48*48-1)])
    masterv1 <- unlist(x[(49+48*48):(49+48*48*2-1)])
    masterv2 <- unlist(x[(49+48*48*2):(49+48*48*3-1)])
    
    masterv2b <- unlist(x[(49+48*48*3):(49+48*48*4-1)])    
    
    time1<-unlist(x[(49+48*48*4)]) #time of intro of Rif resistance
    time2<-unlist(x[(49+48*48*4+1)]) #time of introduction of FQ resistance
    
    # force of infection, density-dependent, for each strain
    # (transmission parameter) * (total # in active TB and failing rx compartments for each strain)
    force_s <- tr*(statev[3]+statev[6]+relinf*statev[4]) #DS
    force_r <- tr*(statev[9]+statev[12]+relinf*statev[10])*fit_r #Rif-R
    force_f <- tr*(statev[15]+statev[18]+relinf*statev[16])*fit_f #FQ-R
    force_p <- tr*(statev[21]+statev[24]+relinf*statev[22])*fit_p #PZA-R
    force_rf <- tr*(statev[27]+statev[30]+relinf*statev[28])*fit_rf #Rif+FQ
    force_rp <- tr*(statev[33]+statev[36]+relinf*statev[34])*fit_rp #Rif+PZA
    force_fp <- tr*(statev[39]+statev[42]+relinf*statev[40])*fit_fp #FQ+PZA
    force_r2 <- tr*(statev[45]+statev[48]+relinf*statev[46])*fit_r2 #Rif+FQ+PZA
    
    # Now, create an "infection matrix"
    infectmat <- array(0, dim=c(6,8,6,8))
    
    #S->L1
    infectmat[1,1,2,1]<-force_s*(1-fp) #Progression from susceptible state to latent infection with DS-TB
    infectmat[1,1,2,2]<-force_r*(1-fp) #Progression from susceptible state to latent infection with Rif-res TB
    infectmat[1,1,2,3]<-force_f*(1-fp) #Progression from susceptible state to latent infection with FQ-res TB
    infectmat[1,1,2,4]<-force_p*(1-fp) #Progression from susceptible state to latent infection with PZA-res TB
    infectmat[1,1,2,5]<-force_rf*(1-fp) #Progression from susceptible state to latent infection with Rif+FQ res
    infectmat[1,1,2,6]<-force_rp*(1-fp) #Progression from susceptible state to latent infection with Rif+PZA res
    infectmat[1,1,2,7]<-force_fp*(1-fp) #Progression from susceptible state to latent infection with FQ+PZA res
    infectmat[1,1,2,8]<-force_r2*(1-fp) #Progression from susceptible state to latent infection with Rif+2 res
    
    #S->A1
    infectmat[1,1,3,1]<-force_s*fp #Rapid progression from susceptible state to active DS-TB
    infectmat[1,1,3,2]<-force_r*fp #Rapid progression from susceptible state to active FQ-res TB
    infectmat[1,1,3,3]<-force_f*fp #Rapid progression from susceptible state to active Rif-res TB
    infectmat[1,1,3,4]<-force_p*fp #Rapid progressionfrom susceptible state to active PZA-res TB
    infectmat[1,1,3,5]<-force_rf*fp #Rapid progression from susceptible state to active Rif+FQ res
    infectmat[1,1,3,6]<-force_rp*fp #Rapid progression from susceptible state to active Rif+PZA res
    infectmat[1,1,3,7]<-force_fp*fp #Rapid progression from susceptible state to active FQ+PZA res
    infectmat[1,1,3,8]<-force_r2*fp #Rapid progression from susceptible state to active Rif+2 res
    
    #L1->L2
    infectmat[2,1,2,2]<-force_r*lp*(1-fp)*cross[1,2] #Reinfection with latent FQ-res TB from latent DS-TB
    infectmat[2,1,2,3]<-force_f*lp*(1-fp)*cross[1,3] #Reinfection with latent Rif-res TB from latent DS-TB
    infectmat[2,1,2,4]<-force_p*lp*(1-fp)*cross[1,4] #Reinfection with latent PZA-res TB from latent DS-TB
    infectmat[2,1,2,5]<-force_rf*lp*(1-fp)*cross[1,5] #Reinfection with latent Rif+FQ TB from latent DS-TB
    infectmat[2,1,2,6]<-force_rp*lp*(1-fp)*cross[1,6] #Reinfection with latent Rif+PZA TB from latent DS-TB
    infectmat[2,1,2,7]<-force_fp*lp*(1-fp)*cross[1,7] #Reinfection with latent FQ+PZA TB from latent DS-TB
    infectmat[2,1,2,8]<-force_r2*lp*(1-fp)*cross[1,8] #Reinfection with latent Rif+2 TB from latent DS-TB
    
    infectmat[2,2,2,1]<-force_s*lp*(1-fp)*cross[2,1] #Reinfection with latent DS-TB from latent Rif-R TB
    infectmat[2,2,2,3]<-force_f*lp*(1-fp)*cross[2,3] 
    infectmat[2,2,2,4]<-force_p*lp*(1-fp)*cross[2,4]
    infectmat[2,2,2,5]<-force_rf*lp*(1-fp)*cross[2,5] 
    infectmat[2,2,2,6]<-force_rp*lp*(1-fp)*cross[2,6]
    infectmat[2,2,2,7]<-force_fp*lp*(1-fp)*cross[2,7]
    infectmat[2,2,2,8]<-force_r2*lp*(1-fp)*cross[2,8] 
    
    infectmat[2,3,2,1]<-force_s*lp*(1-fp)*cross[3,1] #Reinfection with latent DS-TB from latent FQ-R TB
    infectmat[2,3,2,2]<-force_r*lp*(1-fp)*cross[3,2] 
    infectmat[2,3,2,4]<-force_p*lp*(1-fp)*cross[3,4]
    infectmat[2,3,2,5]<-force_rf*lp*(1-fp)*cross[3,5] 
    infectmat[2,3,2,6]<-force_rp*lp*(1-fp)*cross[3,6]
    infectmat[2,3,2,7]<-force_fp*lp*(1-fp)*cross[3,7]
    infectmat[2,3,2,8]<-force_r2*lp*(1-fp)*cross[3,8] 
    
    infectmat[2,4,2,1]<-force_s*lp*(1-fp)*cross[4,1] #Reinfection with latent DS-TB from latent PZA-R TB
    infectmat[2,4,2,2]<-force_r*lp*(1-fp)*cross[4,2] 
    infectmat[2,4,2,3]<-force_f*lp*(1-fp)*cross[4,3]
    infectmat[2,4,2,5]<-force_rf*lp*(1-fp)*cross[4,5] 
    infectmat[2,4,2,6]<-force_rp*lp*(1-fp)*cross[4,6]
    infectmat[2,4,2,7]<-force_fp*lp*(1-fp)*cross[4,7]
    infectmat[2,4,2,8]<-force_r2*lp*(1-fp)*cross[4,8]
    
    infectmat[2,5,2,1]<-force_s*lp*(1-fp)*cross[5,1] #Reinfection with latent DS-TB from latent Rif+FQ TB
    infectmat[2,5,2,2]<-force_r*lp*(1-fp)*cross[5,2] 
    infectmat[2,5,2,3]<-force_f*lp*(1-fp)*cross[5,3]
    infectmat[2,5,2,4]<-force_p*lp*(1-fp)*cross[5,4] 
    infectmat[2,5,2,6]<-force_rp*lp*(1-fp)*cross[5,6]
    infectmat[2,5,2,7]<-force_fp*lp*(1-fp)*cross[5,7]
    infectmat[2,5,2,8]<-force_r2*lp*(1-fp)*cross[5,8] 
    
    infectmat[2,6,2,1]<-force_s*lp*(1-fp)*cross[6,1] #Reinfection with latent DS-TB from latent Rif+PZA TB
    infectmat[2,6,2,2]<-force_r*lp*(1-fp)*cross[6,2] 
    infectmat[2,6,2,3]<-force_f*lp*(1-fp)*cross[6,3]
    infectmat[2,6,2,4]<-force_p*lp*(1-fp)*cross[6,4] 
    infectmat[2,6,2,5]<-force_rf*lp*(1-fp)*cross[6,5]
    infectmat[2,6,2,7]<-force_fp*lp*(1-fp)*cross[6,7]
    infectmat[2,6,2,8]<-force_r2*lp*(1-fp)*cross[6,8] 
    
    infectmat[2,7,2,1]<-force_s*lp*(1-fp)*cross[7,1] #Reinfection with latent DS-TB from latent FQ+PZA TB
    infectmat[2,7,2,2]<-force_r*lp*(1-fp)*cross[7,2]  
    infectmat[2,7,2,3]<-force_f*lp*(1-fp)*cross[7,3] 
    infectmat[2,7,2,4]<-force_p*lp*(1-fp)*cross[7,4]  
    infectmat[2,7,2,5]<-force_rf*lp*(1-fp)*cross[7,5] 
    infectmat[2,7,2,6]<-force_rp*lp*(1-fp)*cross[7,6] 
    infectmat[2,7,2,8]<-force_r2*lp*(1-fp)*cross[7,8]  
    
    infectmat[2,8,2,1]<-force_s*lp*(1-fp)*cross[8,1] #Reinfection with latent DS-TB from latent triple-res TB
    infectmat[2,8,2,2]<-force_r*lp*(1-fp)*cross[8,2]  
    infectmat[2,8,2,3]<-force_f*lp*(1-fp)*cross[8,3] 
    infectmat[2,8,2,4]<-force_p*lp*(1-fp)*cross[8,4]  
    infectmat[2,8,2,5]<-force_rf*lp*(1-fp)*cross[8,5] 
    infectmat[2,8,2,6]<-force_rp*lp*(1-fp)*cross[8,6] 
    infectmat[2,8,2,7]<-force_fp*lp*(1-fp)*cross[8,7]  
    
    #L1->A1
    infectmat[2,1,3,2]<-force_r*lp*fp*cross[1,2] #Reinfection with latent FQ-res TB from latent DS-TB
    infectmat[2,1,3,3]<-force_f*lp*fp*cross[1,3] #Reinfection with latent Rif-res TB from latent DS-TB
    infectmat[2,1,3,4]<-force_p*lp*fp*cross[1,4] #Reinfection with latent PZA-res TB from latent DS-TB
    infectmat[2,1,3,5]<-force_rf*lp*fp*cross[1,5] #Reinfection with latent Rif+FQ TB from latent DS-TB
    infectmat[2,1,3,6]<-force_rp*lp*fp*cross[1,6] #Reinfection with latent Rif+PZA TB from latent DS-TB
    infectmat[2,1,3,7]<-force_fp*lp*fp*cross[1,7] #Reinfection with latent FQ+PZA TB from latent DS-TB
    infectmat[2,1,3,8]<-force_r2*lp*fp*cross[1,8] #Reinfection with latent Rif+2 TB from latent DS-TB
    
    infectmat[2,2,3,1]<-force_s*lp*fp*cross[2,1] #Reinfection with latent DS-TB from latent Rif-R TB
    infectmat[2,2,3,3]<-force_f*lp*fp*cross[2,3] 
    infectmat[2,2,3,4]<-force_p*lp*fp*cross[2,4]
    infectmat[2,2,3,5]<-force_rf*lp*fp*cross[2,5] 
    infectmat[2,2,3,6]<-force_rp*lp*fp*cross[2,6]
    infectmat[2,2,3,7]<-force_fp*lp*fp*cross[2,7]
    infectmat[2,2,3,8]<-force_r2*lp*fp*cross[2,8] 
    
    infectmat[2,3,3,1]<-force_s*lp*fp*cross[3,1] #Reinfection with latent DS-TB from latent FQ-R TB
    infectmat[2,3,3,2]<-force_r*lp*fp*cross[3,2] 
    infectmat[2,3,3,4]<-force_p*lp*fp*cross[3,4]
    infectmat[2,3,3,5]<-force_rf*lp*fp*cross[3,5] 
    infectmat[2,3,3,6]<-force_rp*lp*fp*cross[3,6]
    infectmat[2,3,3,7]<-force_fp*lp*fp*cross[3,7]
    infectmat[2,3,3,8]<-force_r2*lp*fp*cross[3,8] 
    
    infectmat[2,4,3,1]<-force_s*lp*fp*cross[4,1] #Reinfection with latent DS-TB from latent PZA-R TB
    infectmat[2,4,3,2]<-force_r*lp*fp*cross[4,2] 
    infectmat[2,4,3,3]<-force_f*lp*fp*cross[4,3]
    infectmat[2,4,3,5]<-force_rf*lp*fp*cross[4,5] 
    infectmat[2,4,3,6]<-force_rp*lp*fp*cross[4,6]
    infectmat[2,4,3,7]<-force_fp*lp*fp*cross[4,7]
    infectmat[2,4,3,8]<-force_r2*lp*fp*cross[4,8]
    
    infectmat[2,5,3,1]<-force_s*lp*fp*cross[5,1] #Reinfection with latent DS-TB from latent Rif+FQ TB
    infectmat[2,5,3,2]<-force_r*lp*fp*cross[5,2] 
    infectmat[2,5,3,3]<-force_f*lp*fp*cross[5,3]
    infectmat[2,5,3,4]<-force_p*lp*fp*cross[5,4] 
    infectmat[2,5,3,6]<-force_rp*lp*fp*cross[5,6]
    infectmat[2,5,3,7]<-force_fp*lp*fp*cross[5,7]
    infectmat[2,5,3,8]<-force_r2*lp*fp*cross[5,8] 
    
    infectmat[2,6,3,1]<-force_s*lp*fp*cross[6,1] #Reinfection with latent DS-TB from latent Rif+PZA TB
    infectmat[2,6,3,2]<-force_r*lp*fp*cross[6,2] 
    infectmat[2,6,3,3]<-force_f*lp*fp*cross[6,3]
    infectmat[2,6,3,4]<-force_p*lp*fp*cross[6,4] 
    infectmat[2,6,3,5]<-force_rf*lp*fp*cross[6,5]
    infectmat[2,6,3,7]<-force_fp*lp*fp*cross[6,7]
    infectmat[2,6,3,8]<-force_r2*lp*fp*cross[6,8] 
    
    infectmat[2,7,3,1]<-force_s*lp*fp*cross[7,1] #Reinfection with latent DS-TB from latent FQ+PZA TB
    infectmat[2,7,3,2]<-force_r*lp*fp*cross[7,2]  
    infectmat[2,7,3,3]<-force_f*lp*fp*cross[7,3] 
    infectmat[2,7,3,4]<-force_p*lp*fp*cross[7,4]  
    infectmat[2,7,3,5]<-force_rf*lp*fp*cross[7,5] 
    infectmat[2,7,3,6]<-force_rp*lp*fp*cross[7,6] 
    infectmat[2,7,3,8]<-force_r2*lp*fp*cross[7,8]  
    
    infectmat[2,8,3,1]<-force_s*lp*fp*cross[8,1] #Reinfection with latent DS-TB from latent triple-res TB
    infectmat[2,8,3,2]<-force_r*lp*fp*cross[8,2]  
    infectmat[2,8,3,3]<-force_f*lp*fp*cross[8,3] 
    infectmat[2,8,3,4]<-force_p*lp*fp*cross[8,4]  
    infectmat[2,8,3,5]<-force_rf*lp*fp*cross[8,5] 
    infectmat[2,8,3,6]<-force_rp*lp*fp*cross[8,6] 
    infectmat[2,8,3,7]<-force_fp*lp*fp*cross[8,7]  
    
    #L2->L2
    infectmat[5,1,5,2]<-force_r*lp*(1-fp)*cross[1,2] #Reinfection with latent FQ-res TB from latent DS-TB
    infectmat[5,1,5,3]<-force_f*lp*(1-fp)*cross[1,3] #Reinfection with latent Rif-res TB from latent DS-TB
    infectmat[5,1,5,4]<-force_p*lp*(1-fp)*cross[1,4] #Reinfection with latent PZA-res TB from latent DS-TB
    infectmat[5,1,5,5]<-force_rf*lp*(1-fp)*cross[1,5] #Reinfection with latent Rif+FQ TB from latent DS-TB
    infectmat[5,1,5,6]<-force_rp*lp*(1-fp)*cross[1,6] #Reinfection with latent Rif+PZA TB from latent DS-TB
    infectmat[5,1,5,7]<-force_fp*lp*(1-fp)*cross[1,7] #Reinfection with latent FQ+PZA TB from latent DS-TB
    infectmat[5,1,5,8]<-force_r2*lp*(1-fp)*cross[1,8] #Reinfection with latent Rif+2 TB from latent DS-TB
    
    infectmat[5,2,5,1]<-force_s*lp*(1-fp)*cross[2,1] #Reinfection with latent DS-TB from latent Rif-R TB
    infectmat[5,2,5,3]<-force_f*lp*(1-fp)*cross[2,3] 
    infectmat[5,2,5,4]<-force_p*lp*(1-fp)*cross[2,4]
    infectmat[5,2,5,5]<-force_rf*lp*(1-fp)*cross[2,5] 
    infectmat[5,2,5,6]<-force_rp*lp*(1-fp)*cross[2,6]
    infectmat[5,2,5,7]<-force_fp*lp*(1-fp)*cross[2,7]
    infectmat[5,2,5,8]<-force_r2*lp*(1-fp)*cross[2,8] 
    
    infectmat[5,3,5,1]<-force_s*lp*(1-fp)*cross[3,1] #Reinfection with latent DS-TB from latent FQ-R TB
    infectmat[5,3,5,2]<-force_r*lp*(1-fp)*cross[3,2] 
    infectmat[5,3,5,4]<-force_p*lp*(1-fp)*cross[3,4]
    infectmat[5,3,5,5]<-force_rf*lp*(1-fp)*cross[3,5] 
    infectmat[5,3,5,6]<-force_rp*lp*(1-fp)*cross[3,6]
    infectmat[5,3,5,7]<-force_fp*lp*(1-fp)*cross[3,7]
    infectmat[5,3,5,8]<-force_r2*lp*(1-fp)*cross[3,8] 
    
    infectmat[5,4,5,1]<-force_s*lp*(1-fp)*cross[4,1] #Reinfection with latent DS-TB from latent PZA-R TB
    infectmat[5,4,5,2]<-force_r*lp*(1-fp)*cross[4,2] 
    infectmat[5,4,5,3]<-force_f*lp*(1-fp)*cross[4,3]
    infectmat[5,4,5,5]<-force_rf*lp*(1-fp)*cross[4,5] 
    infectmat[5,4,5,6]<-force_rp*lp*(1-fp)*cross[4,6]
    infectmat[5,4,5,7]<-force_fp*lp*(1-fp)*cross[4,7]
    infectmat[5,4,5,8]<-force_r2*lp*(1-fp)*cross[4,8]
    
    infectmat[5,5,5,1]<-force_s*lp*(1-fp)*cross[5,1] #Reinfection with latent DS-TB from latent Rif+FQ TB
    infectmat[5,5,5,2]<-force_r*lp*(1-fp)*cross[5,2] 
    infectmat[5,5,5,3]<-force_f*lp*(1-fp)*cross[5,3]
    infectmat[5,5,5,4]<-force_p*lp*(1-fp)*cross[5,4] 
    infectmat[5,5,5,6]<-force_rp*lp*(1-fp)*cross[5,6]
    infectmat[5,5,5,7]<-force_fp*lp*(1-fp)*cross[5,7]
    infectmat[5,5,5,8]<-force_r2*lp*(1-fp)*cross[5,8] 
    
    infectmat[5,6,5,1]<-force_s*lp*(1-fp)*cross[6,1] #Reinfection with latent DS-TB from latent Rif+PZA TB
    infectmat[5,6,5,2]<-force_r*lp*(1-fp)*cross[6,2] 
    infectmat[5,6,5,3]<-force_f*lp*(1-fp)*cross[6,3]
    infectmat[5,6,5,4]<-force_p*lp*(1-fp)*cross[6,4] 
    infectmat[5,6,5,5]<-force_rf*lp*(1-fp)*cross[6,5]
    infectmat[5,6,5,7]<-force_fp*lp*(1-fp)*cross[6,7]
    infectmat[5,6,5,8]<-force_r2*lp*(1-fp)*cross[6,8] 
    
    infectmat[5,7,5,1]<-force_s*lp*(1-fp)*cross[7,1] #Reinfection with latent DS-TB from latent FQ+PZA TB
    infectmat[5,7,5,2]<-force_r*lp*(1-fp)*cross[7,2]  
    infectmat[5,7,5,3]<-force_f*lp*(1-fp)*cross[7,3] 
    infectmat[5,7,5,4]<-force_p*lp*(1-fp)*cross[7,4]  
    infectmat[5,7,5,5]<-force_rf*lp*(1-fp)*cross[7,5] 
    infectmat[5,7,5,6]<-force_rp*lp*(1-fp)*cross[7,6] 
    infectmat[5,7,5,8]<-force_r2*lp*(1-fp)*cross[7,8]  
    
    infectmat[5,8,5,1]<-force_s*lp*(1-fp)*cross[8,1] #Reinfection with latent DS-TB from latent triple-res TB
    infectmat[5,8,5,2]<-force_r*lp*(1-fp)*cross[8,2]  
    infectmat[5,8,5,3]<-force_f*lp*(1-fp)*cross[8,3] 
    infectmat[5,8,5,4]<-force_p*lp*(1-fp)*cross[8,4]  
    infectmat[5,8,5,5]<-force_rf*lp*(1-fp)*cross[8,5] 
    infectmat[5,8,5,6]<-force_rp*lp*(1-fp)*cross[8,6] 
    infectmat[5,8,5,7]<-force_fp*lp*(1-fp)*cross[8,7]  
    
    #L2->A2
    infectmat[5,1,6,2]<-force_r*lp*fp*cross[1,2] #Reinfection with latent DS-TB from latent Rif-R TB
    infectmat[5,1,6,3]<-force_f*lp*fp*cross[1,3] 
    infectmat[5,1,6,4]<-force_p*lp*fp*cross[1,4]
    infectmat[5,1,6,5]<-force_rf*lp*fp*cross[1,5] 
    infectmat[5,1,6,6]<-force_rp*lp*fp*cross[1,6]
    infectmat[5,1,6,7]<-force_fp*lp*fp*cross[1,7]
    infectmat[5,1,6,8]<-force_r2*lp*fp*cross[1,8] 
    
    infectmat[5,2,6,1]<-force_s*lp*fp*cross[2,1] #Reinfection with latent DS-TB from latent FQ-R TB
    infectmat[5,2,6,3]<-force_f*lp*fp*cross[2,3] 
    infectmat[5,2,6,4]<-force_p*lp*fp*cross[2,4]
    infectmat[5,2,6,5]<-force_rf*lp*fp*cross[2,5] 
    infectmat[5,2,6,6]<-force_rp*lp*fp*cross[2,6]
    infectmat[5,2,6,7]<-force_fp*lp*fp*cross[2,7]
    infectmat[5,2,6,8]<-force_r2*lp*fp*cross[2,8] 
    
    infectmat[5,3,6,1]<-force_s*lp*fp*cross[3,1] #Reinfection with latent DS-TB from latent PZA-R TB
    infectmat[5,3,6,2]<-force_r*lp*fp*cross[3,2] 
    infectmat[5,3,6,4]<-force_p*lp*fp*cross[3,4]
    infectmat[5,3,6,5]<-force_rf*lp*fp*cross[3,5] 
    infectmat[5,3,6,6]<-force_rp*lp*fp*cross[3,6]
    infectmat[5,3,6,7]<-force_fp*lp*fp*cross[3,7]
    infectmat[5,3,6,8]<-force_r2*lp*fp*cross[3,8] 
    
    infectmat[5,4,6,1]<-force_s*lp*fp*cross[4,1] #Reinfection with latent DS-TB from latent Rif+FQ TB
    infectmat[5,4,6,2]<-force_r*lp*fp*cross[4,2] 
    infectmat[5,4,6,3]<-force_f*lp*fp*cross[4,3]
    infectmat[5,4,6,5]<-force_rf*lp*fp*cross[4,5] 
    infectmat[5,4,6,6]<-force_rp*lp*fp*cross[4,6]
    infectmat[5,4,6,7]<-force_fp*lp*fp*cross[4,7]
    infectmat[5,4,6,8]<-force_r2*lp*fp*cross[4,8]
    
    infectmat[5,5,6,1]<-force_s*lp*fp*cross[5,1] #Reinfection with latent DS-TB from latent Rif+PZA TB
    infectmat[5,5,6,2]<-force_r*lp*fp*cross[5,2] 
    infectmat[5,5,6,3]<-force_f*lp*fp*cross[5,3]
    infectmat[5,5,6,4]<-force_p*lp*fp*cross[5,4] 
    infectmat[5,5,6,6]<-force_rp*lp*fp*cross[5,6]
    infectmat[5,5,6,7]<-force_fp*lp*fp*cross[5,7]
    infectmat[5,5,6,8]<-force_r2*lp*fp*cross[5,8] 
    
    infectmat[5,6,6,1]<-force_s*lp*fp*cross[6,1] #Reinfection with latent DS-TB from latent Rif+PZA TB
    infectmat[5,6,6,2]<-force_r*lp*fp*cross[6,2] 
    infectmat[5,6,6,3]<-force_f*lp*fp*cross[6,3]
    infectmat[5,6,6,4]<-force_p*lp*fp*cross[6,4] 
    infectmat[5,6,6,5]<-force_rf*lp*fp*cross[6,5]
    infectmat[5,6,6,7]<-force_fp*lp*fp*cross[6,7]
    infectmat[5,6,6,8]<-force_r2*lp*fp*cross[6,8] 
    
    infectmat[5,7,6,1]<-force_s*lp*fp*cross[7,1] #Reinfection with latent DS-TB from latent FQ+PZA TB
    infectmat[5,7,6,2]<-force_r*lp*fp*cross[7,2]  
    infectmat[5,7,6,3]<-force_f*lp*fp*cross[7,3] 
    infectmat[5,7,6,4]<-force_p*lp*fp*cross[7,4]  
    infectmat[5,7,6,5]<-force_rf*lp*fp*cross[7,5] 
    infectmat[5,7,6,6]<-force_rp*lp*fp*cross[7,6] 
    infectmat[5,7,6,8]<-force_r2*lp*fp*cross[7,8]  
    
    infectmat[5,8,6,1]<-force_s*lp*fp*cross[8,1] #Reinfection with latent DS-TB from latent triple-res TB
    infectmat[5,8,6,2]<-force_r*lp*fp*cross[8,2]  
    infectmat[5,8,6,3]<-force_f*lp*fp*cross[8,3] 
    infectmat[5,8,6,4]<-force_p*lp*fp*cross[8,4]  
    infectmat[5,8,6,5]<-force_rf*lp*fp*cross[8,5] 
    infectmat[5,8,6,6]<-force_rp*lp*fp*cross[8,6] 
    infectmat[5,8,6,7]<-force_fp*lp*fp*cross[8,7]
    
    # convert infectmat into a vector
    infectv<-c(infectmat)
    
    #Rate of scale-up: linear over first 5 years, then flat
    #This will serve as a multplier for the entries in the master matrix
    diff1<-t-time1 #difference from time of Rif res introduction to current time
    diff2<-t-time2 #difference from time of FQ res introduction to current time
    diff3<-t-intro #difference from time of FQ-base regimen introduction to current time
    
    diff1<-ifelse(diff1>0, diff1, 0)
    diff2<-ifelse(diff2>0, diff2, 0)
    diff3<-ifelse(diff3>0, diff3, 0)
    
    scale1<- ifelse(diff1<5 ,diff1/5, 1) #scaling factor over 5 yrs for Rif resistance
    scale2<- ifelse(diff2<5 ,diff2/5, 1) #scaling factor over 5 yrs for FQ resistance
    scale3<- ifelse(diff3<5 ,diff3/5, 1) #scaling factor over 5 yrs for FQ-based regimen
    
    #Choose master matrix by era
      if(t<time1){
        masterv<-masterv0
      }else if(t>=time1&&t<time2){
        masterv<-masterv1*scale1+masterv0*(1-scale1)
      }else if(t>=time2&&t<intro){
        masterv<-masterv2*scale2+masterv1*(1-scale2)
      }else if(t>=intro){
        masterv<-masterv2b*scale3+masterv2*(1-scale3)
      }
      
        
    #Add infectv to masterv (add density-dependent rates to fixed rates)
    masterv3<-masterv+infectv
    
    #Convert vector to matrix
    mastermatr <-array(masterv3,dim=c(48,48))
    
    #Sum deaths at each timestep 
    dth<-array(c(1:48),dim=c(48,2))
    dthb<-c(1:48)
    deaths<-sum(mastermatr[dth]*statev[dthb]) # this computes mastermatr[x,x]*statev[x]; (mortality rate * # in compartment)
    
    births<-deaths #Set births equal to death to maintain stable population
    
    #Rates of change vector
    delta <- rep(0,48) #Start with 0 values throughout
    
    delta[1] <- sum(mastermatr[,1]*statev) -  #rates of movement INTO compartment 1
      sum(mastermatr[1,]*statev[1]) - #rates of movement OUT of compartment 1
      mastermatr[1,1]*statev[1] + #account for double-counting of mastermatr[1,1]*statev[1]
      births #Add in births (all births in susceptible state)
    
    delta[2] <- sum(mastermatr[,2]*statev) - #rates of movement INTO compartment 2
      sum(mastermatr[2,]*statev[2]) - #rates of movement OUT of compartment 2
      mastermatr[2,2]*statev[2]  #account for double-counting of mastermatr[2,2]*statev[2]
    
    delta[3] <- sum(mastermatr[,3]*statev) - 
      sum(mastermatr[3,]*statev[3]) - 
      mastermatr[3,3]*statev[3]
    
    delta[4] <- sum(mastermatr[,4]*statev) - 
      sum(mastermatr[4,]*statev[4]) - 
      mastermatr[4,4]*statev[4]
    
    delta[5] <- sum(mastermatr[,5]*statev) - 
      sum(mastermatr[5,]*statev[5]) - 
      mastermatr[5,5]*statev[5]
    
    delta[6] <- sum(mastermatr[,6]*statev) - 
      sum(mastermatr[6,]*statev[6]) - 
      mastermatr[6,6]*statev[6] 
    
    delta[7] <- sum(mastermatr[,7]*statev) - 
      sum(mastermatr[7,]*statev[7]) - 
      mastermatr[7,7]*statev[7]
    
    delta[8] <- sum(mastermatr[,8]*statev) - 
      sum(mastermatr[8,]*statev[8]) - 
      mastermatr[8,8]*statev[8]
    
    delta[9] <- sum(mastermatr[,9]*statev) - 
      sum(mastermatr[9,]*statev[9]) - 
      mastermatr[9,9]*statev[9]
    
    delta[10] <- sum(mastermatr[,10]*statev) - 
      sum(mastermatr[10,]*statev[10]) - 
      mastermatr[10,10]*statev[10]
    
    delta[11] <- sum(mastermatr[,11]*statev) - 
      sum(mastermatr[11,]*statev[11]) - 
      mastermatr[11,11]*statev[11]
    
    delta[12] <- sum(mastermatr[,12]*statev) - 
      sum(mastermatr[12,]*statev[12]) - 
      mastermatr[12,12]*statev[12]
    
    delta[13] <- sum(mastermatr[,13]*statev) - 
      sum(mastermatr[13,]*statev[13]) - 
      mastermatr[13,13]*statev[13]
    
    delta[14] <- sum(mastermatr[,14]*statev) - 
      sum(mastermatr[14,]*statev[14]) - 
      mastermatr[14,14]*statev[14]
    
    delta[15] <- sum(mastermatr[,15]*statev) - 
      sum(mastermatr[15,]*statev[15]) - 
      mastermatr[15,15]*statev[15]
    
    delta[16] <- sum(mastermatr[,16]*statev) - 
      sum(mastermatr[16,]*statev[16]) - 
      mastermatr[16,16]*statev[16]
    
    delta[17] <- sum(mastermatr[,17]*statev) - 
      sum(mastermatr[17,]*statev[17]) - 
      mastermatr[17,17]*statev[17] 
    
    delta[18] <- sum(mastermatr[,18]*statev) - 
      sum(mastermatr[18,]*statev[18]) - 
      mastermatr[18,18]*statev[18]  
    
    delta[19] <- sum(mastermatr[,19]*statev) - 
      sum(mastermatr[19,]*statev[19]) - 
      mastermatr[19,19]*statev[19]  
    
    delta[20] <- sum(mastermatr[,20]*statev) - 
      sum(mastermatr[20,]*statev[20]) - 
      mastermatr[20,20]*statev[20] 
    
    delta[21] <- sum(mastermatr[,21]*statev) - 
      sum(mastermatr[21,]*statev[21]) - 
      mastermatr[21,21]*statev[21]  
    
    delta[22] <- sum(mastermatr[,22]*statev) - 
      sum(mastermatr[22,]*statev[22]) - 
      mastermatr[22,22]*statev[22]  
    
    delta[23] <- sum(mastermatr[,23]*statev) - 
      sum(mastermatr[23,]*statev[23]) - 
      mastermatr[23,23]*statev[23]   
    
    delta[24] <- sum(mastermatr[,24]*statev) - 
      sum(mastermatr[24,]*statev[24]) - 
      mastermatr[24,24]*statev[24]
    
    delta[25] <- sum(mastermatr[,25]*statev) - 
      sum(mastermatr[25,]*statev[25]) - 
      mastermatr[25,25]*statev[25]   
    
    delta[26] <- sum(mastermatr[,26]*statev) - 
      sum(mastermatr[26,]*statev[26]) - 
      mastermatr[26,26]*statev[26]
    
    delta[27] <- sum(mastermatr[,27]*statev) - 
      sum(mastermatr[27,]*statev[27]) - 
      mastermatr[27,27]*statev[27]
    
    delta[28] <- sum(mastermatr[,28]*statev) - 
      sum(mastermatr[28,]*statev[28]) - 
      mastermatr[28,28]*statev[28]   
    
    delta[29] <- sum(mastermatr[,29]*statev) - 
      sum(mastermatr[29,]*statev[29]) - 
      mastermatr[29,29]*statev[29]   
    
    delta[30] <- sum(mastermatr[,30]*statev) - 
      sum(mastermatr[30,]*statev[30]) - 
      mastermatr[30,30]*statev[30]
    
    delta[31] <- sum(mastermatr[,31]*statev) - 
      sum(mastermatr[31,]*statev[31]) - 
      mastermatr[31,31]*statev[31]
    
    delta[32] <- sum(mastermatr[,32]*statev) - 
      sum(mastermatr[32,]*statev[32]) - 
      mastermatr[32,32]*statev[32]
    
    delta[33] <- sum(mastermatr[,33]*statev) - 
      sum(mastermatr[33,]*statev[33]) - 
      mastermatr[33,33]*statev[33]
    
    delta[34] <- sum(mastermatr[,34]*statev) - 
      sum(mastermatr[34,]*statev[34]) - 
      mastermatr[34,34]*statev[34]
    
    delta[35] <- sum(mastermatr[,35]*statev) - 
      sum(mastermatr[35,]*statev[35]) - 
      mastermatr[35,35]*statev[35]
    
    delta[36] <- sum(mastermatr[,36]*statev) - 
      sum(mastermatr[36,]*statev[36]) - 
      mastermatr[36,36]*statev[36]
    
    delta[37] <- sum(mastermatr[,37]*statev) - 
      sum(mastermatr[37,]*statev[37]) - 
      mastermatr[37,37]*statev[37]
    
    delta[38] <- sum(mastermatr[,38]*statev) - 
      sum(mastermatr[38,]*statev[38]) - 
      mastermatr[38,38]*statev[38]
    
    delta[39] <- sum(mastermatr[,39]*statev) - 
      sum(mastermatr[39,]*statev[39]) - 
      mastermatr[39,39]*statev[39]
    
    delta[40] <- sum(mastermatr[,40]*statev) - 
      sum(mastermatr[40,]*statev[40]) - 
      mastermatr[40,40]*statev[40]
    
    delta[41] <- sum(mastermatr[,41]*statev) - 
      sum(mastermatr[41,]*statev[41]) - 
      mastermatr[41,41]*statev[41]
    
    delta[42] <- sum(mastermatr[,42]*statev) - 
      sum(mastermatr[42,]*statev[42]) - 
      mastermatr[42,42]*statev[42]
    
    delta[43] <- sum(mastermatr[,43]*statev) - 
      sum(mastermatr[43,]*statev[43]) - 
      mastermatr[43,43]*statev[43]
    
    delta[44] <- sum(mastermatr[,44]*statev) - 
      sum(mastermatr[44,]*statev[44]) - 
      mastermatr[44,44]*statev[44]
    
    delta[45] <- sum(mastermatr[,45]*statev) - 
      sum(mastermatr[45,]*statev[45]) - 
      mastermatr[45,45]*statev[45]
    
    delta[46] <- sum(mastermatr[,46]*statev) - 
      sum(mastermatr[46,]*statev[46]) - 
      mastermatr[46,46]*statev[46]
    
    delta[47] <- sum(mastermatr[,47]*statev) - 
      sum(mastermatr[47,]*statev[47]) - 
      mastermatr[47,47]*statev[47]
    
    delta[48] <- sum(mastermatr[,48]*statev) - 
      sum(mastermatr[48,]*statev[48]) - 
      mastermatr[48,48]*statev[48]
    
    # return the rates of change
    list(c(delta))
  })
}





#############################################################
#FUNCTION III: COUNTER FUNCTION TO COMPILE SIMULATION OUTPUT#
#############################################################

#Outputs: (1) Overall incidence, (2) Incidence by resistance state and (3) New vs. retreatment
#Notified cases: count relapses and failures; count new and retreatment separately

Counter<-function(output,parameters,step){ 
  #function arguments are (1) output from ODE solver, (2) state/transition vector, (3) step size
  
  x<-as.list(c(output,parameters,step))
  with(x,{
    
    lstate<-length(output) #length of each row in "output"  = (number of state variables) + 1 
    lmaster<-length(parameters)-2 #length of transition rate vector

    time<-unlist(x[1]) #timestep
    statev2<-unlist(x[2:lstate]) # state vector
    
    masterv30<-unlist(x[(lstate+1):(lstate+lmaster/4)])  # rate of change vectors, by era
    masterv31<-unlist(x[(lstate+lmaster/4+1):(lstate+lmaster*2/4)])
    masterv32<-unlist(x[(lstate+lmaster*2/4+1):(lstate+lmaster*3/4)])
    masterv32b<-unlist(x[(lstate+lmaster*3/4+1):(lstate+lmaster)])
    
    time1<-unlist(x[(lstate+lmaster+1)]) #time of intro of Rif res
    time2<-unlist(x[(lstate+lmaster+2)])  #time of intro of FQ res 
    step<-unlist(x[(lstate+lmaster+3)]) #step size
    
    statematr<-array(statev2,dim=c(6,8)) #convert state vector to matrix
    
    #Need to replicate rates of change 
    #Next section of code is same as TBDx function
   
    
###########################################################################    
    # force of infection, density-dependent
    force_s <- tr*(statev2[3]+statev2[6]+relinf*statev2[4])
    force_r <- tr*(statev2[9]+statev2[12]+relinf*statev2[10])*fit_r
    force_f <- tr*(statev2[15]+statev2[18]+relinf*statev2[16])*fit_f
    force_p <- tr*(statev2[21]+statev2[24]+relinf*statev2[22])*fit_p
    force_rf <- tr*(statev2[27]+statev2[30]+relinf*statev2[28])*fit_rf
    force_rp <- tr*(statev2[33]+statev2[36]+relinf*statev2[34])*fit_rp
    force_fp <- tr*(statev2[39]+statev2[42]+relinf*statev2[40])*fit_fp
    force_r2 <- tr*(statev2[45]+statev2[48]+relinf*statev2[46])*fit_r2
    
    # now, create an "infection matrix"
    infectmat <- array(0, dim=c(6,8,6,8))
    
    #S->L1
    infectmat[1,1,2,1]<-force_s*(1-fp) 
    infectmat[1,1,2,2]<-force_r*(1-fp) 
    infectmat[1,1,2,3]<-force_f*(1-fp) 
    infectmat[1,1,2,4]<-force_p*(1-fp)
    infectmat[1,1,2,5]<-force_rf*(1-fp) 
    infectmat[1,1,2,6]<-force_rp*(1-fp) 
    infectmat[1,1,2,7]<-force_fp*(1-fp)
    infectmat[1,1,2,8]<-force_r2*(1-fp)
    
    #S->A1
    infectmat[1,1,3,1]<-force_s*fp 
    infectmat[1,1,3,2]<-force_r*fp 
    infectmat[1,1,3,3]<-force_f*fp 
    infectmat[1,1,3,4]<-force_p*fp
    infectmat[1,1,3,5]<-force_rf*fp 
    infectmat[1,1,3,6]<-force_rp*fp
    infectmat[1,1,3,7]<-force_fp*fp
    infectmat[1,1,3,8]<-force_r2*fp 
    
    #L1->L2
    infectmat[2,1,2,2]<-force_r*lp*(1-fp)*cross[1,2] 
    infectmat[2,1,2,3]<-force_f*lp*(1-fp)*cross[1,3] 
    infectmat[2,1,2,4]<-force_p*lp*(1-fp)*cross[1,4]
    infectmat[2,1,2,5]<-force_rf*lp*(1-fp)*cross[1,5] 
    infectmat[2,1,2,6]<-force_rp*lp*(1-fp)*cross[1,6]
    infectmat[2,1,2,7]<-force_fp*lp*(1-fp)*cross[1,7]
    infectmat[2,1,2,8]<-force_r2*lp*(1-fp)*cross[1,8] 
    
    infectmat[2,2,2,1]<-force_s*lp*(1-fp)*cross[2,1] 
    infectmat[2,2,2,3]<-force_f*lp*(1-fp)*cross[2,3] 
    infectmat[2,2,2,4]<-force_p*lp*(1-fp)*cross[2,4]
    infectmat[2,2,2,5]<-force_rf*lp*(1-fp)*cross[2,5] 
    infectmat[2,2,2,6]<-force_rp*lp*(1-fp)*cross[2,6]
    infectmat[2,2,2,7]<-force_fp*lp*(1-fp)*cross[2,7]
    infectmat[2,2,2,8]<-force_r2*lp*(1-fp)*cross[2,8] 
    
    infectmat[2,3,2,1]<-force_s*lp*(1-fp)*cross[3,1] 
    infectmat[2,3,2,2]<-force_r*lp*(1-fp)*cross[3,2] 
    infectmat[2,3,2,4]<-force_p*lp*(1-fp)*cross[3,4]
    infectmat[2,3,2,5]<-force_rf*lp*(1-fp)*cross[3,5] 
    infectmat[2,3,2,6]<-force_rp*lp*(1-fp)*cross[3,6]
    infectmat[2,3,2,7]<-force_fp*lp*(1-fp)*cross[3,7]
    infectmat[2,3,2,8]<-force_r2*lp*(1-fp)*cross[3,8] 
    
    infectmat[2,4,2,1]<-force_s*lp*(1-fp)*cross[4,1] 
    infectmat[2,4,2,2]<-force_r*lp*(1-fp)*cross[4,2] 
    infectmat[2,4,2,3]<-force_f*lp*(1-fp)*cross[4,3]
    infectmat[2,4,2,5]<-force_rf*lp*(1-fp)*cross[4,5] 
    infectmat[2,4,2,6]<-force_rp*lp*(1-fp)*cross[4,6]
    infectmat[2,4,2,7]<-force_fp*lp*(1-fp)*cross[4,7]
    infectmat[2,4,2,8]<-force_r2*lp*(1-fp)*cross[4,8]
    
    infectmat[2,5,2,1]<-force_s*lp*(1-fp)*cross[5,1] 
    infectmat[2,5,2,2]<-force_r*lp*(1-fp)*cross[5,2] 
    infectmat[2,5,2,3]<-force_f*lp*(1-fp)*cross[5,3]
    infectmat[2,5,2,4]<-force_p*lp*(1-fp)*cross[5,4] 
    infectmat[2,5,2,6]<-force_rp*lp*(1-fp)*cross[5,6]
    infectmat[2,5,2,7]<-force_fp*lp*(1-fp)*cross[5,7]
    infectmat[2,5,2,8]<-force_r2*lp*(1-fp)*cross[5,8] 
    
    infectmat[2,6,2,1]<-force_s*lp*(1-fp)*cross[6,1] 
    infectmat[2,6,2,2]<-force_r*lp*(1-fp)*cross[6,2] 
    infectmat[2,6,2,3]<-force_f*lp*(1-fp)*cross[6,3]
    infectmat[2,6,2,4]<-force_p*lp*(1-fp)*cross[6,4] 
    infectmat[2,6,2,5]<-force_rf*lp*(1-fp)*cross[6,5]
    infectmat[2,6,2,7]<-force_fp*lp*(1-fp)*cross[6,7]
    infectmat[2,6,2,8]<-force_r2*lp*(1-fp)*cross[6,8] 
    
    infectmat[2,7,2,1]<-force_s*lp*(1-fp)*cross[7,1] 
    infectmat[2,7,2,2]<-force_r*lp*(1-fp)*cross[7,2]  
    infectmat[2,7,2,3]<-force_f*lp*(1-fp)*cross[7,3] 
    infectmat[2,7,2,4]<-force_p*lp*(1-fp)*cross[7,4]  
    infectmat[2,7,2,5]<-force_rf*lp*(1-fp)*cross[7,5] 
    infectmat[2,7,2,6]<-force_rp*lp*(1-fp)*cross[7,6] 
    infectmat[2,7,2,8]<-force_r2*lp*(1-fp)*cross[7,8]  
    
    infectmat[2,8,2,1]<-force_s*lp*(1-fp)*cross[8,1] 
    infectmat[2,8,2,2]<-force_r*lp*(1-fp)*cross[8,2]  
    infectmat[2,8,2,3]<-force_f*lp*(1-fp)*cross[8,3] 
    infectmat[2,8,2,4]<-force_p*lp*(1-fp)*cross[8,4]  
    infectmat[2,8,2,5]<-force_rf*lp*(1-fp)*cross[8,5] 
    infectmat[2,8,2,6]<-force_rp*lp*(1-fp)*cross[8,6] 
    infectmat[2,8,2,7]<-force_fp*lp*(1-fp)*cross[8,7]  
    
    #L1->A1
    infectmat[2,1,3,2]<-force_r*lp*fp*cross[1,2] 
    infectmat[2,1,3,3]<-force_f*lp*fp*cross[1,3] 
    infectmat[2,1,3,4]<-force_p*lp*fp*cross[1,4]
    infectmat[2,1,3,5]<-force_rf*lp*fp*cross[1,5] 
    infectmat[2,1,3,6]<-force_rp*lp*fp*cross[1,6]
    infectmat[2,1,3,7]<-force_fp*lp*fp*cross[1,7]
    infectmat[2,1,3,8]<-force_r2*lp*fp*cross[1,8] 
    
    infectmat[2,2,3,1]<-force_s*lp*fp*cross[2,1] 
    infectmat[2,2,3,3]<-force_f*lp*fp*cross[2,3] 
    infectmat[2,2,3,4]<-force_p*lp*fp*cross[2,4]
    infectmat[2,2,3,5]<-force_rf*lp*fp*cross[2,5] 
    infectmat[2,2,3,6]<-force_rp*lp*fp*cross[2,6]
    infectmat[2,2,3,7]<-force_fp*lp*fp*cross[2,7]
    infectmat[2,2,3,8]<-force_r2*lp*fp*cross[2,8] 
    
    infectmat[2,3,3,1]<-force_s*lp*fp*cross[3,1] 
    infectmat[2,3,3,2]<-force_r*lp*fp*cross[3,2] 
    infectmat[2,3,3,4]<-force_p*lp*fp*cross[3,4]
    infectmat[2,3,3,5]<-force_rf*lp*fp*cross[3,5] 
    infectmat[2,3,3,6]<-force_rp*lp*fp*cross[3,6]
    infectmat[2,3,3,7]<-force_fp*lp*fp*cross[3,7]
    infectmat[2,3,3,8]<-force_r2*lp*fp*cross[3,8] 
    
    infectmat[2,4,3,1]<-force_s*lp*fp*cross[4,1] 
    infectmat[2,4,3,2]<-force_r*lp*fp*cross[4,2] 
    infectmat[2,4,3,3]<-force_f*lp*fp*cross[4,3]
    infectmat[2,4,3,5]<-force_rf*lp*fp*cross[4,5] 
    infectmat[2,4,3,6]<-force_rp*lp*fp*cross[4,6]
    infectmat[2,4,3,7]<-force_fp*lp*fp*cross[4,7]
    infectmat[2,4,3,8]<-force_r2*lp*fp*cross[4,8]
    
    infectmat[2,5,3,1]<-force_s*lp*fp*cross[5,1] 
    infectmat[2,5,3,2]<-force_r*lp*fp*cross[5,2] 
    infectmat[2,5,3,3]<-force_f*lp*fp*cross[5,3]
    infectmat[2,5,3,4]<-force_p*lp*fp*cross[5,4] 
    infectmat[2,5,3,6]<-force_rp*lp*fp*cross[5,6]
    infectmat[2,5,3,7]<-force_fp*lp*fp*cross[5,7]
    infectmat[2,5,3,8]<-force_r2*lp*fp*cross[5,8] 
    
    infectmat[2,6,3,1]<-force_s*lp*fp*cross[6,1] 
    infectmat[2,6,3,2]<-force_r*lp*fp*cross[6,2] 
    infectmat[2,6,3,3]<-force_f*lp*fp*cross[6,3]
    infectmat[2,6,3,4]<-force_p*lp*fp*cross[6,4] 
    infectmat[2,6,3,5]<-force_rf*lp*fp*cross[6,5]
    infectmat[2,6,3,7]<-force_fp*lp*fp*cross[6,7]
    infectmat[2,6,3,8]<-force_r2*lp*fp*cross[6,8] 
    
    infectmat[2,7,3,1]<-force_s*lp*fp*cross[7,1] 
    infectmat[2,7,3,2]<-force_r*lp*fp*cross[7,2]  
    infectmat[2,7,3,3]<-force_f*lp*fp*cross[7,3] 
    infectmat[2,7,3,4]<-force_p*lp*fp*cross[7,4]  
    infectmat[2,7,3,5]<-force_rf*lp*fp*cross[7,5] 
    infectmat[2,7,3,6]<-force_rp*lp*fp*cross[7,6] 
    infectmat[2,7,3,8]<-force_r2*lp*fp*cross[7,8]  
    
    infectmat[2,8,3,1]<-force_s*lp*fp*cross[8,1] 
    infectmat[2,8,3,2]<-force_r*lp*fp*cross[8,2]  
    infectmat[2,8,3,3]<-force_f*lp*fp*cross[8,3] 
    infectmat[2,8,3,4]<-force_p*lp*fp*cross[8,4]  
    infectmat[2,8,3,5]<-force_rf*lp*fp*cross[8,5] 
    infectmat[2,8,3,6]<-force_rp*lp*fp*cross[8,6] 
    infectmat[2,8,3,7]<-force_fp*lp*fp*cross[8,7]  
    
    #L2->L2
    infectmat[5,1,5,2]<-force_r*lp*(1-fp)*cross[1,2] 
    infectmat[5,1,5,3]<-force_f*lp*(1-fp)*cross[1,3] 
    infectmat[5,1,5,4]<-force_p*lp*(1-fp)*cross[1,4]
    infectmat[5,1,5,5]<-force_rf*lp*(1-fp)*cross[1,5] 
    infectmat[5,1,5,6]<-force_rp*lp*(1-fp)*cross[1,6]
    infectmat[5,1,5,7]<-force_fp*lp*(1-fp)*cross[1,7]
    infectmat[5,1,5,8]<-force_r2*lp*(1-fp)*cross[1,8] 
    
    infectmat[5,2,5,1]<-force_s*lp*(1-fp)*cross[2,1] 
    infectmat[5,2,5,3]<-force_f*lp*(1-fp)*cross[2,3] 
    infectmat[5,2,5,4]<-force_p*lp*(1-fp)*cross[2,4]
    infectmat[5,2,5,5]<-force_rf*lp*(1-fp)*cross[2,5] 
    infectmat[5,2,5,6]<-force_rp*lp*(1-fp)*cross[2,6]
    infectmat[5,2,5,7]<-force_fp*lp*(1-fp)*cross[2,7]
    infectmat[5,2,5,8]<-force_r2*lp*(1-fp)*cross[2,8] 
    
    infectmat[5,3,5,1]<-force_s*lp*(1-fp)*cross[3,1] 
    infectmat[5,3,5,2]<-force_r*lp*(1-fp)*cross[3,2] 
    infectmat[5,3,5,4]<-force_p*lp*(1-fp)*cross[3,4]
    infectmat[5,3,5,5]<-force_rf*lp*(1-fp)*cross[3,5] 
    infectmat[5,3,5,6]<-force_rp*lp*(1-fp)*cross[3,6]
    infectmat[5,3,5,7]<-force_fp*lp*(1-fp)*cross[3,7]
    infectmat[5,3,5,8]<-force_r2*lp*(1-fp)*cross[3,8] 
    
    infectmat[5,4,5,1]<-force_s*lp*(1-fp)*cross[4,1] 
    infectmat[5,4,5,2]<-force_r*lp*(1-fp)*cross[4,2] 
    infectmat[5,4,5,3]<-force_f*lp*(1-fp)*cross[4,3]
    infectmat[5,4,5,5]<-force_rf*lp*(1-fp)*cross[4,5] 
    infectmat[5,4,5,6]<-force_rp*lp*(1-fp)*cross[4,6]
    infectmat[5,4,5,7]<-force_fp*lp*(1-fp)*cross[4,7]
    infectmat[5,4,5,8]<-force_r2*lp*(1-fp)*cross[4,8]
    
    infectmat[5,5,5,1]<-force_s*lp*(1-fp)*cross[5,1] 
    infectmat[5,5,5,2]<-force_r*lp*(1-fp)*cross[5,2] 
    infectmat[5,5,5,3]<-force_f*lp*(1-fp)*cross[5,3]
    infectmat[5,5,5,4]<-force_p*lp*(1-fp)*cross[5,4] 
    infectmat[5,5,5,6]<-force_rp*lp*(1-fp)*cross[5,6]
    infectmat[5,5,5,7]<-force_fp*lp*(1-fp)*cross[5,7]
    infectmat[5,5,5,8]<-force_r2*lp*(1-fp)*cross[5,8] 
    
    infectmat[5,6,5,1]<-force_s*lp*(1-fp)*cross[6,1] 
    infectmat[5,6,5,2]<-force_r*lp*(1-fp)*cross[6,2] 
    infectmat[5,6,5,3]<-force_f*lp*(1-fp)*cross[6,3]
    infectmat[5,6,5,4]<-force_p*lp*(1-fp)*cross[6,4] 
    infectmat[5,6,5,5]<-force_rf*lp*(1-fp)*cross[6,5]
    infectmat[5,6,5,7]<-force_fp*lp*(1-fp)*cross[6,7]
    infectmat[5,6,5,8]<-force_r2*lp*(1-fp)*cross[6,8] 
    
    infectmat[5,7,5,1]<-force_s*lp*(1-fp)*cross[7,1] 
    infectmat[5,7,5,2]<-force_r*lp*(1-fp)*cross[7,2]  
    infectmat[5,7,5,3]<-force_f*lp*(1-fp)*cross[7,3] 
    infectmat[5,7,5,4]<-force_p*lp*(1-fp)*cross[7,4]  
    infectmat[5,7,5,5]<-force_rf*lp*(1-fp)*cross[7,5] 
    infectmat[5,7,5,6]<-force_rp*lp*(1-fp)*cross[7,6] 
    infectmat[5,7,5,8]<-force_r2*lp*(1-fp)*cross[7,8]  
    
    infectmat[5,8,5,1]<-force_s*lp*(1-fp)*cross[8,1] 
    infectmat[5,8,5,2]<-force_r*lp*(1-fp)*cross[8,2]  
    infectmat[5,8,5,3]<-force_f*lp*(1-fp)*cross[8,3] 
    infectmat[5,8,5,4]<-force_p*lp*(1-fp)*cross[8,4]  
    infectmat[5,8,5,5]<-force_rf*lp*(1-fp)*cross[8,5] 
    infectmat[5,8,5,6]<-force_rp*lp*(1-fp)*cross[8,6] 
    infectmat[5,8,5,7]<-force_fp*lp*(1-fp)*cross[8,7]  
    
    #L2->A2
    infectmat[5,1,6,2]<-force_r*lp*fp*cross[1,2] 
    infectmat[5,1,6,3]<-force_f*lp*fp*cross[1,3] 
    infectmat[5,1,6,4]<-force_p*lp*fp*cross[1,4]
    infectmat[5,1,6,5]<-force_rf*lp*fp*cross[1,5] 
    infectmat[5,1,6,6]<-force_rp*lp*fp*cross[1,6]
    infectmat[5,1,6,7]<-force_fp*lp*fp*cross[1,7]
    infectmat[5,1,6,8]<-force_r2*lp*fp*cross[1,8] 
    
    infectmat[5,2,6,1]<-force_s*lp*fp*cross[2,1] 
    infectmat[5,2,6,3]<-force_f*lp*fp*cross[2,3] 
    infectmat[5,2,6,4]<-force_p*lp*fp*cross[2,4]
    infectmat[5,2,6,5]<-force_rf*lp*fp*cross[2,5] 
    infectmat[5,2,6,6]<-force_rp*lp*fp*cross[2,6]
    infectmat[5,2,6,7]<-force_fp*lp*fp*cross[2,7]
    infectmat[5,2,6,8]<-force_r2*lp*fp*cross[2,8] 
    
    infectmat[5,3,6,1]<-force_s*lp*fp*cross[3,1] 
    infectmat[5,3,6,2]<-force_r*lp*fp*cross[3,2] 
    infectmat[5,3,6,4]<-force_p*lp*fp*cross[3,4]
    infectmat[5,3,6,5]<-force_rf*lp*fp*cross[3,5] 
    infectmat[5,3,6,6]<-force_rp*lp*fp*cross[3,6]
    infectmat[5,3,6,7]<-force_fp*lp*fp*cross[3,7]
    infectmat[5,3,6,8]<-force_r2*lp*fp*cross[3,8] 
    
    infectmat[5,4,6,1]<-force_s*lp*fp*cross[4,1] 
    infectmat[5,4,6,2]<-force_r*lp*fp*cross[4,2] 
    infectmat[5,4,6,3]<-force_f*lp*fp*cross[4,3]
    infectmat[5,4,6,5]<-force_rf*lp*fp*cross[4,5] 
    infectmat[5,4,6,6]<-force_rp*lp*fp*cross[4,6]
    infectmat[5,4,6,7]<-force_fp*lp*fp*cross[4,7]
    infectmat[5,4,6,8]<-force_r2*lp*fp*cross[4,8]
    
    infectmat[5,5,6,1]<-force_s*lp*fp*cross[5,1] 
    infectmat[5,5,6,2]<-force_r*lp*fp*cross[5,2] 
    infectmat[5,5,6,3]<-force_f*lp*fp*cross[5,3]
    infectmat[5,5,6,4]<-force_p*lp*fp*cross[5,4] 
    infectmat[5,5,6,6]<-force_rp*lp*fp*cross[5,6]
    infectmat[5,5,6,7]<-force_fp*lp*fp*cross[5,7]
    infectmat[5,5,6,8]<-force_r2*lp*fp*cross[5,8] 
    
    infectmat[5,6,6,1]<-force_s*lp*fp*cross[6,1] 
    infectmat[5,6,6,2]<-force_r*lp*fp*cross[6,2] 
    infectmat[5,6,6,3]<-force_f*lp*fp*cross[6,3]
    infectmat[5,6,6,4]<-force_p*lp*fp*cross[6,4] 
    infectmat[5,6,6,5]<-force_rf*lp*fp*cross[6,5]
    infectmat[5,6,6,7]<-force_fp*lp*fp*cross[6,7]
    infectmat[5,6,6,8]<-force_r2*lp*fp*cross[6,8] 
    
    infectmat[5,7,6,1]<-force_s*lp*fp*cross[7,1] 
    infectmat[5,7,6,2]<-force_r*lp*fp*cross[7,2]  
    infectmat[5,7,6,3]<-force_f*lp*fp*cross[7,3] 
    infectmat[5,7,6,4]<-force_p*lp*fp*cross[7,4]  
    infectmat[5,7,6,5]<-force_rf*lp*fp*cross[7,5] 
    infectmat[5,7,6,6]<-force_rp*lp*fp*cross[7,6] 
    infectmat[5,7,6,8]<-force_r2*lp*fp*cross[7,8]  
    
    infectmat[5,8,6,1]<-force_s*lp*fp*cross[8,1] 
    infectmat[5,8,6,2]<-force_r*lp*fp*cross[8,2]  
    infectmat[5,8,6,3]<-force_f*lp*fp*cross[8,3] 
    infectmat[5,8,6,4]<-force_p*lp*fp*cross[8,4]  
    infectmat[5,8,6,5]<-force_rf*lp*fp*cross[8,5] 
    infectmat[5,8,6,6]<-force_rp*lp*fp*cross[8,6] 
    infectmat[5,8,6,7]<-force_fp*lp*fp*cross[8,7] 
    
    # convert infectmat into a vector
    infectv2<-c(infectmat)
    
    #Rate of scale-up: linear over first 5 years, then flat
    #This will serve as a multplier for the entries in the master matrix
    diff1<-time-time1 #difference from time of Rif res introduction to current time
    diff2<-time-time2 #difference from time of FQ res introduction to current time
    diff3<-time-intro #difference from time of FQ-base regimen introduction to current time
    
    diff1<-ifelse(diff1>0, diff1, 0)
    diff2<-ifelse(diff2>0, diff2, 0)
    diff3<-ifelse(diff3>0, diff3, 0)
    
    scale1<- ifelse(diff1<5 ,diff1/5, 1) #scaling factor over 5 yrs for Rif resistance
    scale2<- ifelse(diff2<5 ,diff2/5, 1) #scaling factor over 5 yrs for FQ resistance
    scale3<- ifelse(diff3<5 ,diff3/5, 1) #scaling factor over 5 yrs for FQ-based regimen
    
    #Choose master matrix by era
    if(time<time1){
      masterv3<-masterv30
    }else if(time>=time1&&time<time2){
      masterv3<-masterv31*scale1+masterv30*(1-scale1)
    }else if(time>=time2&&time<intro){
      masterv3<-masterv32*scale2+masterv31*(1-scale2)
    }else if(time>=intro){
      masterv3<-masterv32b*scale3+masterv32*(1-scale3)
    }
        
    
    # add infectv to masterv
    masterv4<-masterv3+infectv2
    
    # Recreate the change matrix
    mastermatr2 <-array(masterv4,dim=c(6,8,6,8)) #Full change matrix with infections
########################################################################################### 

    i <- rep(1:6,each=8) #Rows; This creates 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4, etc...
    j <- rep(1:8,6) #Columns; This creates 1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6, etc...
    k <- array(c(i,j),dim=c(48,2)) #This creates combo of row/column for every possible state
    m <- array(c(k,k),dim=c(48,4)) #This creates combo for every  transition of state x-->state x (for mortality)

        
    #Sum deaths
    n <- rep(c(3,4,6),each=8) # Rows 
    o <- rep(1:8,3) # Columns
    p <- array(c(n,o),dim=c(24,2)) #This creates combo of row;column for every possible state, going down the columns of population matrix
    q <- array(c(p,p),dim=c(24,4)) #This creates combo for transition of state x-->state x
    
    #Compute # of total and TB-specific deaths
    Deaths<-sum(mastermatr2[m]*statematr[k]) #This is # of deaths that would occur over 1 yr given rate of change
    Deaths<-Deaths*step #multiply by duration of timestep to obtain actual # of total deaths during timestep
    TBmort<-mastermatr2[q] #Create vector of mortality rates for active TB/Rx states
    TBmort<-TBmort-m1 #subtract background mortality
    TBdeaths<-sum(TBmort*statematr[p]) #TB-specific deaths that would occur over 1 yr
    TBdeaths<-TBdeaths*step #multiply by timestep to get actual # of deaths in timestep
    
    #Proportion of defaulters in the active TB rx-experiences compartment
    propdef<-sum(mastermatr2[3,,6,],mastermatr2[4,,6,])/  
      sum(mastermatr2[3,,6,],mastermatr2[4,,6,],mastermatr2[5,,6,]) 
    
    #Incidence by treatment and drug resistance state
    n1<-rep(c(3,3,6),each=8) 
    o<-array(rep(1:8,each=24),dim=c(24,8))
    
    p0<-array(c(rep(c(1,2,5),each=8),rep(1:8,3)),dim=c(24,2))
    
    p1<-array(c(n1,o[,1]),dim=c(24,2))
    p2<-array(c(n1,o[,2]),dim=c(24,2))
    p3<-array(c(n1,o[,3]),dim=c(24,2))
    p4<-array(c(n1,o[,4]),dim=c(24,2))
    p5<-array(c(n1,o[,5]),dim=c(24,2))
    p6<-array(c(n1,o[,6]),dim=c(24,2))
    p7<-array(c(n1,o[,7]),dim=c(24,2))
    p8<-array(c(n1,o[,8]),dim=c(24,2))
    
    q1<-array(c(p0,p1),dim=c(24,4))
    q2<-array(c(p0,p2),dim=c(24,4))
    q3<-array(c(p0,p3),dim=c(24,4))
    q4<-array(c(p0,p4),dim=c(24,4))
    q5<-array(c(p0,p5),dim=c(24,4))
    q6<-array(c(p0,p6),dim=c(24,4))
    q7<-array(c(p0,p7),dim=c(24,4))
    q8<-array(c(p0,p8),dim=c(24,4))
    
    #Incidence of new cases
    inc1<-sum(mastermatr2[q1][1:16]*statematr[p0][1:16])
    inc2<-sum(mastermatr2[q2][1:16]*statematr[p0][1:16])
    inc3<-sum(mastermatr2[q3][1:16]*statematr[p0][1:16])
    inc4<-sum(mastermatr2[q4][1:16]*statematr[p0][1:16])
    inc5<-sum(mastermatr2[q5][1:16]*statematr[p0][1:16])
    inc6<-sum(mastermatr2[q6][1:16]*statematr[p0][1:16])
    inc7<-sum(mastermatr2[q7][1:16]*statematr[p0][1:16])
    inc8<-sum(mastermatr2[q8][1:16]*statematr[p0][1:16])
    
    #Incidence of retreatment cases
    reinc1<-sum(mastermatr2[q1][17:24]*statematr[p0][17:24])
    reinc2<-sum(mastermatr2[q2][17:24]*statematr[p0][17:24])
    reinc3<-sum(mastermatr2[q3][17:24]*statematr[p0][17:24])
    reinc4<-sum(mastermatr2[q4][17:24]*statematr[p0][17:24])
    reinc5<-sum(mastermatr2[q5][17:24]*statematr[p0][17:24])
    reinc6<-sum(mastermatr2[q6][17:24]*statematr[p0][17:24])
    reinc7<-sum(mastermatr2[q7][17:24]*statematr[p0][17:24])
    reinc8<-sum(mastermatr2[q8][17:24]*statematr[p0][17:24])
    
    #Retreatment case notification (add A2->F, A2->L, F->L, F->A2)
    #Need to include not just new infections but cases of retreatment after failure/relapse
    #So rather than count infections, count treatment outcomes as these will be ~equal to # dx'd, notified
    n21<-rep(c(5,6,4,5),each=8)
    o2<-array(rep(1:8,each=32),dim=c(32,8))
    
    p20<-array(c(rep(c(4,6),each=16),rep(1:8,4)),dim=c(32,2))
    
    p21<-array(c(n21,o2[,1]),dim=c(32,2))
    p22<-array(c(n21,o2[,2]),dim=c(32,2))
    p23<-array(c(n21,o2[,3]),dim=c(32,2))
    p24<-array(c(n21,o2[,4]),dim=c(32,2))
    p25<-array(c(n21,o2[,5]),dim=c(32,2))
    p26<-array(c(n21,o2[,6]),dim=c(32,2))
    p27<-array(c(n21,o2[,7]),dim=c(32,2))
    p28<-array(c(n21,o2[,8]),dim=c(32,2))
    
    q21<-array(c(p20,p21),dim=c(32,4))
    q22<-array(c(p20,p22),dim=c(32,4))
    q23<-array(c(p20,p23),dim=c(32,4))
    q24<-array(c(p20,p24),dim=c(32,4))
    q25<-array(c(p20,p25),dim=c(32,4))
    q26<-array(c(p20,p26),dim=c(32,4))
    q27<-array(c(p20,p27),dim=c(32,4))
    q28<-array(c(p20,p28),dim=c(32,4))
    
    #Notified cases
    not1<-sum(mastermatr2[q21]*statematr[p20])
    not2<-sum(mastermatr2[q22]*statematr[p20])
    not3<-sum(mastermatr2[q23]*statematr[p20])
    not4<-sum(mastermatr2[q24]*statematr[p20])
    not5<-sum(mastermatr2[q25]*statematr[p20])
    not6<-sum(mastermatr2[q26]*statematr[p20])
    not7<-sum(mastermatr2[q27]*statematr[p20])
    not8<-sum(mastermatr2[q28]*statematr[p20])
    
    #Output of the function: new cases, retreatment cases, & notified cases by strain, + % defaulters
    Count<-c(inc1,inc2,inc3,inc4,inc5,inc6,inc7,inc8,
             reinc1,reinc2,reinc3,reinc4,reinc5,reinc6,reinc7,reinc8,
             not1,not2,not3,not4,not5,not6,not7,not8,
             propdef) 
  })
}


