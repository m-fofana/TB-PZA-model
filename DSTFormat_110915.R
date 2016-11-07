#Mariam Fofana
#November 2015
#Modeling analysis of drivers of TB drug resistance trajectories
#Data formatting code

###########################################################
#I: create a random subsample of all projected trajectories
###########################################################

#Change working directory as necessary
setwd("~/DST2")

#Enter total number of sets simulated; change as necessary
totset<-20

#Read in data files for selected trajectories, "totset" sets 
#Randomly sample # of trajectories to be sampled from each set
tmp<-sample((5000*totset),5000,replace=F)
pick<-vector(length=totset)
for(set in 1:totset){
  tmp2<-which(tmp<=5000)
  pick[set]<-length(tmp2)
  tmp<-tmp[-tmp2]
  tmp<-tmp-5000
}

#pick is now a vector of # of trajectories to sample from each of "totset" sets
#Check that sums up to 5000
head(pick)
sum(pick)

#Remove temporary variables
remove(tmp,tmp2)

#Compile input variables for all and selected trajectories
#Compile trajectory data for subsample of all trajectories 
#Compile trajectory data for full set of selected trajectories

for(set in 1:totset){
  #Read in input values 
  #tempa<-read.table(paste("Baseline/Init/dstInit110915_Baseline",set,".csv",sep=""),
  #                  sep=",",header=T)
  tempa<-read.table(paste("NoProtection/Init/dstInit110915_NoProtect",set,".csv",sep=""),
                    sep=",",header=T)

  #Identify selected simulations (must have score of 12 on selection criteria)
  keepa<-which(rowSums(tempa[,63:74])==12 &
                 tempa[,59]<tempa[,58] &
                 tempa[,60]<tempa[,58])
    
  
  #Read in trajectory data
  #tempb<-read.table(paste("Baseline/Data/dstData110915_Baseline",set,".csv",sep=""),
  #                  sep=",",header=T)
  tempb<-read.table(paste("NoProtection/Data/dstData110915_NoProtect",set,".csv",sep=""),
                    sep=",",header=T)
  
  #Extract input values
  assign(paste("init",set,sep=""),tempa) #All trajectories
  assign(paste("initsel",set,sep=""),tempa[keepa,]) #Selected trajectories
  
  #Extract trajectory data
  sub<-sample(5000,pick[set],replace=F)
  assign(paste("data",set,sep=""),tempb[sub,-1]) #Random subsample of all trajectories
  assign(paste("datasel",set,sep=""),tempb[keepa,-1]) #Selected trajectories
  
  remove(tempa,tempb,keepa)
  print(set) #Counter to monitor progress; this can take a while...
}

##########################################
#II: Combine data from all files together
##########################################

#Input variables
form<-paste("init1")
for(tmp in 2:totset){
  form<-paste(form,",","init",tmp, sep="")
}

initALL<-eval(parse(text=paste("rbind(",form,")",sep="")))
dim(initALL) #should be 100K by 74

eval(parse(text=paste("remove(",form,")",sep=""))) #remove init arrays


form<-paste("initsel1")
for(tmp in 2:totset){
  form<-paste(form,",","initsel",tmp, sep="")
}

init<-eval(parse(text=paste("rbind(",form,")",sep="")))
dim(init)

eval(parse(text=paste("remove(",form,")",sep=""))) #remove initsel arrays

#Remove input variables that are not used in simulations
#Not all of the variables are sampeld actually matter depending on the scenario that is simulated
#In this analysis, variables related to treatment regimens other than HRZE and STR do not matter
#Need to remove these from regression and correlation analyses
#Remove: variables 30,35-49
#Also need to remove variables that are set equal to each other
#Remove: variables 24, 25, 28; final # of variables is 31
## init<-init[,-c(24,25,28,30,35:49)]
## initALL<-initALL[,-c(24,25,28,30,35:49)]

#Trajectories
form<-paste("data1")
for(tmp in 2:totset){
  form<-paste(form,",","data",tmp, sep="")
}

dataALL<-eval(parse(text=paste("rbind(",form,")",sep="")))
dim(dataALL) #should be 5000 by 6292

eval(parse(text=paste("remove(",form,")",sep=""))) #remove data arrays


form<-paste("datasel1")
for(tmp in 2:totset){
  form<-paste(form,",","datasel",tmp, sep="")
}

data<-eval(parse(text=paste("rbind(",form,")",sep="")))
dim(data)

eval(parse(text=paste("remove(",form,")",sep=""))) #remove datasel arrays
remove(form)

##########################################
#III: Write to csv for future use
##########################################

#write.csv(initALL,"Baseline/Init/InputAllBsl_110915.csv")
#write.csv(init,"Baseline/Init/InputSelBsl_110915.csv")
#write.csv(dataALL,"Baseline/Data/DataAllBsl_110915.csv")
#write.csv(data,"Baseline/Data/DataSelBsl_110915.csv")

write.csv(initALL,"NoProtection/Init/InputAllNP_110915.csv")
write.csv(init,"NoProtection/Init/InputSelNP_110915.csv")
write.csv(dataALL,"NoProtection/Data/DataAllNP_110915.csv")
write.csv(data,"NoProtection/Data/DataSelNP_110915.csv")
