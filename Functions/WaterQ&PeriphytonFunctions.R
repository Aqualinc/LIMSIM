###################################################################################
###################################################################################

#                   WATER QUALITY AND PERIPHYTON FUNCTIONS                        #

###################################################################################
###################################################################################
# Functions indcluded in this script:



###################################################################################
###################################################################################
#setwd(dir$fun)
#source(paste(dir$fun,"PredMcDowell.r",paste="/"))
#setwd(dir$fun)
#source("REC_McDowell.r")
#setwd(dir$fun)
#load("FilamentModel.RData")

load(paste(dir$fun,"MeanMatsModel.RData",sep="/"))
load(paste(dir$fun,"MeanFilamentModel.RData",sep="/"))
load(paste(dir$fun,"MaxMatsModel.RData",sep="/"))
load(paste(dir$fun,"MaxFilamentModel.RData",sep="/"))

#load("MCIMcD.Rdata")  # MCIWQRF.Rdata  this model was based on the Unwin data
load(paste(dir$fun,"WqRFMods.Rdata",sep="/"))
load(paste(dir$fun,"hbMCIrf.Rdata",sep="/"))

############################################################################################################
#     Assign information for functions in this script
############################################################################################################
MyVars <<- c("NO3N", "DRP", "NH4N", "TN", "TP", "CLAR") # the variables needed for peri and MCI
varconv<<-c(3,3,3,3,3,0)  #This adjusts the units in log space.
names(varconv)<-MyVars

myCI<-list(Good=-1.96,Fair=0,Poor=1.96)

###################################################################################
#                            CALLING FUNCTION                                     #
###################################################################################Q

Run_WQandPeri<-function(MyREC=NULL,LandManagement="Fair",IrrigableAreaTarget=NULL,AllocQ=0,GWAlloc=0,PropIrri=0,SCCS=FALSE,dRegime=TRUE,WQModel="Unwin"){
  
  if (IrrigableAreaTarget>0){
    MyREC$usPastoral<-pmax(MyREC$usPastoral,pmin(1,MyREC$usIrriArea/MyREC$CATCHAREA*MyREC$IrriAreaIrrigated*IrrigableAreaTarget/100))}  
  
  MyREC$usPastoral[is.nan(MyREC$usPastoral)==T|is.na(MyREC$usPastoral)==T]<-0
      
      if (WQModel=="McDowell"){
          print("Calculating McDowell WQ estimates...")
          
          RECPastWQ <-read.csv(paste(dir$data,"BaseLineWQ.csv",sep="/"))
                               # store as an array
          Vars <- as.character(unique(RECPastWQ$Var)) # the water quality variables that we can model
          RECclass <- names(RECPastWQ)[c(4:21)]       # the REC classes that
                               
          meanWQ <<- array(0, dim=c(21,18,length(Vars)), dimnames = list(seq(0,100,5),RECclass,Vars))
          seWQ <<- array(0, dim=c(21,18,length(Vars)), dimnames = list(seq(0,100,5),RECclass,Vars))
                               
          RECPastWQ$Type <- paste(RECPastWQ$Var, RECPastWQ$est.or.se, sep=".")
                               
                               # load the data into the array!
          for(i in 1:length(Vars)) meanWQ[,,Vars[i]] <- as.matrix(RECPastWQ[RECPastWQ$Var == Vars[i] &  RECPastWQ$est.or.se == "est", c(4:21)])
          for(i in 1:length(Vars)) seWQ[,,Vars[i]] <- as.matrix(RECPastWQ[RECPastWQ$Var == Vars[i] &  RECPastWQ$est.or.se == "se", c(4:21)])
                               
          
          if (AllocQ>0&PropIrri==1){
            MyREC<-REC_McDowell(MyREC=MyREC,myCI=myCI[[LandManagement]],meanWQ=meanWQ)
          } else{
            MyREC<-REC_McDowell_scaled(MyREC=MyREC,myCI=myCI[[LandManagement]],PropIrri=PropIrri,meanWQ=meanWQ)
          }
      }else{
          print("Calculating Unwin WQ estimates...")
        
          MyREC<-REC_TonsWQ(MyREC=MyREC,myCI=myCI[[LandManagement]],PropIrri)
          MyREC$ECOLI<-10^(as.data.frame(predict(WqRFMods[["ECOLI"]], newdata = (MyREC))))
      }
     # browser()
  
  if (dRegime==TRUE){
  print("Calculating Mean and Annual Maximum, Filaments and Mats...")
      MyRECtemp<-MyREC;MyRECtemp$LowFlow<-MyRECtemp$BFIorig;MyRECtemp$FRE2<-MyRECtemp$FRE2orig
      MeanFils <-predict(MeanFilsModel, newdata = MyRECtemp);MeanFils[MeanFils<0]<-0
      MyREC$MeanFils <-MeanFils^2
      MaxFils<-predict(MaxFilsModel, newdata = MyRECtemp);MaxFils[MaxFils<0]<-0
      MyREC$MaxFils<- MaxFils^2
      MeanMats<-(predict(MeanMatsModel,newdata=MyRECtemp));MeanMats[MeanMats<0]<-0
      MyREC$MeanMats<-MeanMats^2
      MaxMats<-(predict(MaxMatsModel,newdata=MyRECtemp));MaxMats[MaxMats<0]<-0
      MyREC$MaxMats<-MaxMats^2}
  print("Calculating Mean and Annual Maximum, Filaments and Mats with changed flow regime...")
      MyRECtemp2<-MyREC;MyRECtemp2$LowFlow<-MyREC$BFI;
      MyREC$MeanFilsDeltaQ <-(predict(MeanFilsModel, newdata = MyRECtemp2))^2 
      MyREC$MaxFilsDeltaQ<-(predict(MaxFilsModel, newdata = MyRECtemp2))^2
      MyREC$MeanMatsDeltaQ<-(predict(MeanMatsModel,newdata=MyRECtemp2))^2
      MyREC$MaxMatsDeltaQ<-(predict(MaxMatsModel,newdata=MyRECtemp2))^2 
  
  #Calculate MCI Scores
  print("Calculating MCI scores...")
      #MCIUseVars <- row.names(importance(MCI.WQMcDreduce)) # to avoid passing unused predictors in MyREC that have missing values
      #MyREC$MCI <-  predict(MCI.WQMcDreduce, newdata = MyREC[,MCIUseVars])
  MyRECtemp3<-MyREC;  MyRECtemp3$usAveSlope<-MyRECtemp3$usSlope;  MyRECtemp3$FRE3.Count<-MyRECtemp3$FRE3
  MyREC$MCI<-REC_TonsMCI(Data=MyRECtemp3,myCI=myCI[[LandManagement]],PropIrri,AllocQ)
  
  MyRECtemp4<-MyREC;  MyRECtemp4$usAveSlope<-MyRECtemp4$usSlope;  MyRECtemp4$FRE3.Count<-MyRECtemp4$FRE3
  MyREC$MCIdeltaQ<-REC_TonsMCI(Data=MyRECtemp4,myCI=myCI[[LandManagement]],PropIrri,AllocQ)
  return(MyREC)
}


###################################################################################
###################################################################################
REC_TonsMCI<-function(Data=NULL,myCI=0,PropIrri=0,AllocQ=NULL){
  #browser()
  if (PropIrri==1){
    myCItemp<-myCI*(min(Data$usPastoral,Data$usIrriArea/Data$CATCHAREA))/Data$usPastoral
    myCItemp[is.na(myCItemp)]<-0
    ind<-which(MyREC$TheTake>0)
    myCI<-matrix(data=myCI,nrow=length(Data[,1]),ncol=1)
    myCI[ind]<-myCItemp[ind]
  }else{
    myCI<-matrix(data=myCI,nrow=length(Data[,1]),ncol=1)
    ind<-which(MyREC$TheTake==0)
    myCI[ind]<-0
  }
    
  temp<-predict(hbMCIrf, newdata = Data[,row.names(importance(hbMCIrf))])
  se<-hbMCIrf$RMSD*myCI
  MCI<-(temp+se)
  
  return(MCI)
}



REC_TonsWQ<-function(MyREC=NULL,myCI=0,PropIrri=0,AllocQ=0){
#browser()
  if (PropIrri==1){
    myCItemp<-myCI*(min(MyREC$usPastoral,MyREC$usIrriArea/MyREC$CATCHAREA))/MyREC$usPastoral
    myCItemp[is.na(myCItemp)]<-0
    ind<-which(MyREC$TheTake>0)
    myCI<-matrix(data=myCI,nrow=length(MyREC[,1]),ncol=1)
    myCI[ind]<-myCItemp[ind]
  }else{
    myCI<-matrix(data=myCI,nrow=length(MyREC[,1]),ncol=1)
    ind<-which(MyREC$TheTake==0)
    myCI[ind]<-0
  }
  print(c(max(myCI),min(myCI)))
  #browser()
  temp<-Doug.cbind.list(sapply(MyVars,predTonsWQ,myCI=myCI,Data=MyREC))
  names(temp)<-paste("log10",MyVars,sep="")
  temp2<-(10^temp)
  colnames(temp2)<-MyVars
  MyREC<-cbind(MyREC,temp,temp2)
  MyREC$log10DINDRP <- log10((MyREC$NO3N + MyREC$NH4N)/MyREC$DRP)
  MyREC$log10TNTP<-log10(MyREC$TN/MyREC$TP)
  MyREC$SIN<-(MyREC$NO3N + MyREC$NH4N)
  return(MyREC)
}

predTonsWQ<-function(myWQvar = "DRP", myCI = 0,Data=MyREC) {

  temp<-as.data.frame(predict(WqRFMods[[myWQvar]], newdata = (Data)))
  se<-WqRFMods[[myWQvar]]$RMSD*myCI
  WQ<-(temp+se)+varconv[myWQvar]
  return(WQ)
}


REC_McDowell<-function(MyREC=NULL,myCI=0,meanWQ=NA){
  MyREC$NO3N<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "NO3N", myCI = 0,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)) * 1000
  MyREC$DRP<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "DRP", myCI = 0,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)) * 1000
  MyREC$log10DRP<-log10(MyREC$DRP)
  MyREC$NH4N<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "NH4N", myCI = 0,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)) * 1000
  
  MyREC$log10DINDRP <- log10((MyREC$NO3N + MyREC$NH4N)/MyREC$DRP)
  MyREC$TN <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "TN", myCI = myCI,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)) * 1000
  MyREC$log10TN  <-  log10(MyREC$TN)
  MyREC$TP  <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "TP", myCI = myCI,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)) * 1000
  MyREC$log10TNTP<-log10(MyREC$TN/MyREC$TP)
  MyREC$SIN  <- 1000 * (unlist(lapply(pick, predMcDowellWQ, myWQvar = "NO3N", myCI = myCI,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC))
                        + unlist(lapply(pick, predMcDowellWQ, myWQvar = "NH4N", myCI = myCI,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC)))  
  MyREC$CLAR  <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "CLAR", myCI = -myCI,RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC))
  return(MyREC)
}

REC_McDowell_scaled<-function(MyREC=NULL,myCI=0,PropIrri=0,meanWQ=NA){
  
  if (PropIrri==1){
    myCItemp<-myCI*(min(Data$usPastoral,Data$usIrriArea/Data$CATCHAREA))/Data$usPastoral
    myCItemp[is.na(myCI)]<-0
    ind<-which(MyREC$AllocQ>0)
    myCI<-matrix(data=myCI,nrow=length(Data[,1]),ncol=1)
    myCI[ind]<-myCItemp[ind]
  }
  
  MyREC$NO3N<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NO3N", myCI = 0,Data=MyREC)) * 1000
  MyREC$DRP<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "DRP", myCI = 0,Data=MyREC)) * 1000
  MyREC$log10DRP<-log10(MyREC$DRP)
  MyREC$NH4N<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NH4N", myCI = 0,Data=MyREC)) * 1000
  
  MyREC$log10DINDRP <- log10((MyREC$NO3N + MyREC$NH4N)/MyREC$DRP)
  MyREC$TN <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "TN", myCI = myCI,Data=MyREC)) * 1000
  MyREC$log10TN  <-  log10(MyREC$TN)
  MyREC$TP  <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "TP", myCI = myCI,Data=MyREC)) * 1000
  MyREC$log10TNTP<-log10(MyREC$TN/MyREC$TP)
  MyREC$SIN  <- 1000 * (unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NO3N", myCI = myCI,Data=MyREC)) + unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NH4N", myCI = myCI)))  
  MyREC$CLAR  <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "CLAR", myCI = -myCI,Data=MyREC))
  return(MyREC)
}

############################################################################################################
# function to compute any WQ concentration + or - error for any NZReach based on McDowell et al.
predMcDowellWQ <- function(ThisNZReach = 13524724,  myWQvar = "DRP", myCI = 0, RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC) {
  mySoF = as.character(Data[Data$NZReach == ThisNZReach, "SoF"])
  if(is.na(match(mySoF, RECclasses))) { ThisPred <- NA # a class may not have predictions available in the array
  } else {
    percPasture <- Data[Data$NZReach == ThisNZReach, "usPastoral"]
    log10Median <- approx(y=meanWQ[,mySoF,myWQvar], x=seq(0,1,0.05), xout=percPasture)$y # interpolate median value
    log10SE <- approx(y=seWQ[,mySoF,myWQvar], x=seq(0,1,0.05), xout=percPasture)$y # interpolate SE value
    # if myCI = 0 then the median is returned otherwise the median + myCI * the SE
    if(myCI == 0|percPasture<0.01) { ThisPred <- 10^log10Median
    } else {
      upSE <- 10^(log10Median + log10SE) - 10^log10Median
      downSE <- 10^(log10Median) - 10^(log10Median - log10SE)
      if(myCI > 0)  { ThisPred <- 10^log10Median + myCI*upSE
      } else {  ThisPred <- 10^log10Median + myCI*downSE
      }
    }
  }
  return(ThisPred)
}

# function to compute any WQ concentration + or - error for any NZReach based on McDowell et al.
predMcDowellWQ_scaled <- function(ThisNZReach = 13524724,  myWQvar = "DRP", myCI = 0, RECclasses = unlist(dimnames(meanWQ)[2]),Data=MyREC) {
  mySoF = as.character(Data[Data$NZReach == ThisNZReach, "SoF"])
  Data$usIrriArea[is.nan(Data$usIrriArea)==T|is.na(Data$usIrriArea)==T]<-0
  myCI<-myCI*(min(Data$usPastoral[ThisNZReach],Data$usIrriArea[ThisNZReach]/Data$CATCHAREA[ThisNZReach]))/Data$usPastoral[ThisNZReach]
  if (is.na(myCI)==T) myCI<-0
  if(is.na(match(mySoF, RECclasses))) { ThisPred <- NA # a class may not have predictions available in the array
  } else {
    percPasture <- Data[Data$NZReach == ThisNZReach, "usPastoral"]
    log10Median <- approx(y=meanWQ[,mySoF,myWQvar], x=seq(0,1,0.05), xout=percPasture)$y # interpolate median value
    log10SE <- approx(y=seWQ[,mySoF,myWQvar], x=seq(0,1,0.05), xout=percPasture)$y # interpolate SE value
    # if myCI = 0 then the median is returned otherwise the median + myCI * the SE
    if(myCI == 0) { ThisPred <- 10^log10Median
    } else {
      upSE <- 10^(log10Median + log10SE) - 10^log10Median
      downSE <- 10^(log10Median) - 10^(log10Median - log10SE)
      if(myCI > 0)  { ThisPred <- 10^log10Median + myCI*upSE
      } else {  ThisPred <- 10^log10Median + myCI*downSE
      }
    }
  }
  return(ThisPred)
}