###################################################################################
###################################################################################

#                   FLOW AND ABSTRACTION FUNCTIONS                                #

###################################################################################
###################################################################################
# Functions indcluded in this script:
#       - Run_FlowandAbstraction   Controlling function that adds all flow and abstraction variables to MyREC
#       - GenFDC        Generates flow duraction curve
#       - TakeRate      Determines the abstraction for a given reach 
#       - freq.restrict Calculates the frequency of restrictions and the reliability


###################################################################################
###################################################################################


###################################################################################
#CONTROLLING FUNCTION THAT CREATES ALL OF THE FLOW AND ABSTRACTION VARIABLES IN MY REC
###################################################################################
Run_FlowandAbstraction<-function(MyREC=MyREC,AllocQ=AllocQ,MinQ=MinQ,pick=pick,TakeAll=0,SysCap=0.58,EffIrriArea=0.8){
  
  
  #setwd(dir$fun)
  load(paste(dir$fun,"DeltaInd_v2.RData",sep="/"))
  
  # generate FDC for all NZReach.
  Perc <<- seq(0,1, length = 101) # Perc sets percentiles for which flows are to be generated at - needs to be set globally for other functions (this is the default!)
  
  print("Calculating flow duration curves....")
  MyFDC <<- Doug.rbind.list(lapply(pick, GenFDC, Data = MyREC)) # generate FDC for all locations (Function inside "LIMSIM Functions.R")
  names(MyFDC) <- paste("P", 100*Perc, sep="") #
  
  print("Calculating irrigation takes and reliability...")
      if (length(AllocQ==1)){  
        MyREC$MinQ<-matrix(MinQ,length(pick),1); MyREC$AllocQ<-matrix(AllocQ,length(pick),1)
      }else{
        MyREC$MinQ<-MinQ; MyREC$AllocQ<-AllocQ
      }  
  
  MyREC$BFIorig<-MyREC$BFI;MyREC$FRE2orig<-MyREC$FRE2;MyREC$FRE3orig<-MyREC$FRE3;MyREC$nNegorig<-MyREC$nNeg

  
  MyREC$TheTake <-0; MyREC$IrriAreaIrrigated<-NA; MyREC$freq.man <-NA; MyREC$freq.min <-NA; MyREC$freq.reliability <-NA;MyREC$AllocQUse<-NA
  
  ind1<-which(MyREC$AllocQ>0)  #Determine those reaches where there is some allocaiton
  if (length(ind1)>0){
  MyRECtemp<-MyREC[ind1,];pick2<-pick[ind1]
  MyRECtemp$TheTake <- sapply(pick2, TakeRate, RegAllocation=NULL, Data = MyRECtemp,TakeAll=TakeAll)  
  MyRECtemp$AllocQUse<-MyRECtemp$TheTake/(MyRECtemp$MALF*MyRECtemp$AllocQ)
  MyRECtemp$AllocQUse[MyRECtemp$TheTake==0]<-NA
  MyRECtemp$IrriAreaIrrigated<-NA
  ind<-which(MyRECtemp$TheTake>0)
  
  MyRECtemp$IrriAreaIrrigated[ind]<-MyRECtemp$TheTake[ind]*1000*10000/SysCap/(MyRECtemp[ind,"usIrriArea"]*EffIrriArea)*100
  
  Restriction2 <- Doug.rbind.list(lapply(pick2, freq.restrict, prop=NULL, RegAllocation=NULL,Plot=F, Data = MyRECtemp,Data1=MyREC)) # turn Plot=F for speed.
  MyRECtemp$freq.man <- Restriction2$freq.man; MyRECtemp$freq.min <- Restriction2$freq.min; MyRECtemp$freq.reliability <- Restriction2$freq.reliability
  MyRECtemp$freq.man[MyRECtemp$TheTake==0]<-NA; MyRECtemp$freq.min[MyRECtemp$TheTake==0]<-NA; MyRECtemp$freq.reliability[MyRECtemp$TheTake==0]<-NA;
  
  MyRECtemp$BFI<-MyRECtemp$BFI+predict(dBFImod,newdata=MyRECtemp)
  MyRECtemp$BFI[MyRECtemp$TheTake==0]<-MyRECtemp$BFIorig[MyRECtemp$TheTake==0]
  
  MyRECtemp$FRE2<-MyRECtemp$FRE2+predict(dFRE2mod,newdata=MyRECtemp)
  MyRECtemp$FRE2[MyRECtemp$TheTake==0]<-MyRECtemp$FRE2orig[MyRECtemp$TheTake==0]
  
  MyRECtemp$FRE3<-MyRECtemp$FRE3+predict(dFRE3mod,newdata=MyRECtemp)
  MyRECtemp$FRE3[MyRECtemp$TheTake==0]<-MyRECtemp$FRE3orig[MyRECtemp$TheTake==0]
  
  MyRECtemp$nNeg<-MyRECtemp$nNeg+predict(dnNegMOD,newdata=MyRECtemp)
  MyRECtemp$nNeg[MyRECtemp$TheTake==0]<-MyRECtemp$nNegorig[MyRECtemp$TheTake==0]
  
  MyREC[ind1,]<-MyRECtemp
  }
  #browser()
return(MyREC)
}
###################################################################################
###################################################################################
# function for generating FDC (m3/s) from GEV distribution (from R package nsRFA)
GenFDC <- function(ThisNZReach = 13524724, P = Perc, Data = MyREC) {
  
  if((length(P))>1) {P[length(P)] <- 0.995}  # re-set highest probability (to avoid generating infinite Q) (Note highest Q is extremely sensitive to this value)
  ThisRow <- which(Data$NZReach == ThisNZReach) # get the parameters from the FDC dataset is this faster???
  xi  <- Data[ThisRow, "XI"] # the parameters  xi, alfa, k will be found in MyFDCPars
  alfa <- Data[ThisRow, "ALFA"]
  k    <- Data[ThisRow, "K"]
  
  #browser()
  if (is.na(k)==F){
  if ((k > -1e-07) & (k < 1e-07)) {
    Q <- xi - alfa * log(-log(P))
  }
  else {
    Q <- xi + alfa * (1 - (-log(P))^k)/k
  }
  Q[Q<0] <- 0   # no values should be negative
  Q <- Q * Data[ThisRow, "MeanFlow"]   # multiply by estimated mean flow
  }else{
    Q<-P*NA
  }
  return(Q)
}

###############################################################################
# compute the amount of take from any location
###############################################################################
TakeRate <- function(ThisNZReach = 13524724, RegAllocation=0.5, Data = 0, UseRate = 0.58,TakeAll=TakeAll) {  # use rate = irrigation use l/s/ Ha
  #This function determines the Take rate or a certain reach, which is the smaller of the Allocation and the maximum demand 
  ThisRow <- which(Data$NZReach == ThisNZReach) # the row in the dataset for this NZREach
  #browser()
  if (is.null(RegAllocation)==T){RegAllocation<-Data$AllocQ[ThisRow]}
  if(Data[ThisRow,"usIrriArea"] > 10000&!is.na(Data[ThisRow,"usIrriArea"])==T) {   
    MaxDemand <-  Data[ThisRow, "usIrriArea"]/10000 * UseRate/1000*0.8 #The upstream irrigable area multiplied by the     
    RegulatedTake <-  RegAllocation * Data[ThisRow, "MALF"]  # the "Regulatory take" i.e. the allowed according to a % of MALF rule (MALf from the hyd predictions data)
    Take <- min(MaxDemand, RegulatedTake)
    if (TakeAll==1) {Take<-RegulatedTake}
  } else {
    Take <- 0  # we assume the landscape is too steep to irrigate No water  is taken
  }
  return(Take)}

#######################################################################################
#   calculate  freq.restrictions based on MIN FLOW and ALLOCATION rate for a GAUGE    #
#######################################################################################

freq.restrict <- function(ThisNZReach = 13524724, FDC = MyFDC,  minFlow=NULL, prop=0.8, Qref=NULL, RegAllocation=0.5, Plot = TRUE, Data = 0,Data1=0) {
  # the function will asume that allocations are based on MALF unless given
  # prop is the proportion of MALF for Qmin and allocate is the proportion of MALF to allocate
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach) # get the parameters from the FDC dataset is this faster???
  ThisRow_ALL <- which(Data1$NZReach == ThisNZReach)
  if (is.null(prop)==TRUE){prop<-Data$MinQ[ThisRow]}
  if (is.null(RegAllocation)==T){RegAllocation<-Data$AllocQ[ThisRow]}
  #if(Data$TheTake[ThisRow] > 0) 
  #browser()
  if(Data$TheTake[ThisRow] == 0|is.na(Data$TheTake[ThisRow])==T) {  # if theere is no take skip entirely frequency of restriction is ZERO NA, is.na(Data$TheTake[ThisRow])==T
    freq.man <- NA
    freq.min <- NA
    freq.reliability<-NA
  } else {
   # browser()
    if (any(is.na(FDC[ThisRow, ]))==T){
      freq.man <- NA
      freq.min <- NA
      freq.reliability<-NA
    }else{
    #
    if (is.null(Qref)) Qref <- Data[ThisRow, "MALF"]  # the reference flow is the MALf from the hyd predictions data
    if (is.null(minFlow)) minFlow <- Qref * prop   # minflow = flow at which there is TOTAL restriction
    
    if(!is.null(Data$TheTake[ThisRow]))   {
      allocation <-  Data$TheTake[ThisRow] # if a take for each reach is defined by using the function "TakeRate" then use
    } else {
      allocation <- RegAllocation*Qref   #  OTHERWISE a regulatory rule of thumb allocation is used Qref*allocate. The DEFAULT for Qref= MALF, (NES type rule)
    }
    
    if(allocation == 0) minFlow <- 0 # if there are no takes there should be no restrictions (min=0)                     
    manFlow <- minFlow+allocation #  minflow+allocation = "management flow"  the flow at which restrictions begin
    FDC.data <- FDC[ThisRow_ALL, ]   # THIS is the FDC
    freqs <- 100*Perc # the percentiles the FDC is estimated for from the global value of Perc              
    freq.diff <-  approx(y=freqs, x=FDC.data, xout=c(manFlow, minFlow))$y # interploate W vs Q data
    freq.man <- freq.diff[1]
    freq.min <- freq.diff[2]     # minflow= flow at which there is TOTAL restriction   
    FDC_altered<-FDC.data
    FDC.take<-approx(x=c(100,freq.diff,0),y=c(allocation,allocation,0,0),xout=freqs)
    freq.reliability<-mean(FDC.take$y)/allocation*100  #NOTE THIS ONLY WORKS IF FREQS iS EVENLY SAMPLED
    #if(is.na(freq.min)) freq.min <- 100  # the minFlow can be OFF the end of the FDC curve! (especially when MALF is annual but FDC is for summer.
    #if(is.na(freq.man)) freq.man <- 100  # the FREQ manFlow can be NA in MALF = 0                                        
  }
  }
  
  if(Plot==TRUE)   {
    par(mfrow=c(2,1), bg="grey90")
    plot(freqs, log10(FDC.data), col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="log10 Flow")
    abline(h=log10(minFlow), col = "green")
    abline(v=freq.min, col = "green", lty=2)
    abline(h=log10(manFlow), col = "red")
    abline(v=freq.man, col = "red", lty=2)
    legend("top", cex=0.8,
           legend=c(paste("Minimum Flow",round(minFlow,2)),
                    paste("Freq Min Flow", round(freq.min, 2)),
                    paste("Management Flow", round(manFlow,2)),
                    paste("Freq Man Flow",round(freq.man,2))),
           col = c("green","green", "red", "red"),   text.col = "black", lty = c(1, 2,1,2), pch = c(-1, -1),   lwd = c(1,1,1,1),
           merge = TRUE, bg = 'white', inset = .05)
    plot(freqs, FDC.data, col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="Flow",
         xlim=c(freq.man*0.9, 100), ylim=c(0, manFlow*1.1))
    abline(h=minFlow, col = "green")
    abline(v=freq.min, col = "green", lty=2)
    abline(h=manFlow, col = "red")
    abline(v=freq.man, col = "red", lty=2)
  } #end plot
  
  return(cbind(freq.man, freq.min,freq.reliability))   #returns two values 1 = % of time restrition 2 % of time TOTAL restriction (note the variable names (freq) is poor use of words these are not frequencies!
}# end

#######################################################################################
#######################################################################################
MonthString  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
GetDaily3 <- function(Daily,Date=NA) {  # this has the natural median as a column
  # make factors for time-series

  if (is.na(Date[1])==F){Daily$Date<-Date}
  Daily$chrono <- chron(as.character(Daily$Date), format = c(dates = "dd-mmm-yy"), out.format = c(dates = "d-m-y"))
  Daily$Dateformat <- as.Date(Daily$chrono, "%d/%m/%y")
  Daily$Julian <- julian(Daily$Dateformat)

  DayMonthYear <- unlist(sapply(as.character(Daily$Date), strsplit, split = "-"))
  Daily$day <- as.numeric(DayMonthYear[seq(1, length(DayMonthYear), by = 3)])
  Daily$month <- factor(DayMonthYear[seq(2, length(DayMonthYear), by = 3)], levels = MonthString)
  Daily$year <- as.numeric(DayMonthYear[seq(3, length(DayMonthYear), by = 3)])
  Daily$YYYY <- Daily$year
  
  Daily$Date1Jan <- as.Date(paste(as.character(Daily$YYYY), "01", "01", sep = "-") )
  Daily$Julian1Jan <- julian(Daily$Date1Jan)
  Daily$Julian <- julian(Daily$Dateformat)
  Daily <- Daily[order(Daily$Julian), ]
  Daily$DayOfYear <- Daily$Julian - Daily$Julian1Jan + 1
  
  names(Daily)[which(names(Daily) == "Data")] <- "Daily"
  
  # Get day to day changes in flow
  Daily$Daily <- as.numeric(Daily$Daily)
  Daily$DailyChange <- c(diff(Daily$Daily), NA )
  Daily$NaturalMedian <- median(Daily$Daily, na.rm=T)  # put the natural median into the frame for later use
  
  return(Daily)
}

#######################################################################################
#######################################################################################
GetPeriHydPreds2  <- function(Daily, Durations=7, LowThreshold = NA, HighThreshold = NA, # version with nNeg and NO low Flow
                              StartDay = c(min = 183, max = 1),NaturalMean=mean(Daily$Daily,na.rm=T),
                              # WaterYearDefinition = c(JulianMin = "StartsJune", JulianMax = "StartsDecember"),
                              ...) {
  
  #print(Daily$Date[1])
  
  Daily$YYYY <- as.factor(Daily$YYYY)
  
  ##### FreX (X = 3?) in each year with a Y (7 day) windows
  FREMedianMultiples = c(2,3,4,1.5)
  names(FREMedianMultiples) <- FREMedianMultiples
  FREXList <- lapply(FREMedianMultiples, GetLengthsBelowThreshold2, myFrame = Daily[, c("Daily","YYYY", "NaturalMedian")]) # use myWindow to change window size NOTE median set to natural!!!
  names(FREXList) <- paste("FRE", names(FREXList), sep = "")
  FREX <- Doug.cbind.list(FREXList)
  
  
  RatesOfChange <- Doug.rbind.list(tapply(Daily$DailyChange, Daily$YYYY, GetPositive, ...))
  Reversals <- tapply(Daily$DailyChange, Daily$YYYY, GetReversal)
  
  # get number of days without data in each year
  #GapDays <- tapply(Daily$Daily, Daily$YYYY, function(x)sum(is.na(x)))
  # length of record, number of consecutive days
  nDays <- tapply(Daily$DayOfYear, Daily$YYYY, function(x) max(x) - min(x) + 1)
  
  
  ##############################################################
  #              Group1: Get monthly values                    #
  ##############################################################
  
  # get running means with these windows
  RunningMeans <- sapply(Durations, GetRunningMean, myDaily = Daily$Daily, na.rm = T)
  colnames(RunningMeans) <- paste("Mean", as.character(Durations), "DayFlow", sep = "")
  RunningMeans <- data.frame(RunningMeans)
  # get max and min for each year of running means for different poeriods
  Yearly <- t(Doug.rbind.list(lapply(RunningMeans, GetExtremes, myYear = Daily$YYYY)))
  colnames(Yearly) <- paste(rep(colnames(RunningMeans), rep(2, ncol(RunningMeans))), colnames(Yearly), sep = "")
  
  #ZeroFlowDays <- tapply(Daily$Daily, Daily$YYYY, function(x)sum(x[!is.na(x)] == 0))
  MeanFlow <- tapply(Daily$Daily, Daily$YYYY, mean, na.rm = T)
  
  BFI <- Yearly[ ,"Mean7DayFlowMins"] / NaturalMean # NOTE this function sets the median to be the natural median for the calculation of FRE3
  
  ##############################################################
  #               Get summaries of all RVA stats               #
  ##############################################################
  
  # all RVA data for each year        #
  #ZeroAndBFI <- cbind(BFI)
  GapsAndNum <-  cbind(RatesOfChange, Reversals, nDays)
  FRE2 <- (cbind(FREX[, "FRE2.Count"]))
  FRE3 <- (cbind(FREX[, "FRE3.Count"]))
  FRE4 <- (cbind(FREX[, "FRE4.Count"]))
  FRE1_5 <- (cbind(FREX[, "FRE1.5.Count"]))
  AnnualRVA <- cbind.data.frame(FRE2,FRE3,FRE4,FRE1_5, GapsAndNum,BFI)
  names(AnnualRVA) <- c( "FRE2","FRE3","FRE4","FRE1_5", "nPos", "nNeg", "Reversals", "nDays","BFI")
  return(AnnualRVA)
  # return(list(Monthly, Yearly, ExtremeTiming, DurPulses, nPulses, RatesOfChange, GapDays, nDays))
}
#################################################################################
#################################################################################
GetBFI  <- function(Daily, Durations=7, LowThreshold = NA, HighThreshold = NA,
                    StartDay = c(min = 183, max = 1)){
  
  Daily$YYYY <- as.factor(Daily$YYYY)
  
}
#################################################################################
#################################################################################
# function to calc running mean with a x day window
GetRunningMean <- function(myWindow, myDaily, ...) {
  myRunningMean <- running(myDaily, Y=NULL, fun=mean, width=myWindow, allow.fewer = T, ...)
  return(myRunningMean)
}

#################################################################################
#################################################################################
# function to calculate annual min and max for each time-series
GetExtremes <- function(myRunningMean, myYear) {
  Mins <- tapply(myRunningMean, myYear, DougMin, na.rm = T)
  Maxs <- tapply(myRunningMean, myYear, DougMax, na.rm = T)
  myExtremes <- t(cbind(Mins, Maxs))
  
  return(myExtremes)
}
#################################################################################
#################################################################################
GetLengthsBelowThreshold2 <- function(myMedianMultiple, myFrame, myWindow = 5, myMedianFlow){
  RealOnly <- !is.na(myFrame$Daily)                                                                                                        
  DaysBetweenFRE3list <- tapply(myFrame$Daily[RealOnly], myFrame$YYYY[RealOnly], GetLengthsBelowThreshold, 
                                myThreshold = myMedianMultiple * myFrame$NaturalMedian[1]) #  THIS takes the NATURAL MedianFlow for computing FRE2
  
  # dont count periods between that are less the an filter size
  DaysBetweenFRE3list <- lapply(DaysBetweenFRE3list, function(x) x[x>myWindow])
  # print(DaysBetweenFRE3list)  # count number of events
  Count <- sapply(DaysBetweenFRE3list, function(x) ifelse(is.null(x),NA,length(x)))
  
  # what was the flow at the start of the year
  StartOfYearFlow <- tapply(myFrame$Daily[RealOnly], myFrame$YYYY[RealOnly], function(x) x[1])
  # if all of year is missing then StartOfYearFlow will be NA
  StartOfYearFlow[is.na(StartOfYearFlow)] <- 0
  # which years start with flows lower than the threshold
  #LowStarts <-  StartOfYearFlow  < myMedianMultiple *  median(myFrame$Daily, na.rm = T)# THIS NEEDS TO BE REPLACED BY NATURAL MedianFlow for scenarios
  LowStarts <-  StartOfYearFlow  < myMedianMultiple *  myFrame$NaturalMedian[1]
  # if starting below threshold, then subtract one from count
  Count[LowStarts] <- Count[LowStarts] - 1
  
  MeanDurBetween <- sapply(DaysBetweenFRE3list, function(x) ifelse(is.null(x),NA,mean(x, na.rm = T)))          # mean duration between events
  MaxDurBetween <- sapply(DaysBetweenFRE3list,  function(x) ifelse(length(x)==0,NA,max (x, na.rm = T)))          # max duration between events
  StDevDurBetween <- sapply(DaysBetweenFRE3list,function(x) ifelse(is.null(x),NA,sd  (x, na.rm = T)))          # st dev duration between events
  myOut <- data.frame(Count, MeanDurBetween, MaxDurBetween, StDevDurBetween)
  return(myOut)
}
#################################################################################
#################################################################################
# get number of times direction of rate of flow change changes
GetReversal <- function(x){
  myCount <- 0
  x <- x[!is.na(x)]
  if (length(x) > 2) {
    myTrend <- x[1] > 0
    for (n in 2:length(x)){
      myDir <- x[n] > 0
      if (myDir != myTrend) {
        myCount <- myCount + 1
        myTrend <- myDir
      }
    }
  }
  return(myCount)
}

#################################################################################
#################################################################################
# function to compute alf for stations
GetALF7  <- function(Daily, Durations=7, LowThreshold = NA, HighThreshold = NA,
                     StartDay = c(min = 183, max = 1),
                     # WaterYearDefinition = c(JulianMin = "StartsJune", JulianMax = "StartsDecember"),
                     ...) {
  #browser()
  #print(Daily$Date[1])
  Daily$YYYY <- as.factor(Daily$YYYY)
  
  # get running means with these windows
  RunningMeans <- sapply(Durations, GetRunningMean, myDaily = Daily$Daily, na.rm = T)
  colnames(RunningMeans) <- paste("Mean", as.character(Durations), "DayFlow", sep = "")
  RunningMeans <- data.frame(RunningMeans)
  # get max and min for each year of running means for different poeriods
  Yearly <- t(Doug.rbind.list(lapply(RunningMeans, GetExtremes, myYear = Daily$YYYY)))
  colnames(Yearly) <- paste(rep(colnames(RunningMeans), rep(2, ncol(RunningMeans))), colnames(Yearly), sep = "")
  
  ALF <- Yearly[ ,"Mean7DayFlowMins"]
  
  # get number of days  data in each year
  #nDays <- tapply(Daily$DayOfYear, Daily$YYYY, function(x) max(x) - min(x) + 1)
  nDays<-sapply(unique(Daily$YYYY),nDayscalc,data=Daily)
  names(nDays)<-unique(Daily$YYYY)
  
  AllALF <- cbind(ALF)
  GapsAndNum <-  cbind(nDays)
  AnnualRVA <- cbind.data.frame(AllALF, GapsAndNum)
  names(AnnualRVA) <- c("ALF", "nDays")
  
  return(AnnualRVA)
}
#################################################################################
#########

nDayscalc<-function(YY=2000,data=Daily){
  NoD<-length(which(is.na(data$Daily)==FALSE&data$YYYY==YY))
  
}

###################################################################################
DougMin <- function(myVec, ...){
  if (all(is.na(myVec))) myMin <- NA
  else myMin <- min(myVec, ...)
  return(myMin)
}
DougMax <- function(myVec, ...){
  if (all(is.na(myVec))) myMin <- NA
  else myMin <- max(myVec, ...)
  return(myMin)
}

#################################################################################
#################################################################################
GetPositive <- function(x, SplitSteadyFlows = F){
  mySteady <- ifelse(SplitSteadyFlows, sum(x == 0, na.rm = T)/2, 0)
  myPos <- (x > 0)
  nPos <- sum(myPos, na.rm = T) + mySteady
  #meanPos <- mean(x[myPos], na.rm = T)
  #medianPos <- median(x[myPos], na.rm = T)
  #ShapiroPos <- GetShapiro(x[myPos])[2]
  
  myNeg <- (x <= 0)
  nNeg <- sum(myNeg, na.rm = T) + mySteady
  #meanNeg <- mean(x[myNeg], na.rm = T)
  #medianNeg <- median(x[myNeg], na.rm = T)
  #ShapiroNeg <- GetShapiro(x[myNeg])[2]
  
  #QuantCheck = (nPos * meanPos) + (nNeg * meanNeg) # check on continuity
  
  #return(data.frame(nPos, meanPos, medianPos, ShapiroPos, nNeg, meanNeg, medianNeg, ShapiroNeg))
  return(data.frame(nPos,  nNeg))
}
GetLengthsBelowThreshold <- function(myVec, myThreshold) {
  myCount <- 0
  myLength <- 0
  IsNA <- is.na(myVec)
  if (all(IsNA) | (length(myVec[!IsNA]) < 2) ) {
    myCount <- NA
    myLength <- NA
  } else {
    myVec <- myVec[!IsNA]
    if(myVec[1] < myThreshold) {
      myCount <- myCount + 1
      myLength[myCount] <- 1
    }
    for (n in 2:length(myVec)) {
      if( (myVec[n] < myThreshold) & (myVec[n - 1] >= myThreshold) ) {
        myCount <- myCount + 1
        myLength[myCount] <- 0
      }
      if (myVec[n] < myThreshold) {
        myLength[myCount] <- myLength[myCount] + 1
      }
    }
  }
  return(myLength)
}

#################################################################################
#################################################################################
GetFlow <- function(x, Min=20, Alloc=20, Gap=0, AllocB=0) {
  
  A <- Alloc
  B <- AllocB
  
  #browser()
  if(is.na(x)) {
    myOut <- NA
  } else {
    
    if(x <= Min) myOut <- x        # below minium
    if(x > Min & x < (A+Min)) {   # within a block
      # print("within a block")
      myAvailable <- x - Min
      myOut <- x - myAvailable
    }
    if(x > (A+Min) & x <= (A+Min+Gap) ) {            # within gap
      #print("within a gap")
      myAvailable <-  A
      myOut <- x - myAvailable
    }
    if( x > (A+Min+Gap) & x <= (A+Min+Gap+B)) {            # within b
      #print("within b")
      myAvailable <-  x - Min - A - Gap # A + B
      myOut <- x - A - myAvailable
    }
    if(x >= (A+Min+Gap+B))  {
      # print("over b")
      myAvailable <-  B
      myOut <- x - A - myAvailable
    }
  }
  if (exists("myOut")==FALSE){browser()}
  return(as.numeric(myOut))
}

#################################################################################
#################################################################################
GetFlowTS <- function(x, Min=20, Alloc=20, Gap=0, AllocB=0,UnDem=NULL) {
  #Same as get flow, but can insert a unit demand timeseries
  #browser()
  A <- Alloc*UnDem
  B <- AllocB*UnDem
  
  #browser()
  if(is.na(x)) {
    myOut <- NA
  } else {
    
    if(x <= Min) myOut <- x        # below minium
    if(x > Min & x < (A+Min)) {   # within a block
      # print("within a block")
      myAvailable <- x - Min
      myOut <- x - myAvailable
    }
    if(x > (A+Min) & x <= (A+Min+Gap) ) {            # within gap
      #print("within a gap")
      myAvailable <-  A
      myOut <- x - myAvailable
    }
    if( x > (A+Min+Gap) & x <= (A+Min+Gap+B)) {            # within b
      #print("within b")
      myAvailable <-  x - Min - A - Gap # A + B
      myOut <- x - A - myAvailable
    }
    if(x >= (A+Min+Gap+B))  {
      # print("over b")
      myAvailable <-  B
      myOut <- x - A - myAvailable
    }
  }
  if (exists("myOut")==FALSE){browser()}
  #if (is.na(myOut)==TRUE){browser()}
  return(as.numeric(myOut))
}