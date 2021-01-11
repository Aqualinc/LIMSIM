###################################################################################
###################################################################################

#                   FLOW AND ABSTRACTION FUNCTIONS                                #

###################################################################################
###################################################################################
# Functions included in this script:
#       - Run_FlowandAbstraction   Controlling function that adds all flow and abstraction variables to MyREC
#       - GenFDC        Generates flow duraction curve
#       - TakeRate      Determines the abstraction for a given reach 
#       - freq.restrict Calculates the frequency of restrictions and the reliability
#
#
###################################################################################
###################################################################################

###################################################################
#' CONTROLLING FUNCTION THAT CREATES ALL OF THE FLOW AND ABSTRACTION VARIABLES IN MY REC
#'
#' This function calculates how much water should be taken from a reach based on which is the smaller of the Allocation and the maximum demand 
#' @param MyREC The REC data table with all the necessary attributes
#' @param AllocQ the allocation flow as a proportion of 7 day mean annual low flow
#' @param MinQ the minimum flow requirement as a proportion of the 7 dy mean annual low flow
#' @param GWAlloc the ground water allocation as a fraction of the mean annual aquifer recharge
#' @param pick a vector of the NZReach numbers of the river reaches of interest
#' @param TakeAll Set to 1 if all allocation is to be taken regardless of demand. Defaults to 0
#' @param SysCap the maximum amount of irrigation that could be applied in litres/second/hectare. Default to 0.58, equivalent of 5 mm per day.
#' @param EffIrriArea the fraction of irrigable land that is actually able to be irrigated (e.g. accounting for roads and houses etc). Defaults to 0.8
#' @return a dataframe of the MyREC parameter with a bunch of additional columns
#' @keywords REC
#' @export
#' @examples
#' Run_FlowandAbstraction()
Run_FlowandAbstraction<-function(MyREC=MyREC,AllocQ=AllocQ,MinQ=MinQ,GWAlloc=GWAlloc,pick=pick,TakeAll=0,SysCap=0.58,EffIrriArea=0.8){
  browser()
  load(file.path(dir$fun,"DeltaInd_v2.RData"))
  if(!exists("Doug.rbind.list", mode="function")) source(file.path(dir$fun,"GeneralFunctions.R"))
  
  #add the water take rules as attributes to the local "MyREC" dataframe (which is eventually returned as the function output)
  if (length(AllocQ==1)){  
    MyREC$MinQ<-matrix(MinQ,length(pick),1); MyREC$AllocQ<-matrix(AllocQ,length(pick),1); MyREC$GWAlloc<-matrix(GWAlloc,length(pick),1)
  }else{
    MyREC$MinQ<-MinQ; MyREC$AllocQ<-AllocQ; MyREC$GWAlloq<-GWAlloc
  }  
  
  # generate FDC for all NZReach.
  Perc <<- seq(0,1, length = 101) # Perc sets percentiles for which flows are to be generated at - needs to be set globally (hence the "<<-") for other functions (this is the default!)
  print("Calculating flow duration curves....")
  MyFDC <<- Doug.rbind.list(lapply(pick, GenFDC, Data = MyREC)) # generate FDC for all locations - also needs to be set globally (hence the "<<-")
  names(MyFDC) <- paste("P", 100*Perc, sep="") 
  
  print("Calculating irrigation takes and reliability...")
  #Save some of the original reach attributes  
  MyREC$BFIorig<-MyREC$BFI;MyREC$FRE2orig<-MyREC$FRE2;MyREC$FRE3orig<-MyREC$FRE3;MyREC$nNegorig<-MyREC$nNeg
  #Initialise some new attributes to NA
  MyREC$TheTake <-0; MyREC$IrriAreaIrrigated<-NA; MyREC$freq.man <-NA; MyREC$freq.min <-NA; MyREC$freq.reliability <-NA;MyREC$AllocQUse<-NA
   
  ind1<-which(MyREC$AllocQ>0)  #Determine those reaches where there is some allocation
  if (length(ind1)>0){
    MyRECtemp           <- MyREC[ind1,]                                               #Find the subset of the reaches that have an allocation greater than 0
    pick2               <- pick[ind1]                                                 #Get the NZReach number of the reaches that have an allocation greater than 0
    MyRECtemp$TheTake   <- sapply(pick2, TakeRate, RegAllocation=NULL, SysCap = SysCap, EffIrriArea = EffIrriArea, Data = MyRECtemp,TakeAll=TakeAll)  
    MyRECtemp$AllocQUse <-MyRECtemp$TheTake/(MyRECtemp$MALF*MyRECtemp$AllocQ)                             #TheTake / allocation
    MyRECtemp$AllocQUse[MyRECtemp$TheTake==0]<-NA
    MyRECtemp$IrriAreaIrrigated<-NA
    ind                 <-which(MyRECtemp$TheTake>0)
    
    MyRECtemp$IrriAreaIrrigated[ind]<-MyRECtemp$TheTake[ind]*1000*10000/SysCap/(MyRECtemp[ind,"usIrriArea"]*EffIrriArea)*100
    
    Restriction2 <- Doug.rbind.list(lapply(pick2, freq.restrict, FDC = MyFDC, prop=NULL, RegAllocation=NULL,Plot=F, Data = MyRECtemp,Data1=MyREC)) # turn Plot=F for speed.
    MyRECtemp$freq.man <- Restriction2$freq.man; MyRECtemp$freq.min <- Restriction2$freq.min; MyRECtemp$freq.reliability <- Restriction2$freq.reliability
    MyRECtemp$freq.man[MyRECtemp$TheTake==0]<-NA; MyRECtemp$freq.min[MyRECtemp$TheTake==0]<-NA; MyRECtemp$freq.reliability[MyRECtemp$TheTake==0]<-NA;
    
    MyRECtemp$BFI<-MyRECtemp$BFI+predict(dBFImod,newdata=MyRECtemp)                     #Update the BFI based on the dBFImod model
    MyRECtemp$BFI[MyRECtemp$TheTake==0]<-MyRECtemp$BFIorig[MyRECtemp$TheTake==0]        #Reset all the zero Take reaches back to 0
    
    MyRECtemp$FRE2<-MyRECtemp$FRE2+predict(dFRE2mod,newdata=MyRECtemp)                  #Update the FRE2 based on the dFRE2mod model
    MyRECtemp$FRE2[MyRECtemp$TheTake==0]<-MyRECtemp$FRE2orig[MyRECtemp$TheTake==0]      #Reset all the zero Take reaches back to 0
    
    MyRECtemp$FRE3<-MyRECtemp$FRE3+predict(dFRE3mod,newdata=MyRECtemp)                  #Update the FRE3 based on the dFRE3mod model
    MyRECtemp$FRE3[MyRECtemp$TheTake==0]<-MyRECtemp$FRE3orig[MyRECtemp$TheTake==0]      #Reset all the zero Take reaches back to 0
    
    MyRECtemp$nNeg<-MyRECtemp$nNeg+predict(dnNegMOD,newdata=MyRECtemp)                  #Update the nNeg based on the dnNegmod model
    MyRECtemp$nNeg[MyRECtemp$TheTake==0]<-MyRECtemp$nNegorig[MyRECtemp$TheTake==0]      #Reset all the zero Take reaches back to 0
    
    MyREC[ind1,]<-MyRECtemp
  }
  #browser()
return(MyREC)
}

###################################################################
#' function for generating FDC (m3/s) from GEV distribution (from R package nsRFA)
#'
#' This function calculates the exceedance flows for the given percentiles, based on the flow duration curve parameters of the reach 
#' @param ThisNZReach the NZReach number of the river reach on which to apply the function
#' @param P a vector of percentiles
#' @param Data the REC attribute table
#' @return a vector of flow rates of the same length as the P parameter
#' @keywords REC
#' @export
#' @examples
#' GenFDC()
GenFDC <- function(ThisNZReach = 13524724, P = Perc, Data = MyREC) {
  if((length(P))>1) {P[length(P)] <- 0.995}                       # re-set highest probability (to avoid generating infinite Q) (Note highest Q is extremely sensitive to this value)
  ThisRow <- which(Data$NZReach == ThisNZReach)                   # get the parameters from the FDC dataset is this faster???
  xi      <- Data[ThisRow, "XI"]                                  # the parameters  xi, alfa, k will be found in MyFDCPars
  alfa    <- Data[ThisRow, "ALFA"]
  k       <- Data[ThisRow, "K"]
  
  #browser()
  if (is.na(k)==F){
    if ((k > -1e-07) & (k < 1e-07)) {
      Q <- xi - alfa * log(-log(P))
    } else {
      Q <- xi + alfa * (1 - (-log(P))^k)/k
    }
    Q[Q<0] <- 0                                                     # no values should be negative
    Q <- Q * Data[ThisRow, "MeanFlow"]                              # multiply by estimated mean flow
    #about here I need to subtract from Q the GWAlloc * baseflowfrom Aquifers or do I?
    #Q  <- Q - Data[ThisRow, "GWAlloc"] * Data[ThisRow, "AqBaseFlow"] * Data[ThisRow, "BaseFlow"]
  } else {
    Q<-P*NA
  }
  return(Q)
}


###############################################################################
#' A function to compute the amount of take from any location
#'
#' This function calculates how much water should be taken from a reach based on which is the smaller of the Allocation and the maximum demand 
#' @param ThisNZReach the NZReach number of the river reach on which to apply the function
#' @param RegAllocation the allocation as a fraction of MALF. If NULL then it is set to the "AllocQ" attribute of the reach. Defaults to 0.5
#' @param Data the REC attribute table
#' @param SysCap the maximum amount of irrigation that could be applied in litres/second/hectare. Default to 0.58, equivalent of 5 mm per day.
#' @param TakeAll Set to 1 if all allocation is to be taken regardless of demand.
#' @param EffIrriArea the fraction of irrigable land that is actually able to be irrigated (e.g. accounting for roads and houses etc). Defaults to 0.8
#' @return one number that is the take.
#' @keywords REC
#' @export
#' @examples
#' TakeRate()
TakeRate <- function(ThisNZReach = 13524724, RegAllocation=0.5, Data = 0, SysCap = 0.58,TakeAll=TakeAll,EffIrriArea=0.8) {  
  
  ThisRow <- which(Data$NZReach == ThisNZReach) # the row in the dataset for this NZREach
  #browser()
  if (is.null(RegAllocation)){RegAllocation<-Data$AllocQ[ThisRow]}                  
  if(Data[ThisRow,"usIrriArea"] > 10000 & !is.na(Data[ThisRow,"usIrriArea"])) {            #Check if there is more than a hectare upstream
    MaxDemand     <- Data[ThisRow, "usIrriArea"] / 10000 * SysCap / 1000 * EffIrriArea     #The upstream irrigable area (converted from m2 to hectares) multiplied by the use rate (converted from l/s to m3/s) multiplied by an efficiency factor     
    SWTake <- RegAllocation * Data[ThisRow, "MALF"]                                        #The maximum surface water take in cumecs
    GWTake <- Data[ThisRow, "AqBaseFlow"] * Data[ThisRow,"TrueBFI"] * Data[ThisRow, "MeanFlow"] * Data[ThisRow, "GWAlloc"] #The maximum groundwater take
    RegulatedTake <- SWTake + GWTake
    Take          <- min(MaxDemand, RegulatedTake)
    if (TakeAll==1) {Take<-RegulatedTake}
  } else {
    Take <- 0                                                                       #If the up stream irrigable area is less than a hectare then we assume the landscape is too steep to irrigate No water is taken
  }
  return(Take)}


#######################################################################################
#' A function to calculate the frequency of restrictions
#'
#' This function calculates the amount of time that the allocation will be restricted because of lack of water
#' @param ThisNZReach the NZReach number of the river reach on which to apply the function
#' @param FDC The flow duration curve
#' @param minFlow the minimum flow in cumecs that is allowed. NULL means it will be calculated from "prop". Defaults to NULL
#' @param prop the minimum flow in fraction of the "Qref". prop stands for proportion. Defaults to 0.8
#' @param Qref the reference flow in cumecs which the "prop" and "RegAllocation" are fractions of. If it is NULL, then it is set to the "MALF" attribute of the reach. Defaults to NULL
#' @param RegAllocation the allocation as a fraction of "Qref". If NULL then it is set to the "AllocQ" attribute of the reach. Defaults to 0.5
#' @param Plot whether to plot the flow duration curve with lines showing the minimum and managed flows in normal and log space. Defaults to TRUE
#' @param Data the REC attribute table which includes "minQ", "TheTake" and "AllocQ" attributes
#' @param Data1 the REC attribute table which includes the flow duration curve attributes
#' @return two values: 1/ The percentage of time restritions apply 2/ The percentage of time of TOTAL restriction
#' @keywords REC
#' @export
#' @examples
#' freq.restrict()
freq.restrict <- function(ThisNZReach = 13524724, FDC = MyFDC,  minFlow=NULL, prop=0.8, Qref=NULL, RegAllocation=0.5, Plot = TRUE, Data = MyREC,Data1=MyREC) {
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach)                          # get the parameters from the FDC dataset is this faster???
  ThisRow_ALL <- which(Data1$NZReach == ThisNZReach)
  if (is.null(prop)){prop<-Data$MinQ[ThisRow]}
  if (is.null(RegAllocation)){RegAllocation<-Data$AllocQ[ThisRow]}
  if (is.null(Qref)) Qref <- Data[ThisRow, "MALF"]                       # the reference flow is the MALf from the hyd predictions data
  if (is.null(minFlow)) minFlow <- Qref * prop                           # minflow = flow at which there is TOTAL restriction
  #if(Data$TheTake[ThisRow] > 0) 
  #browser()
  
  if(Data$TheTake[ThisRow] == 0 | is.na(Data$TheTake[ThisRow]) | any(is.na(FDC[ThisRow, ]))) {        # if the take is 0 or unknown or if the FDC is unknown then the restrictions are unknown
    freq.man <- NA
    freq.min <- NA
    freq.reliability<-NA
  }else{
    #browser()
    if(!is.null(Data$TheTake[ThisRow]))   {
      allocation <-  Data$TheTake[ThisRow]                                        # if a take for each reach is defined then use it
    } else {
      allocation <- RegAllocation*Qref                                            #  OTHERWISE a regulatory rule of thumb allocation is used Qref*allocate. The DEFAULT for Qref= MALF, (NES type rule)
    }
    
    if(allocation == 0) minFlow <- 0                                              # if there are no takes there should be no restrictions (min=0)                     
    manFlow          <- minFlow+allocation                                        #  minflow+allocation = "management flow"  the flow at which restrictions begin
    FDC.data         <- FDC[ThisRow_ALL, ]                                        # THIS is the FDC
    freqs            <- 100*Perc                                                    # the percentiles the FDC is estimated for from the global value of Perc              
    
    #Adjust the flow duration curve for any groundwater allocation
    GWTakeFDCShift <- Data$GWAlloc[ThisRow] * Data$AqBaseFlow[ThisRow] * Data$BaseFlow[ThisRow] #Calculate the shift in the FDC resulting from the Groundwater allocation
    #GWTakeFDCShift <- 0 
    FDCGWTakeAffected.data <- FDC.data - GWTakeFDCShift
    #Set any -ves to 0
    FDCGWTakeAffected.data[FDCGWTakeAffected.data <0] <- 0
    
    #browser()
    if (minFlow < FDCGWTakeAffected.data[1]){minFlow <- FDCGWTakeAffected.data[1]}              # Make sure the minimum flow is at least as large as the smallest value in the flow duration curve
    
    freq.diff <-  approx(y=freqs, x=FDCGWTakeAffected.data, xout=c(minFlow+allocation, minFlow))$y # estimate the percentiles for which the "restricted takes" and "stopped takes" occur
    freq.man         <- freq.diff[1]
    freq.min         <- freq.diff[2]                                              # minflow= flow at which there is TOTAL restriction   
    
    #if (length(freq.diff)!=2) {browser()}
    
    
    #AlteredFDC<-FDCGWTakeAffected.data - approx(x=c(100,freq.diff,0),y=c(allocation+GWTakeFDCShift,allocation+GWTakeFDCShift,GWTakeFDCShift,GWTakeFDCShift),xout=freqs)$y
    
    #browser()
    managed.freqs <- rev(subset(freqs, freqs < freq.diff[1] & freqs > freq.diff[2]))  #These are the percentiles that disappear, in reverse order
    managed.freq.indices <- rev(which(freqs %in% managed.freqs))                      #These are the indices of the percentiles that disappear, in reverse order             
    
    #The new flow duration curve has flows reduced by the allocation flow from percentiles 100 to the percentile of the allocation flow + minimum flow (i.e. freq.diff[1])
    #The new flow duration curve is unchanged below the minimum flow, i.e. from the percentile of the miniumum flow (freq.diff[2]) down to 0
    #Between the upper and lower managed flow percentiles, the flows all go to the minimum flow.
    #Note that in a previous version, the flow reduction in the managed flow section was a linear interpoaltion between allocation and 0. This results
    #in a strange flow duration curve when the increase in flows between percentiles is less than the interpolated change. You can end up with a hollow in the flow duration curve.
    if(length(managed.freqs > 0)) managed.flows <- FDCGWTakeAffected.data[managed.freq.indices]-minFlow else managed.flows <- c()   #If there are no percentiles between the minimum and managed percentiles, then set the related flows to an empty set
    FDC.take<-approx(x=c(100,freq.diff[1],managed.freqs,freq.diff[2],0),y=c(allocation,allocation,managed.flows,0,0),xout=freqs)
    
    #browser()
    
    #FDC.take         <- approx(x=c(100,freq.diff,0),y=c(allocation,allocation,0,0),xout=freqs)
    freq.reliability <- mean(FDC.take$y)/allocation*100  #NOTE THIS ONLY WORKS IF FREQS iS EVENLY SAMPLED
    #if(is.na(freq.min)) freq.min <- 100  # the minFlow can be OFF the end of the FDC curve! (especially when MALF is annual but FDC is for summer.
    #if(is.na(freq.man)) freq.man <- 100  # the FREQ manFlow can be NA in MALF = 0                                        
  }
    
  if(Plot==TRUE)   {
    par(mfrow=c(2,1), bg="grey90")
    #plot(freqs, log10(FDC.data), col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="log10 Flow")
    plot(freqs, log10(FDCGWTakeAffected.data), col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="log10 Flow")
    abline(h=log10(minFlow), col = "green")
    abline(v=freq.min, col = "green", lty=2)
    abline(h=log10(manFlow), col = "red")
    abline(v=freq.man, col = "red", lty=2)
    legend("topleft", cex=0.8,
           legend=c(paste("Minimum Flow",round(minFlow,2)),
                    paste("Freq Min Flow", round(freq.min, 2)),
                    paste("Management Flow", round(manFlow,2)),
                    paste("Freq Man Flow",round(freq.man,2))),
           col = c("green","green", "red", "red"),   text.col = "black", lty = c(1, 2,1,2), pch = c(-1, -1),   lwd = c(1,1,1,1),
           merge = TRUE, bg = 'white', inset = .05)
    #plot(freqs, FDC.data, col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="Flow",
    #     xlim=c(freq.man*0.9, 100), ylim=c(0, manFlow*1.1))
    plot(freqs, FDCGWTakeAffected.data, col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="Flow",
         xlim=c(freq.man*0.9, 100), ylim=c(0, manFlow*1.1))
    abline(h=minFlow, col = "green")
    abline(v=freq.min, col = "green", lty=2)
    abline(h=manFlow, col = "red")
    abline(v=freq.man, col = "red", lty=2)
  } #end plot
  
  return(cbind(freq.man, freq.min,freq.reliability))   #returns two values 1 = % of time restrition 2 % of time TOTAL restriction (note the variable names (freq) is poor use of words these are not frequencies!
}# end

#######################################################################################
#' A function to calculate the reliability of supply for each band of multiband restrictions
#'
#' This function calculates the percentage of maximum allocation that is available for each band, assuming lower bands are fully allocated
#' @param ThisNZReach the NZReach number of the river reach on which to apply the function
#' @param FDC The flow duration curve
#' @param minFlow a vector of minimum flows in cumecs that is allowed
#' @param allocation a vector of allocation flows in cumecs that is allowed
#' @param Data the REC attribute table which includes "minQ", "TheTake" and "AllocQ" attributes
#' @param Data1 the REC attribute table which includes the flow duration curve attributes
#' @return a three column dataframe of 1/ The percentage of time restritions apply 2/ The percentage of time of TOTAL restriction 3/ The fraction of full allocation
#' @keywords REC
#' @export
#' @examples
#' freq.restrict()
freq.restrict.multiband <- function(ThisNZReach = 13524724, FDC = MyFDC,  minFlow=NULL, prop=0.8, Qref=NULL, RegAllocation=0.5, Plot = TRUE, Data = MyREC,Data1=MyREC) {
  ThisRow <- which(Data$NZReach == ThisNZReach)                          # get the parameters from the FDC dataset is this faster???
  ThisRow_ALL <- which(Data1$NZReach == ThisNZReach)
  FDC.data         <- FDC[ThisRow_ALL, ]                                        # THIS is the FDC
  
  
  freqs            <- 100*Perc                                                    # the percentiles the FDC is estimated for from the global value of Perc
    
  #Adjust the flow duration curve for any groundwater allocation
    GWTakeFDCShift <- Data$GWAlloc[ThisRow] * Data$AqBaseFlow[ThisRow] * Data$BaseFlow[ThisRow] #Calculate the shift in the FDC resulting from the Groundwater allocation
    #GWTakeFDCShift <- 0 
    FDCGWTakeAffected.data <- FDC.data - GWTakeFDCShift
    #Set any -ves to 0
    FDCGWTakeAffected.data[FDCGWTakeAffected.data <0] <- 0
    
    if (minFlow < FDCGWTakeAffected.data[1]){minFlow <- FDCGWTakeAffected.data[1]}              # Make sure the minimum flow is at least as large as the smallest value in the flow duration curve
    
  #Need to loop through each band
  for (Band in seq(1,length(minFlow))) {
    
  
    if(allocation[band] == 0) minFlow[band] <- 0                                              # if there are no takes there should be no restrictions (min=0)                     
    manFlow          <- minFlow[band]+allocation[band]                                        #  minflow+allocation = "management flow"  the flow at which restrictions begin
    
    freq.diff <-  approx(y=freqs, x=FDCGWTakeAffected.data, xout=c(minFlow[band]+allocation[band], minFlow[band]))$y # estimate the percentiles for which the "restricted takes" and "stopped takes" occur
    freq.man         <- freq.diff[1]
    freq.min         <- freq.diff[2]                                              # minflow= flow at which there is TOTAL restriction   
    
    managed.freqs <- rev(subset(freqs, freqs < freq.man & freqs > freq.min))     #These are the percentiles that disappear, in reverse order
    managed.freq.indices <- rev(which(freqs %in% managed.freqs))                      #These are the indices of the percentiles that disappear, in reverse order             
    
    #The new flow duration curve has flows reduced by the allocation flow from percentiles 100 to the percentile of the allocation flow + minimum flow (i.e. freq.diff[1])
    #The new flow duration curve is unchanged below the minimum flow, i.e. from the percentile of the miniumum flow (freq.diff[2]) down to 0
    #Between the upper and lower managed flow percentiles, the flows all go to the minimum flow.
    #Note that in a previous version, the flow reduction in the managed flow section was a linear interpoaltion between allocation and 0. This results
    #in a strange flow duration curve when the increase in flows between percentiles is less than the interpolated change. You can end up with a hollow in the flow duration curve.
    if(length(managed.freqs > 0)) managed.flows <- FDCGWTakeAffected.data[managed.freq.indices]-minFlow[band] else managed.flows <- c()   #If there are no percentiles between the minimum and managed percentiles, then set the related flows to an empty set
    FDC.take<-approx(x=c(100,freq.man,managed.freqs,freq.min,0),y=c(allocation,allocation[band],managed.flows,0,0),xout=freqs)
    freq.reliability[band] <- mean(FDC.take$y)/allocation[band]*100  #NOTE THIS ONLY WORKS IF FREQS iS EVENLY SAMPLED
    
    #Now need to update the FDC to account for the takes
    #Original FDC from 0 to minimum flow, then step down to % at managed flows, then retain offset to maximum flows
    #FDCGWTakeAffected.data <- 
    #Need to figure out how to handle % allocation share. This means that between onset of restrictions and full restrictions, the available allocation might be made fully available for allocation, or maybe just 50 % is available for allocation (or some other percentage). This avoids flat lining the stream.
}                              
 

  return(cbind(freq.man, freq.min,freq.reliability))   #returns three columns freq.man = % of time on restrition freq.min = % of time TOTAL restriction, freq.reliability = % of max allocation available (note the variable names (freq) is poor use of words these are not frequencies!
}# end