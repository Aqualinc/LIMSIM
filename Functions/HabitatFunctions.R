###################################################################################
###################################################################################

#                         HABITAT FUNCTIONS                                       #

###################################################################################
###################################################################################
# Functions indcluded in this script:
#       - HabitatFunctions runs all the other functions for the reaches of interest within an REC set.
#       - IntegratedWidth Calculates the reduction in width over the whole FDC
#       - habitat         Calculates the available habitat for a prescribed species
#       - delta.hab       Calcualtes the change in habitat


###################################################################################
###################################################################################

#need to make sure that the functions within the Flow&AbstractionFunctions.R are available.
if(!exists("Run_FlowandAbstraction", mode="function")) source(file.path(dir$fun,"FlowAbstractionFunctions.R"))


###################################################################
#' Function to run all the habitat functions
#'
#' @description this function runs all the LimSim habitat functions
#' @references Booker, D. j., 2010. Predicting wetted width in any river at any discharge. Earth Surf. Process. Landforms 35, 828–841. doi:10.1002/esp.1981
#' @references Jowett, I.G., 1998. Hydraulic geometry of New Zealand rivers and its use as a preliminary method of habitat assessment. Regul. Rivers: Res. Mgmt. 14, 451–466. doi:10.1002/(SICI)1099-1646(1998090)14:5<451::AID-RRR512>3.0.CO;2-1
#' @param pick The REC reach numbers of interest
#' @param method either "Booker" based on Booker 2009 (default) or "Jowett" based on Jowett 1998
#' @param MyREC the REC attribute table
#' @param MinQ the minimum flow rule as a fraction of the seven day MALF
#' @return A dataframe of flow (Q) and width (W)
#' @keywords REC
#' @export
#' @examples
#' Run_HabitatFunctions()

Run_HabitatFunctions<-function(MyREC=NULL,pick=pick,MinQ=0,method="Jowett"){

  if (method=="Booker"){
    print("Calculating Booker QW relationships...")
    load(paste(dir$fun,"BookerModel.RData",sep="/"))  # this loads Doug's fitted models to estimate width using catchment area and climate category 
    QW <<- lapply(pick, GetWidth, method="Booker",Data=MyREC,IntModelCl=IntModelCl,SlopeModelCl=SlopeModelCl,QuadModelCl=SlopeModelCl) # this gets the width flow relationships for all reaches
  }else{
    print("Calculating Jowett QW relationships...")
    QW <<- lapply(pick, GetWidth, method="Jowett",Data=MyREC) 
  }
#QW2 <<- lapply(pick, GetWidth, method="booker") # this gets the width flow relationships for all reaches.
names(QW) <- pick   # names the list elements
# the mean reduction in width for flows up to the mean (the largest value in QW
print("Calculating reduction in width...")
#browser()
MyREC$RedWidth<- sapply(pick, IntegratedWidth, FDC = MyFDC, FlowWidth = QW, minFlow=NULL,
                          Qref=NULL, Plot = F, Data = MyREC)
#RedWidth(is.na(RedWidth))<-0
#MyREC$RedWidth[which(is.na(MyREC$RedWidth)==T)]<-0

MyREC$NaturalWidth <- sapply(pick,NaturalWidth,FDC = MyFDC, FlowWidth = QW, Data = MyREC)

print("Calculating change in habitat...")
#browser()
QWH <- lapply(pick, habitat,sp="LongEel",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaEEL<-delta.EEL

QWH <- lapply(pick, habitat,sp="ShortEel",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaShortEEL<-delta.EEL

QWH <- lapply(pick, habitat,sp="Brown trout adult",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaTrout<-delta.EEL

QWH <- lapply(pick, habitat,sp="Bluegill Bully",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaBully<-delta.EEL

QWH <- lapply(pick, habitat,sp="Inanga",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaInanga<-delta.EEL

QWH <- lapply(pick, habitat,sp="Torrent",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaTorrent<-delta.EEL

QWH <- lapply(pick, habitat,sp="Kokopu",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaKokopu<-delta.EEL

QWH <- lapply(pick, habitat,sp="ComBully",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaComBully<-delta.EEL

QWH <- lapply(pick, habitat,sp="Fry",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaTroutFry<-delta.EEL

QWH <- lapply(pick, habitat,sp="Spawn",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaTroutSpawn<-delta.EEL

QWH <- lapply(pick, habitat,sp="Upland Bully",FDC = MyFDC,FlowWidth = QW,Data=MyREC)  #  , sp="Fry"
names(QWH) <- pick  # names the list elements
delta.EEL <- sapply(pick, delta.hab, Plot = F,prop=MinQ,FlowHab = QWH, Hab = "WUA",Data = MyREC)
MyREC$deltaUplandBully<-delta.EEL

return(MyREC)
}


###################################################################
#' Function to compute the flow to width relationship for a reach
#'
#' @description this function prepares a two column dataframe of flow (Q) and width (W) for flows
#'  from 0.0001 of the mean flow through to the mean flow in a  specified number of steps using one of two methods
#' @references Booker, D. j., 2010. Predicting wetted width in any river at any discharge. Earth Surf. Process. Landforms 35, 828–841. doi:10.1002/esp.1981
#' @references Jowett, I.G., 1998. Hydraulic geometry of New Zealand rivers and its use as a preliminary method of habitat assessment. Regul. Rivers: Res. Mgmt. 14, 451–466. doi:10.1002/(SICI)1099-1646(1998090)14:5<451::AID-RRR512>3.0.CO;2-1
#' @param ThisNZReach The REC reach number of the reach of interest
#' @param method either "Booker" based on Booker 2009 (default) or "Jowett" based on Jowett 1998
#' @param n.pairs The number of flow divisions to provide widths for
#' @param Data the REC attribute table
#' @param IntModelCl a model to calculate intercept of the flow to width relationship, from Booker 2010, with R model in the BookerModel.RData found in the "Functions" sub directory
#' @param SlopeModelCl a model to calculate the slope of the flow to width relationship, from Booker 2010, with R model in the BookerModel.RData found in the "Functions" sub directory
#' @param QuadModelCl a model to calculate quadratics coefficient of the flow to width relationship, from Booker 2010, with R model in the BookerModel.RData found in the "Functions" sub directory
#' @return A dataframe of flow (Q) and width (W)
#' @keywords REC
#' @export
#' @examples
#' GetWidth()

GetWidth <- function(ThisNZReach = 13524724, method="Booker", n.pairs = 10, Data=MyREC,IntModelCl=NA,SlopeModelCl=NA,QuadModelCl=NA) {   
  ThisRow <- which(Data$NZReach == ThisNZReach) # find which row of the REC data has the reach of interest
  this.seg <- Data[ThisRow, ]                   # get all the data for the reach of interest
  #browser()
  #Qbar <- this.seg$Flow_L_s/1000
  Qbar <- Data[ThisRow, "MeanFlow"]             # get the mean flow for the reach of interest
  if(is.na(Qbar)==F){                           # Only do the calculation if there is a value for the mean flow 
    Q <- seq(Qbar*0.0001, Qbar, length.out=n.pairs)  # prepare a vector of flows from 0.0001 through to the mean flow in n.pairs steps
    if (method=="Booker") {                     # Check if the "Booker" method is to be used
      this.seg$LogCatchmentArea <- log10(this.seg$usArea/1000000)    #Find the log10 of the upstream area in metres squared
      this.seg$LogFlow <- log10(Qbar)                                #Find the log10 of the mean flow
      d0<- predict(IntModelCl, newdata = this.seg)                   #Find three parameters, d0,d1,d2 using three models  
      d1 <- predict(SlopeModelCl, newdata = this.seg)
      d2 <- predict(QuadModelCl, newdata = this.seg)
      LogQ <- log10(Q)
      W  <-  10^(d0 +  LogQ*d1  +  d2*(LogQ^2))                      #Calculate the width 
      QW <- data.frame(cbind(Q, W))
    } else {                              # The other option ids the "Jowett" method
      Wbar <- 7.76*(Qbar^0.488)           #  mean width and mean flow from Downstream Hydraulic Geometry  Jowett 1998
      b    <- Wbar/(Qbar^0.176)           # Compute b far at-a-station hydraulic geometry (Width) Jowett 1998
      W    <- b*(Q^0.176)
      QW   <- data.frame(cbind(Q, W))
    }
  }else{
    QW<-NA
  }
  return(QW)
  #browser()
}


###################################################################
#' Function to compute loss of width over the whole hydrograph (FDC)
#'
#' this function evaluates the reduction in width over the whole FDC for the natural and altered flows 
#' @param ThisNZReach The REC reach number of the reach of interest
#' @param FDC a dataframe of flow rates for different percentiles (columns) for different reaches (rows)
#' @param FlowWidth dataframe of flow vs width
#' @param minFlow the minimum flow in cumecs. If this is NULL then it is calculated from the MinQ attribute in the REC table. Default is NULL
#' @param Qref the reference flow in cumecs which the "prop" and "allocate" are fractions of. If it is NULL, then it is set to the "MALF" attribute of the reach. Defaults to NULL
#' @param Plot whether to plot the flow duration curve with lines showing the minimum and managed flows in normal and log space. Defaults to TRUE
#' @param Data the REC attribute table
#' @return The average width loss resulting from the allocations
#' @keywords REC
#' @export
#' @examples
#' IntegratedWidth()

IntegratedWidth <- function(ThisNZReach = 13524724, FDC = MyFDC, FlowWidth = QW, minFlow=NULL, Qref=NULL, Plot = TRUE, Data = MyREC) {
  
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach)                        # Get the row number of the REC attribute table for the reach of interest 
  if(Data$TheTake[ThisRow] == 0 | anyNA(FlowWidth[[ThisRow]])) {  # if there is no take or there is an NA in the FlowWidth table, skip entirely width loss is ZERO  -> NA
    AveWidthLoss <- NA
  } else {
    
    #Sort out the values for the minimum flow and allocation flow in cumeces
    if (is.null(Qref)) Qref <- Data$MALF[ThisRow]                      # If the "Qref" parameter is NULL, set it to the MALf from the REC attributes table
    if (is.null(minFlow)) minFlow <- Qref * Data$MinQ[ThisRow]         # If the "minflow" parameter is NULL, set it based on the REC attribute tables "MinQ" and "MALF"
    
    if(!is.null(Data$TheTake))   {                                     # Check to see if take has been established
      allocation <-  Data$TheTake[ThisRow]                             # If it has, then use it to set the allocation
    } else {
      allocation <- Data$AllocQ[ThisRow] * Qref                        # if not then the allocation will be Qref * allocate. The DEFAULT for Qref= MALF, this allows for rules of thumb to be used.
    }
    
    if (anyNA(FDC[ThisRow, ])){
      AveWidthLoss <-NA
    }else{
      # obtain FDC data for this REACH 
      FDC.data <- FDC[ThisRow, ]                                      # THIS is the FDC
      freqs    <- 100*Perc                                            # the percentiles the FDC is estimated for from the global value of Perc
      
      #Adjust the flow duration curve for any groundwater allocation
      GWTakeFDCShift <- Data$GWAlloc[ThisRow] * Data$AqBaseFlow[ThisRow] * Data$BaseFlow[ThisRow] #Calculate the shift in the FDC resulting from the Groundwater allocation
      #GWTakeFDCShift <- 0 
      FDCGWTakeAffected.data <- FDC.data - GWTakeFDCShift
      #Set any -ves to 0
      FDCGWTakeAffected.data[FDCGWTakeAffected.data <0] <- 0
      
      #browser()
      if (minFlow < FDCGWTakeAffected.data[1]){minFlow <- FDCGWTakeAffected.data[1]}              # Make sure the minimum flow is at least as large as the smallest value in the flow duration curve
      
      freq.diff <-  approx(y=freqs, x=FDCGWTakeAffected.data, xout=c(minFlow+allocation, minFlow))$y # estimate the percentiles for which the "restricted takes" and "stopped takes" occur
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
      AlteredFDC<-FDCGWTakeAffected.data - approx(x=c(100,freq.diff[1],managed.freqs,freq.diff[2],0),y=c(allocation,allocation,managed.flows,0,0),xout=freqs)$y
     
      #browser()
      Qa <- FlowWidth[[ThisRow]]$Q                                    # the QW calculations corresponding to the reach of interest
      W  <- FlowWidth[[ThisRow]]$W
        
      NaturalWidths <- approx(x=c(0,Qa), y=c(0,W), xout=FDC.data)$y   # interploate W vs Q data
      AlteredWidths <- approx(x=c(0,Qa), y=c(0,W), xout=AlteredFDC)$y # interploate W vs Q data
      
      AveWidthLoss <-   mean((AlteredWidths - NaturalWidths)/NaturalWidths, na.rm=T) *100 # mean reduction in width (percentage of natural flow) 
    }
    #browser()
    if(Plot==TRUE)   {
      x11(); par(mfrow=c(2,1), bg="grey90")
      freq.diff <-  approx(y=freqs, x=FDC.data, xout=c(minFlow+allocation, minFlow))$y # interploate W vs Q data
      freq.man <- freq.diff[1]
      freq.min <- freq.diff[2]     # minflow= flow at which there is TOTAL restriction
      plot(freqs, FDC.data, log="y", col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="Flow")
      points(freqs, AlteredFDC,col="red", lwd=2, type="l")
      abline(v=freq.min, col = "green", lty=2)
      abline(v=freq.man, col = "green", lty=2)
      plot(freqs, NaturalWidths, log="y", col="blue", lwd=2, type="l", xlab="Time flow is equalled or exceeded (%)", ylab="Width")
      points(freqs, AlteredWidths,col="red", lwd=2, type="l")
    } #end plot
  #browser()
  }  
  #browser()
  return(AveWidthLoss)   #returns the loss of width
}# end

NaturalWidth <- function(ThisNZReach = 13524724, FDC = MyFDC, FlowWidth = QW, Data = MyREC) {
  
  ThisRow <- which(Data$NZReach == ThisNZReach) # 
  if (anyNA(FDC[ThisRow, ])==T){
    NaturalWidths<-NA 
    
  }else{
    Qa <- FlowWidth[[ThisRow]]$Q   # the QW calculations corresponding to the NZReach
    W <- FlowWidth[[ThisRow]]$W
    NaturalWidths <- approx(x=c(0,Qa), y=c(0,W), xout=Data$MeanFlow[ThisRow])$y # interploate W vs Q data
  }
  return(NaturalWidths)
}
###############################################################################                                                                            #
#         Compute habitat vs flow for each NZreach                            #                                                                             #
###############################################################################

habitat <- function(ThisNZReach = 13524724, FlowWidth = QW, Data = MyREC, sp="Brown trout adult",FDC = MyFDC) {  #
  #browser()
  if(any(sp%in%"Brown trout adult")) hab.par <- c(1.17, 4.35)    # hab.par = parameters for generalised habitat models C and K
  if(any(sp%in%"Redfin"))            hab.par <- c(0.26, 7.39)    # redfin
  if(any(sp%in%"ComBully"))          hab.par <- c(0.39, 6.51)    # common bully
  if(any(sp%in%"Fry"))               hab.par <- c(0.86, 10.21)   # brown trout fry
  if(any(sp%in%"Spawn"))             hab.par <- c(1.24, 9.89)    # brown trout spawning
  if(any(sp%in%"Deli"))              hab.par <- c(0.33, 1.92)    # deleatidium
  if(any(sp%in%"LongEel"))           hab.par <- c(0.07, 2.07)    # longfin eel
  if(any(sp%in%"Torrent"))           hab.par <- c(0.88, 4.05)    # torrent fish
  if(any(sp%in%"ShortEel"))          hab.par <- c(0.13, 2.32)    # shortfin eel
  if(any(sp%in%"Bluegill Bully"))    hab.par <- c(1.01, 6.13)    # Bluegill bully
  if(any(sp%in%"Inanga"))            hab.par <- c(0.19, 19.74)   # INanaga
  if(any(sp%in%"Kokopu"))            hab.par <- c(0.03, 2.29)    # Kokopu
  if(any(sp%in%"Upland Bully"))      hab.par <- c(0.19, 13.13)   # upland bully
  
  ThisRow <- which(Data$NZReach == ThisNZReach) #
  if (anyNA(FDC[ThisRow, ])==T){
    Q<-NA; HV<-NA; WUA<-NA
  }else{
  Q <- FlowWidth[[ThisRow]]$Q
  W <- FlowWidth[[ThisRow]]$W
  HV.raw <- (Q/W)^hab.par[1] *  exp(-hab.par[2] * (Q/W))
  # HV.raw[1] <- 0   # replace the first column with zero NaN produced above
  HV <- HV.raw/max(HV.raw) # normalise the HV values
  WUA <- HV*W   # convert to WUA by x by width estimates
  }
  
  habMat <- data.frame(cbind(Q, HV, WUA))
  return(habMat)
}

###################################################################
#' Function to calculate the change in habitat for a change in flow for a reach
#'
#' this function evaluates the percentage change in habitat based on the change in flow 
#' @param ThisNZReach The REC reach number of the reach of interest
#' @param Data The REC attribute table
#' @param FlowHab A list for each reach giving a lookup table of flow, habitat value (HV) and weighted usable area (WUA)
#' @param Hab The type of habitat measure to return. Options are "WUA" or "HV". "WUA" is weighted usable area. The percentage of habitat area available at a particular flow.
#' "HV" is the habitat value. Value from 0 to 1 to represent what proportion of the maximum possible habitat value at a particular flow.
#' @param Qref the reference flow in cumecs which the "MinQ" and "AllocQ" are fractions of. If it is NULL, then it is set to the "MALF" attribute of the reach. Defaults to NULL
#' @param prop Unused
#' @param Plot Whether to plot the flow vs habitat curve with locations of original MALF and new MALF
#' @return The new habitat area as a percetage of the original.
#' @keywords REC
#' @export
#' @examples
#' delta.hab()

delta.hab <- function(ThisNZReach = 13524724, Data = MyREC, FlowHab = QWH, Hab = "WUA", Qref=NULL, prop=0.8, Plot = TRUE) {  # returns the widths at Qref and prop times Qref and plots this QW curve if plot=TRUE
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach)
  if (is.null(Qref)) { Qref <- Data[ThisRow, "MALF"] } # the reference flow is the MALf from the hyd predictions data
  QrefNew <- Data[ThisRow, "MinQ"]*Qref               #This is the new minimum flow after allocation
  #browser()
  
  #With groundwater takes it may be possible to end up with a minimum flow less than the minimum flow rule (e.g. if the reach is all aquifer sourced)
  #So we need to estimate what the MALF might be just from GW takes. Do this by finding the frequency of MALF under natural flows, then
  # find the flow at the same frequency from the flow duration curve calculated accounting for groundwater takes

  #Find the FDC for the reach of interest
  Perc <- seq(0,1, length = 101)
  FDC.data <- GenFDC(ThisNZReach = ThisNZReach, P = Perc, Data = Data)
  #Find frequency of Qref
  QrefFreq <- approx(x=FDC.data,y=Perc, xout=Qref)$y
  #Find the Altered FDC accounting for the Groundwater abstraction
  #Adjust the flow duration curve for any groundwater allocation
  GWTakeFDCShift <- Data$GWAlloc[ThisRow] * Data$AqBaseFlow[ThisRow] * Data$BaseFlow[ThisRow] #Calculate the shift in the FDC resulting from the Groundwater allocation
  #GWTakeFDCShift <- 0 
  FDCGWTakeAffected.data <- FDC.data - GWTakeFDCShift
  #Set any -ves to 0
  FDCGWTakeAffected.data[FDCGWTakeAffected.data <0] <- 0
  #Find the new flow for the Altered FDC
  QrefNewGW <- approx(x=Perc, y=FDCGWTakeAffected.data,xout=QrefFreq)$y
  #If this is less than the reaches minimum flow (as given by the product of its "MinQ" and  attribute)
  if (QrefNewGW < QrefNew) QrefNew <- QrefNewGW
  
  Q <- FlowHab[[ThisRow]]$Q   #  get flow vs habitat
  H <- unlist(FlowHab[[ThisRow]][Hab])

  
  
  if(Qref==0|is.na(Qref)==T| Data$TheTake[ThisRow]==0) {  # if Qref is zero the habitat is ZERO and so is delta-habitat
    hab.Qref <- NA
    hab.Qnew <- NA
    DeltaHab <- NA
  } else {
    
    #The next line is where props is supposed to be used, but instead it uses MinQ
    #QrefNew <- Data[ThisRow, "MinQ"]*Qref
    habs <-  approx(x=Q, y=H, xout=c(Qref, QrefNew))$y # interploate W vs Q data
    #habs <-  approx(x=Q, y=H, xout=c(Qref, Data[ThisRow, "MinQ"]*Qref))$y # interploate W vs Q data
    hab.Qref <- habs[1]
    hab.Qnew <- habs[2]
    DeltaHab <- hab.Qnew/hab.Qref*100
  } # end if
  
  if (Plot==TRUE) {
    par(bg="grey90")
    plot(Q, H, type="l", lwd = 2,
         main=paste("Difference in", Hab, "for", Data[ThisRow, "MinQ"], "times 7DayMALF\nand groundwater allocation of", Data[ThisRow, "GWAlloc"] ),
         xlab="Flow (m3/s)", ylab=Hab)
    abline(v=Qref, col = "green")
    abline(h=hab.Qref, col = "green", lty=2)
    abline(h=hab.Qnew, col = "orange", lty=2)
    #abline(v=Data[ThisRow, "MinQ"]*Qref, col = "orange")
    abline(v=QrefNew, col = "orange")
    legend("topright",xpd=TRUE,
           legend=c(paste("Habitat", Hab), paste("Qref = ",round(Qref, 3)), paste("Habitat at Qref =", round(hab.Qref, 2)),
                    paste("Qnew = ",round(QrefNew, 3)), paste("Habitat at Qnew = ", round(hab.Qnew, 2))),
           col = c("black", "green", "green","orange","orange"),   text.col = "black", lty = c(1, 1, 2,1,2), pch = c(-1, -1),   lwd = c(2,1,1,1,1,1),
           merge = TRUE, bg = 'white', inset = 0.05)
    #browser()
  }
  #browser()
  return(DeltaHab)
}

##############################################################################
