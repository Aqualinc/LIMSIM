###################################################################################
###################################################################################

#                         HABITAT FUNCTIONS                                       #

###################################################################################
###################################################################################
# Functions indcluded in this script:
#       - IntegratedWidth Calculates the reduction in width over the whole FDC
#       - habitat         Calculates the available habitat for a prescribed species
#       - delta.hab       Calcualtes the change in habitat


###################################################################################
###################################################################################

Run_HabitatFunctions<-function(MyREC=NULL,pick=pick,MinQ=0,method="Jowett"){

  if (method=="Booker"){
    print("Calculating BOOKER QW relationships...")
    load(paste(dir$fun,"BookerModel.RData",sep="/"))  # this loads Doug's fitted models to estimate width using catchment area and climate category 
    QW <<- lapply(pick, GetWidth, method="Booker",Data=MyREC,IntModelCl=IntModelCl,SlopeModelCl=SlopeModelCl,QuadModelCl=SlopeModelCl) # this gets the width flow relationships for all reaches
  }else{
    print("Calculating JOWETT QW relationships...")
    QW <<- lapply(pick, GetWidth, method="Jowett",Data=MyREC) 
  }
#QW2 <<- lapply(pick, GetWidth, method="booker") # this gets the width flow relationships for all reaches.
names(QW) <- pick   # names the list elements
# the mean reduction in width for flows up to the mean (the largest value in QW
print("Calculating reduction in width...")
#browser()
MyREC$RedWidth<- sapply(pick, IntegratedWidth, FDC = MyFDC, FlowWidth = QW, minFlow=NULL, prop=NULL,
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


###############################################################################
#   Compute QW curve for an NZReach  Using Booker 2009                        #
#   Compute QW curve for an NZReach  Using Jowett 1998                        #
###############################################################################
GetWidth <- function(ThisNZReach = 13524724, method="Booker", n.pairs = 10, Data=MyREC,IntModelCl=NA,SlopeModelCl=NA,QuadModelCl=NA) {   # Pass data in one line (=NZReach)
  # this function calculates Q and W in n.pairs from 0 to Qbar(mean flow)
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach) # get the parameters from the FDC dataset is this faster???
  this.seg <- Data[ThisRow, ]
  #browser()
  #Qbar <- this.seg$Flow_L_s/1000
  Qbar <- Data[ThisRow, "MeanFlow"]  # OR take mean flow from Doug's estimates
  if(is.na(Qbar)==F){
    Q <- seq(Qbar*0.0001, Qbar, length.out=n.pairs)  # better not to have flow of zer0
    if (method=="Booker") {
      this.seg$LogCatchmentArea <- log10(this.seg$usArea/1000000)
      this.seg$LogFlow <- log10(Qbar)    
      d0<- predict(IntModelCl, newdata = this.seg)
      d1 <- predict(SlopeModelCl, newdata = this.seg)
      d2 <- predict(QuadModelCl, newdata = this.seg)
      LogQ <- log10(Q)
      W  <-  10^(d0 +  LogQ*d1  +  d2*(LogQ^2))
      QW <- data.frame(cbind(Q, W))
    } else {
      Wbar <- 7.76*(Qbar^0.488)           #  mean width and meam flow from Downstream Hydraulic Geometry  Jowett 1998
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
#' @param FDC a dataframe of flow rates for different percentiles columns) for different reaches (rows)
#' @param FlowWidth dataframe of flow vs width
#' @param minFlow the minimum flow in cumecs. If this is NULL then it is calculated from the MinQ attribute in the REC table. Default is NULL
#' @param prop the minimum flow in fraction of the "Qref". prop stands for proportion. Defaults to 0.8
#' @param Qref the reference flow in cumecs which the "prop" and "allocate" are fractions of. If it is NULL, then it is set to the "MALF" attribute of the reach. Defaults to NULL
#' @param allocate the allocation as a fraction of "Qref". If NULL then it is set to the "AllocQ" attribute of the reach. Defaults to 0.5
#' @param Plot whether to plot the flow duration curve with lines showing the minimum and managed flows in normal and log space. Defaults to TRUE
#' @param Data the REC attribute table
#' @return The average width loss resulting from the allocations
#' @keywords REC
#' @export
#' @examples
#' IntegratedWidth()

IntegratedWidth <- function(ThisNZReach = 13524724, FDC = MyFDC, FlowWidth = QW, minFlow=NULL, prop=0.8,
                            Qref=NULL, allocate=0.5, Plot = TRUE, Data = MyREC) {
  
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach)                        # Get the row number of the REC attribute table for the reach of interest 
  if(Data$TheTake[ThisRow] == 0 | any(is.na(FlowWidth[[ThisRow]]))) {  # if there is no take or there is an NA in the FlowWidth table, skip entirely width loss is ZERO  -> NA
    AveWidthLoss <- NA
  } else {
    
    #Sort out the values for the minimum flow and allocation flow in cumeces
    if (is.null(Qref)) Qref <- Data$MALF[ThisRow]                      # If the "Qref" parameter is NULL, set it to the MALf from the REC attributes table
    if (is.null(minFlow)) minFlow <- Qref * Data$MinQ[ThisRow]         # If the "minflow" parameter is NULL, set it based on the REC attribute tables "MinQ" and "MALF"
    
    if(!is.null(Data$TheTake))   {                                     # Check to see if take has been established
      allocation <-  Data$TheTake[ThisRow]                             # If it has, tehn use it to set the allocation
    } else {
      allocation <- Data$AllocQ[ThisRow] * Qref                        # if not then the allocation will be Qref * allocate. The DEFAULT for Qref= MALF, this allows for rules of thumb to be used.
    }
    
    if (is.na(FDC[ThisRow, ])){
      AveWidthLoss <-NA
    }else{
      # obtain FDC data for this REACH from Gauges
      FDC.data <- FDC[ThisRow, ]                                      # THIS is the FDC
      freqs    <- 100*Perc                                            # the percentiles the FDC is estimated for from the global value of Perc
      
      #browser()
      if (minFlow < FDC.data[1]){minFlow <- FDC.data[1]}              # Make sure the minimum flow is at least as large as the smallest value in the flow duration curve
      #y=freqs; x=FDC.data
      
      freq.diff <-  approx(y=freqs, x=FDC.data, xout=c(minFlow+allocation, minFlow))$y # estimate the percentiles for which the "restricted takes" and "stopped takes" occur
      #if (length(freq.diff)!=2) {browser()}
      GWTakeFDCShift <- Data$GWAlloc[ThisRow] * Data$TrueBFI[ThisRow] * Data$AqBaseFlow[ThisRow] * Data$BaseFlow[ThisRow]

      AlteredFDC<-FDC.data - approx(x=c(100,freq.diff,0),y=c(allocation+GWTakeFDCShift,allocation+GWTakeFDCShift,GWTakeFDCShift,GWTakeFDCShift),xout=freqs)$y
      #browser()
      Qa <- FlowWidth[[ThisRow]]$Q                                    # the QW calculations corresponding to the reach of interest
      W  <- FlowWidth[[ThisRow]]$W
      
      NaturalWidths <- approx(x=c(0,Qa), y=c(0,W), xout=FDC.data)$y   # interploate W vs Q data
      AlteredWidths <- approx(x=c(0,Qa), y=c(0,W), xout=AlteredFDC)$y # interploate W vs Q data
      
      AveWidthLoss <-   mean((AlteredWidths - NaturalWidths)/NaturalWidths, na.rm=T) *100 # mean reduction in width (percentage of natural flow) 
    }
    
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
  }  
  
  return(AveWidthLoss)   #returns the loss of width
}# end

NaturalWidth <- function(ThisNZReach = 13524724, FDC = MyFDC, FlowWidth = QW, Data = MyREC) {
  
  ThisRow <- which(Data$NZReach == ThisNZReach) # 
  if (is.na(FDC[ThisRow, ])==T){
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
  if(any(sp%in%"Brown trout adult")) hab.par <- c(1.17, 4.35)   #hab.par = parameters for generalised habitat models C and K
  if(any(sp%in%"Redfin"))            hab.par <- c(0.26, 7.39)   # redfin
  if(any(sp%in%"ComBully"))           hab.par <- c(0.39, 6.51)   # common bully
  if(any(sp%in%"Fry"))           hab.par <- c(0.86, 10.21)   # brown trout fry
  if(any(sp%in%"Spawn"))           hab.par <- c(1.24, 9.89)   # brown trout spawning
  if(any(sp%in%"Deli"))           hab.par <- c(0.33, 1.92)   # deleatidium
  if(any(sp%in%"LongEel"))           hab.par <- c(0.07, 2.07)   # longfin eel
  if(any(sp%in%"Torrent"))           hab.par <- c(0.88, 4.05)   # torrent fish
  if(any(sp%in%"ShortEel"))           hab.par <- c(0.13, 2.32)   # shortfin eel
  if(any(sp%in%"Bluegill Bully"))           hab.par <- c(1.01, 6.13)   # Bluegill bully
  if(any(sp%in%"Inanga"))           hab.par <- c(0.19, 19.74)   # INanaga
  if(any(sp%in%"Kokopu"))           hab.par <- c(0.03, 2.29)   # Kokopu
  if(any(sp%in%"Upland Bully"))    hab.par <- c(0.19, 13.13)   # upland bully
  
  ThisRow <- which(Data$NZReach == ThisNZReach) #
  if (is.na(FDC[ThisRow, ])==T){
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

###############################################################################
#   calculate delta habitat with delta Q for an NZreach                       #
###############################################################################
# WH can be HV or WUA
delta.hab <- function(ThisNZReach = 13524724, Data = MyREC, FlowHab = QWH, Hab = "WUA", Qref=NULL, prop=0.8, Plot = TRUE) {  # returns the widths at Qref and prop times Qref and plots this QW curve if plot=TRUE
  #browser()
  ThisRow <- which(Data$NZReach == ThisNZReach)
  if (is.null(Qref)) { Qref <- Data[ThisRow, "MALF"] } # the reference flow is the MALf from the hyd predictions data
  #browser()
  Q <- FlowHab[[ThisRow]]$Q   #  get flow vs habitat
  H <- unlist(FlowHab[[ThisRow]][Hab])

  
  
  if(Qref==0|is.na(Qref)==T| Data$TheTake[ThisRow]==0) {  # if Qref is zero the habitat is ZERO and so is delta-habitat
    hab.Qref <- NA
    hab.Qnew <- NA
    DeltaHab <- NA
  } else {
    
    habs <-  approx(x=Q, y=H, xout=c(Qref, Data[ThisRow, "MinQ"]*Qref))$y # interploate W vs Q data
    hab.Qref <- habs[1]
    hab.Qnew <- habs[2]
    DeltaHab <- hab.Qnew/hab.Qref*100
  } # end if
  
  if (Plot==TRUE) {
    par(bg="grey90")
    plot(Q, H, type="l", lwd = 2,
         main=paste("Difference in", Hab, "for", Data[ThisRow, "MinQ"], "times Qref"),
         xlab="Flow (m3/s)", ylab=Hab)
    abline(v=Qref, col = "green")
    abline(h=hab.Qref, col = "green", lty=2)
    abline(h=hab.Qnew, col = "orange", lty=2)
    abline(v=Data[ThisRow, "MinQ"]*Qref, col = "orange")
    
    legend("bottomright",
           legend=c(paste("Habitat", Hab), paste("Qref = ",round(Qref, 3)), paste("Habitat at Qref =", round(hab.Qref, 2)),
                    paste("Qnew = ",round(Data[ThisRow, "MinQ"]*Qref, 3)), paste("Habitat at Qnew = ", round(hab.Qnew, 2))),
           col = c("black", "green", "green","orange","orange"),   text.col = "black", lty = c(1, 1, 2,1,2), pch = c(-1, -1),   lwd = c(2,1,1,1,1,1),
           merge = TRUE, bg = 'white', inset = .05)
  }
  return(DeltaHab)
}

##############################################################################
