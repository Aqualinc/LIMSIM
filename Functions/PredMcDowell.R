#setwd(dir$data)
RECPastWQ <-read.csv(paste(dir$data,"BaseLineWQ.csv",sep="/"))
# store as an array
Vars <- as.character(unique(RECPastWQ$Var)) # the water quality variables that we can model
RECclass <- names(RECPastWQ)[c(4:21)]       # the REC classes that

meanWQ <- array(0, dim=c(21,18,length(Vars)), dimnames = list(seq(0,100,5),RECclass,Vars))
seWQ <- array(0, dim=c(21,18,length(Vars)), dimnames = list(seq(0,100,5),RECclass,Vars))

RECPastWQ$Type <- paste(RECPastWQ$Var, RECPastWQ$est.or.se, sep=".")

# load the data into the array!
for(i in 1:length(Vars)) meanWQ[,,Vars[i]] <- as.matrix(RECPastWQ[RECPastWQ$Var == Vars[i] &  RECPastWQ$est.or.se == "est", c(4:21)])
for(i in 1:length(Vars)) seWQ[,,Vars[i]] <- as.matrix(RECPastWQ[RECPastWQ$Var == Vars[i] &  RECPastWQ$est.or.se == "se", c(4:21)])


# function to compute any WQ concentration + or - error for any NZReach based on McDowell et al.
predMcDowellWQ <- function(ThisNZReach = 13524724,  myWQvar = "DRP", myCI = 0, RECclasses = unlist(dimnames(meanWQ)[2]), Data = MyREC) {
  mySoF = as.character(Data[Data$NZReach == ThisNZReach, "SoF"])
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

# function to compute any WQ concentration + or - error for any NZReach based on McDowell et al.
predMcDowellWQ_scaled <- function(ThisNZReach = 13524724,  myWQvar = "DRP", myCI = 0, RECclasses = unlist(dimnames(meanWQ)[2]), Data = MyREC) {
  mySoF = as.character(Data[Data$NZReach == ThisNZReach, "SoF"])
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