REC_McDowell<-function(MyREC=MyREC,myCI=0){
NO3N<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "NO3N", myCI = 0,Data=MyREC)) * 1000
DRP<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "DRP", myCI = 0,Data=MyREC)) * 1000
NH4N<-unlist(lapply(pick, predMcDowellWQ, myWQvar = "NH4N", myCI = 0,Data=MyREC)) * 1000

MyREC$log10NPratio <- log10((NO3N + NH4N)/DRP)
MyREC$TN <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "TN", myCI = myCI,Data=MyREC)) * 1000
MyREC$log10TN  <-  log10(MyREC$TN)
MyREC$TP  <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "TP", myCI = myCI,Data=MyREC)) * 1000
MyREC$SIN  <- 1000 * (unlist(lapply(pick, predMcDowellWQ, myWQvar = "NO3N", myCI = myCI,Data=MyREC)) + unlist(lapply(pick, predMcDowellWQ, myWQvar = "NH4N", myCI = myCI)))  
MyREC$CLAR  <- unlist(lapply(pick, predMcDowellWQ, myWQvar = "CLAR", myCI = -myCI,Data=MyREC))
return(MyREC)
}

REC_McDowell_scaled<-function(MyREC=MyREC,myCI=0){
  NO3N<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NO3N", myCI = 0,Data=MyREC)) * 1000
  DRP<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "DRP", myCI = 0,Data=MyREC)) * 1000
  NH4N<-unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NH4N", myCI = 0,Data=MyREC)) * 1000
  
  MyREC$log10NPratio <- log10((NO3N + NH4N)/DRP)
  MyREC$TN <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "TN", myCI = myCI,Data=MyREC)) * 1000
  MyREC$log10TN  <-  log10(MyREC$TN)
  MyREC$TP  <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "TP", myCI = myCI,Data=MyREC)) * 1000
  MyREC$SIN  <- 1000 * (unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NO3N", myCI = myCI,Data=MyREC)) + unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "NH4N", myCI = myCI)))  
  MyREC$CLAR  <- unlist(lapply(pick, predMcDowellWQ_scaled, myWQvar = "CLAR", myCI = -myCI,Data=MyREC))
  return(MyREC)
}