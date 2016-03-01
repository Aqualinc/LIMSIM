#________________________________________________________________________________________________________________
# CODE TO IDENITFY REACHES UPSTREAM OF A DEFINED REACH
#________________________________________________________________________________________________________________
FindCatchment<-function(EndR,NewREC,MinOrder=0) {  #Where EndR is the NZreach number at the (sub)catchment outlet
#Start by limiting the values to search to those within the same region
pick <- NewREC$NZReach[NewREC$NZReach > floor(EndR/100000)*100000 & NewREC$NZReach < ceiling(EndR/100000)*100000]  # make sure no NAs in NZReach
pick <- pick[!is.na(pick)]  # the NZReaches in the domain
pick<-NewREC$NZReach[NewREC$ORDER > MinOrder]
MyREC <- NewREC[match(pick, NewREC$NZReach),]

#Identify the row of the end reach
EndR_ind<-which(MyREC$NZReach==EndR)

#Identify the segments upstream of the end reach
UPST<-which(MyREC$Nztnode==MyREC$Nzfnode[EndR_ind])
CatchSel<-MyREC$NZReach[UPST];

#Sequentially work upstream
while (length(UPST)>0) {
  UPSTn<-c()
  for (i in 1:length(UPST)){
      UPSTn<-c(UPSTn,as.vector(which(MyREC$Nztnode==MyREC$Nzfnode[UPST[i]])))
  }
  CatchSel<-c(CatchSel,MyREC$NZReach[UPSTn])
  UPST<-UPSTn
}

return(CatchSel)
}

