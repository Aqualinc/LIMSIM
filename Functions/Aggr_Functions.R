#######################################################################
#             AGGREGATION FUNCTIONS for REC                           #
#######################################################################

setwd(dir$data)
load("Aggr_Data.RData")   #This dataset includes a list called "Reach" that lists the upst reaches from all segments
load("REC_ids.Rdata")
 #~~~~~~~~~~~ see end of file for instructions about how to create "Reach" and "IDSort" ~~~~~~~~~~~~~~~~~~~~~~~~~~~


######################################################################################
#  Downstream Aggregation of a single variable                                       #
######################################################################################
#Example call: MyREC$PloadSUM<-DstAggr_OneVar(IDSort=MyREC$IDSort,Aggrvar=MyREC$Pload,Reach=Reach,IDnew=match(MyREC$NZReach,IDs))
DstAggr_OneVar=function(IDSort=NULL,Aggrvar=NULL,Reach=Reach,IDnew=NULL){
  #Note: the use of IDnew is because my verisons of IDsort and Reaches are for the whole country.  IDnew allows subsets of the country to be processed
  val<-sort(unique(IDSort))                                 #Determines the order to do the calculations
  OUTPUT<-vector(mode="numeric",length=length(Reach)) *NA                                      #Initiate output matrix
  
  OUTPUT[IDnew[which(IDSort==1)]]<-Aggrvar[which(IDSort==1)]  #Fill in the 1st order streams (which are just the local values)
  print("Start accumulation....")
  for (i in 2:length(val)){
    index2<-which(IDSort==val[i])
    index<-IDnew[index2]
    for (j in 1:length(index)){
      OUTPUT[index[j]]<-sum(OUTPUT[Reach[[index[j]]]])+Aggrvar[index2[j]] 
    } 
  }
  OUTPUT<-OUTPUT[IDnew];OUTPUT[OUTPUT==0]<-NA
  return(OUTPUT)
}

######################################################################################
#  Downstrean Aggregation with area weighting (gives UPST average)                   #
######################################################################################
#Example call: 
DstAggr_AreaAvg=function(IDSort=NULL,Aggrvar=NULL,segArea=NULL,upstArea=NULL,Reach=Reach,IDnew=NULL){
  
  val<-sort(unique(IDSort))                                 #Determines the order to do the calculations
  OUTPUT<-vector(mode="numeric",length=length(Reach)) *NA   #Initiate output matrix
  OUTPUT[IDnew[which(IDSort==1)]]<-Aggrvar[which(IDSort==1)] #Fill in the 1st order streams (which are just the local values)
  print("Start accumulation....")
  for (i in 2:length(val)){
    index2<-which(IDSort==val[i])
    index<-IDnew[index2]
    for (j in 1:length(index)){
      OUTPUT[index[j]]<-(sum(OUTPUT[Reach[[index[j]]]]*upstArea[[IDnew[index[j]]]])+Aggrvar[index2[j]]*segArea[index2[j]])/upstArea[index2[j]] 
    } 
  }
  OUTPUT<-OUTPUT[IDnew];OUTPUT[OUTPUT==0]<-NA
  return(OUTPUT)
}

######################################################################################
#  Identify all reaches upstream of a specified location (delineate catchment)       #
######################################################################################
FindCatchment=function(EndR,NewREC) {  #Where EndR is the NZreach number at the (sub)catchment outlet
  #Start by limiting the values to search to those within the same region
  pick <- NewREC$NZReach[NewREC$NZReach > floor(EndR/100000)*100000 & NewREC$NZReach < ceiling(EndR/100000)*100000]  # make sure no NAs in NZReach
  pick <- pick[!is.na(pick)]  # the NZReaches in the domain
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

######################################################################################
#  Identify all reaches downsteram for each reach of REC                            #
######################################################################################
#This must be done before information can be routed upstream
#The code is slow, so recommened running once, then saving for use in the UPSTaggr function

FindDstReaches=function(REC=REC){
StrPath<-matrix(NA,nrow=length(REC[,1]),ncol=length(unique(REC$IDSort)))

IDs<-sort(unique(REC$IDSort),decreasing=TRUE)
for (i in 1:length(IDs)){
  print(IDs[i])
  ids<-which(REC$IDSort==IDs[i])
  for (j in 1:length(ids)){
    UPST<-FindCatchment(EndR=REC$NZReach[ids[j]],NewREC=REC) 
    StrPath[match(UPST,REC$NZReach),length(REC$IDsort)+1-i]<-REC$NZReach[ids[j]]
  }
} 
StrPath2<-cbind(REC$NZReach,StrPath)  #Add in the local reach to the list as well
return(StrPath2)
}

######################################################################################
#  Assign the highest or lowest downstream value to each reach                       #
######################################################################################
#call expample: tracked_CLAR<-apply(StrPath2,1,trackback_maxmat,REC=REC,var="CLAR"))
trackback_maxmat=function(ThisRow=NULL,REC=NULL,var="TN"){
  val<-max(REC[match(ThisRow,REC$NZReach),var],na.rm=T)    
  return(val)
}
trackback_minmat=function(ThisRow=NULL,REC=NULL,var="TN"){
  val<-min(REC[match(ThisRow,REC$NZReach),var],na.rm=T)    
  return(val)
}


######################################################################################
#  Creating "Reaches" and "IDSort"                                                   #
######################################################################################

#Reach<-list()
#LReach<-seq(1,length(pick),1)
#for (i in 1:length(pick)){
#  Reach[[i]]<-which(MyREC$Nztnode==MyREC$Nzfnode[i])    #Note - to avoid some issues that i encoutered, I would actually make Reach as the ID names rather than just the indexing if i did it again!
#  LReach[i]<-length(Reach[[i]])
#}


#~~~~~~~~~~~~~~~~~Making IDSort~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#IDSort describes the order that reaches must be evaluated in, so that downstream accumulation works properly!

##Initialise the output matrix
#IDRiver<-seq(1,length(pick),1)            
#IDSort<-IDRiver*0

##First assign a "1" to the most upstream reaches
#num<-1
#IDSort[LReach==0]<-num
#in1<-sum(IDSort)

##Assign the order (IDSort) that the reaches must be calculated in
#num<-num+1
#index<-which(IDSort==0&MyREC$ORDER<=num)
#index2<-which(IDSort!=0) #These are the most upstream catchments

#while (length(pick)-in1>0) {
#  ind<-seq(1,length(index),1)
#  Lmatch<-unlist(lapply(ind,function(x) length(which(!is.na(match(Reach[[index[x]]],index2))==T))))
#  i2<-which(LReach[index]==Lmatch)
#  IDSort[index[i2]]<-num
#  in1<-in1+length(i2)
  
#  num<-num+1
#  index<-which(IDSort==0&MyREC$ORDER<=num)
#  index2<-which(IDSort!=0) #These are the most upstream catchments
  
#  print(in1/length(pick)*100)
#}