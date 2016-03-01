####################################################################
#                                                                  #
# CODE extend/ make predictions of WT levels for bores with short data  #
#                                                                  #
####################################################################
# Prepared by C. Fraser (2015)
#
####################################################################
#This code:
#     - reads in all ECAN wells data.  
#     - Removes silly data points (999 etc)
#     - Generates monthly average WT levels for all bores
#     - Allows user to:
#         - read in observed WT levles from a new bore
#         - predict WT levels in new bore for user specified dates
#           by correlating with bores from teh ECAN database and then 
#           developing a linear model to make predictions of WT level
#           (in some cases, this will involve two linear models - one between
#           user obs bore and most correlated site (CORR1 site), and then one between
#           most correlated site and another site (CORR2 site) that has observations
#           on the days of interest)
#     - Output are predictions on the days of interest

# To run: select all (ctrl+A) then press "Run"
rm(list=ls()) # clear memory

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   VARIABLES TO CHANGE:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Select months where predictions must be made (format: "MMM YYYY")
predM<-c("Sep 2004","Sep 2005","Apr 2006")
#Define minimum nubmer of concurrent observations required to evalaute correlation between bores (recommend >5)
nobs<-6
#Define WT buffer range for correlations (i.e. range of oBS bore WT levels +/- WTr in m)
WTr<-20
#Define name of csv file with the observation bore data
WTfile<-"ExampleWT.csv"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#################################################################
# 1. Prepare workspace (clear workspace, set directory, load packages, load functions and data)

rootdir<-"G:/Common Data/R_Scripts"
setwd(rootdir); source("SetDirectories.r"); dir<-SetDirectories(rootdir)

require(chron)                                                # times and dates
require(zoo)
#require(lattice)
#require(latticeExtra)

setwd(dir$fun)
source("Functions_For_GW_analysis.r")


#################################################################
# 2A. Tidy up data and make monthly average
#   - uncomment this step if you need to start from the beginning again, otherwise
#     use the preprocessed monthly average matrix

        #setwd(dir$data)
        #load("Canty_Wells2.RData")                  #This is a big matrix with all the well level observations
        
        #Data_n<-Data[which(Data$GWL<900&Data$GWL>-900),] #Remove samples that are crazy!
        #BoresALL<-unique(Data_n$BoreName)             #get a list of all the bore names
        #lBores2<-sapply(BoresALL,function(x) length(which(Data_n$BoreName==x)))  #Count how many observations per bore
        ##Don't keep bores with less than "thold" (here 5 is min number) number of observations
        #BoresIN<-as.character(BoresALL[which(lBores2>4)])
        
        ##Put the data for each bore as a separate item in a list (format easier for later analysis)
        #Data_m<-list()
        #for (i in 1:length(BoresIN)){
        #  Data_m[[i]]<-Data_n[(which(Data_n$BoreName==BoresIN[i])),]
        #}
        
        #Data_m2<-lapply(Data_m,function(x) zoo(x$GWL,x$chrono))  #Change into format for timeseries analysis
        ##Calculate monthly averages
        #MonthlyAvg <- Doug.cbind.list(lapply(Data_m2, function(x) aggregate(x, as.yearmon, mean)))  
  #names(MonthlyAvg)<-BoresIN  
  #save(MonthlyAvg,file="BoresMthAvg.Rdata")

#################################################################
# 2B. Load monthly average WT observation data and select subset for analysis
setwd(dir$data)
load("BoresMthAvg.Rdata")

#Identify the columns where there are data in the prediciton months
IN<-which(is.na(sapply(MonthlyAvg[match(predM,row.names(MonthlyAvg)),],sum))==F)

MonthlyAvg2<-MonthlyAvg[,IN]

#################################################################
#3. Load in data from new bore
Obs<-read.csv(WTfile,colClasses=c("character","numeric"))
#Add a column with a Month Year format
Obs$Mth<-as.character(as.yearmon(as.character(sapply(Obs$Date,function(x) strsplit(x,split=" ")[[1]][1])),format="%d/%m/%Y"))
Obs$D<-chron(as.character(sapply(Obs$Date,function(x) strsplit(x,split=" ")[[1]][1])),format=c(dates="d/m/y"))

#Take averages if multiple observations in a month
OBS2<-as.data.frame(aggregate(zoo(Obs$WT,Obs$D), as.yearmon, mean))
names(OBS2)<-"WT"

#################################################################
# Identify wells with a minimum of "nobs" concurrent observations and within same Wt ranges

iMth<-match(row.names(OBS2),row.names(MonthlyAvg))
# in case there are no matches:
iobs<-which(is.na(iMth)==F)
iMth<-iMth[which(is.na(iMth)==F)]

#Generate subset of obserations
  #Identify the number of concurrent observations
in2<-sapply(MonthlyAvg[iMth,],function(x) length(which(is.na(x)==F)))


#Define an acceptable WT level range
WTmin<-min(OBS2$WT)-WTr
WTmax<-max(OBS2$WT)+WTr
WTin<-which(apply(MonthlyAvg,2,min,na.rm=T)>WTmin & apply(MonthlyAvg,2,max,na.rm=T)<WTmax)

#Take subset of data
MonthlyAvg3<-MonthlyAvg[iMth,intersect(WTin,which(in2>=nobs))]

#################################################################
# calculate correlation between the well of interest and the remaining wells

CorWT<-sort(cor(cbind(OBS2$WT[iobs],MonthlyAvg3),use="pairwise.complete.obs")[-1,1],decreasing=TRUE)
best5<-names(CorWT[1:5])


#Define a function to develop linear models
Mydo.lm<-function(x1,x2,pred=NULL){
  lmdata<-as.data.frame(cbind(x1[which(is.na(x2)==F)],x2[which(is.na(x2)==F)]))
  names(lmdata)<-c("MyBore","corBore")
  WTmod<-lm(MyBore~corBore,data=lmdata)
  
  y<-lmdata$MyBore
  yhat<-predict(WTmod,data=lmdata)
  
  #Calculate performance statistics of the linear model
  NSE <- (1 -  sum((y-yhat)^2)/sum((y-mean(y))^2))    #Nash sutcliffe efficiency
  RMSD <- sqrt(1/(length(y)-1) * (sum((yhat-y)^2)))   #Root mean square deviation
  OUTPUT<-list(mod=WTmod,NSE=NSE,RMSD=RMSD)
}

#Generate a model between the observations and the most correlated site
t1<-Mydo.lm(x1=OBS2$WT,x2=MonthlyAvg3[,best5[1]])

#Check - are there osbervations avaiable with the most correlated site on the prediction dates?
if(length(which(is.na(MonthlyAvg[predM,best5[1]])==F))<length(predM)){
  #Case 1: The correlated bore does not cover all the desired dates
  CorWT2<-sort(cor(cbind(MonthlyAvg[,best5[1]],MonthlyAvg2),use="pairwise.complete.obs")[-1,1],decreasing=TRUE)
  Thisone<-names(CorWT2[1])  #the bore that the corretlated site is most correlated to
  
  ix<-which(is.na(MonthlyAvg[,best5[1]])==F)
  t2<-Mydo.lm(x1=MonthlyAvg[ix,best5[1]],x2=MonthlyAvg2[ix,Thisone])
  
  
  lmdata1<-as.data.frame(MonthlyAvg2[predM,Thisone]); names(lmdata1)<-"corBore"
  #Make predicted WT levels for the dates of interest
  lmdata2<-as.data.frame(predict(t2[["mod"]],newdata=lmdata1));names(lmdata2)<-"corBore"
  #Replace predictions with observations if possible:
  lmdata2[which(is.na(MonthlyAvg[predM,best5[1]])==F),1]<-MonthlyAvg[predM[which(is.na(MonthlyAvg[predM,best5[1]])==F)],best5[1]]
  
  
  #Now Make predictions for teh origial observation bore
  PREDWT<-predict(t1[[1]],newdata=lmdata2)
  names(PREDWT)<-predM
  print(paste("R2 for OBS to Corr1 model", t1[["NSE"]]))
  print(paste("R2 for Corr1 to Corr2 model", t2[["NSE"]]))
  print("The top 5 correlated sites are:")
  print(best5)
  print(paste("Corr1 site is",best5[[1]]))
  print(paste("Corr2 site is",Thisone))
  print(PREDWT)
  
}else{
  lmdata1<-as.data.frame(MonthlyAvg2[predM,best5[1]]); names(lmdata1)<-"corBore"
  PREDWT<-predict(t1[[1]],newdata=lmdata2)
  names(PREDWT)<-predM
  
  print(paste("R2 for OBS to Corr1 model", t1[["NSE"]]))
  print("The top 5 correlated sites are:")
  print(best5)
  print(paste("Corr1 site is",best5[[1]]))
  print(PREDWT)
  
  
}

