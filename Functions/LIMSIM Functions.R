#DEFINE THE VARIABLE NAMES, LABELS and LEGEND NAMES
varn<<-c("TheTake","freq.min","freq.man","freq.reliability","IrriAreaIrrigated","AllocQUse","CLAR","MCI","RedWidth","MeanFilsDeltaQ","MaxFilsDeltaQ","MeanMatsDeltaQ","MaxMatsDeltaQ","deltaEEL","deltaShortEEL","deltaTrout","deltaBully","deltaInanga","deltaTorrent","deltaKokopu","deltaComBully","deltaTroutFry","deltaTroutSpawn","deltaUplandBully")
labeln<<-c("Irrigation Take (m3/s)","Minimum flow restrictions (% of time)", "Management flow restrictions (% of time)","Irrigation Bulk Reliability","% of Irrigable area Irrigated","% of Allocation used",
         "Clarity","MCI","% Reduction in river width","Mean Periphyton Filaments","Max Periphyton Filaments","Mean Periphyton Mats","Max Periphyton Mats","% Longfin Eel habitat","% Shortfin Eel habitat","% Brown Trout habitat","% Bluegill Bully habitat","% Inanga habitat","% Torrent habitat","% Kokopu habitat","% Common Bully habitat","% Brown Trout Fry Habitat","% Brown Trout Spawning Habitat","% Upland Bully Habitat")
legname<<-c("take m3/s","% time","% time","% available","% area","% used","Clarity","MCI","dW (%)","Mean Fils","Max Fils","Mean Mats","Max Mats","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline","% of baseline")
#coldir<<-c(F,F,F,F,F,F,F,F,F,T,T,T,T,F,F,F,F)
coldir<<-c(F,T,T,F,F,T,F,F,F,T,T,T,T,F,F,F,F,F,F,F,F,F,F,F)

###################################################################################################################
####################################################################################################################
PopulateScenarios<-function(MyREC=MyREC,MinQ=MinQ,AllocQ=AllocQ,ShadingON=0,Shading=0,LandManagement=LandManagement,IrrigableAreaTarget=IrrigableAreaTarget,PropIrri=0){

MyREC1<-MyREC;MyREC2<-MyREC;MyREC3<-MyREC
#FIRST Flows and abstractions
MyREC1<-Run_FlowandAbstraction(MyREC=MyREC1,AllocQ=AllocQ[1],MinQ=MinQ[1],pick=pick)
MyREC2<-Run_FlowandAbstraction(MyREC=MyREC2,AllocQ=AllocQ[2],MinQ=MinQ[2],pick=pick)
MyREC3<-Run_FlowandAbstraction(MyREC=MyREC3,AllocQ=AllocQ[3],MinQ=MinQ[3],pick=pick)

if (ShadingON==1){
  #Can optionally recalculate shading
  if(!is.na(Shading[1])==T) MyREC1$segShade<-Shading[1];MyREC1$PAR<-53.2 * MyREC1$L * 0.01 * (100 -  MyREC1$segShade*100)
  if(!is.na(Shading[2])==T) MyREC2$segShade<-Shading[2];MyREC2$PAR<-53.2 * MyREC2$L * 0.01 * (100 -  MyREC2$segShade*100)
  if(!is.na(Shading[3])==T) MyREC3$segShade<-Shading[3];MyREC3$PAR<-53.2 * MyREC3$L * 0.01 * (100 -  MyREC3$segShade*100)
}
#Second Water quality, periphyton and MCI
MyREC1<-Run_WQandPeri(MyREC=MyREC1,LandManagement=LandManagement[1],IrrigableAreaTarget=IrrigableAreaTarget[1],AllocQ=AllocQ[1],PropIrri=PropIrri)
MyREC2<-Run_WQandPeri(MyREC=MyREC2,LandManagement=LandManagement[2],IrrigableAreaTarget=IrrigableAreaTarget[2],AllocQ=AllocQ[2],PropIrri=PropIrri)
MyREC3<-Run_WQandPeri(MyREC=MyREC3,LandManagement=LandManagement[3],IrrigableAreaTarget=IrrigableAreaTarget[3],AllocQ=AllocQ[3],PropIrri=PropIrri)

#Third Habitat Functions
MyREC1<-Run_HabitatFunctions(MyREC=MyREC1,MinQ=MinQ[1])
MyREC2<-Run_HabitatFunctions(MyREC=MyREC2,MinQ=MinQ[2])
MyREC3<-Run_HabitatFunctions(MyREC=MyREC3,MinQ=MinQ[3])

AllMyREC<-list(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3)
return(AllMyREC)
}
###################################################################################################################
# MakePick is a function that returns all the NZReach numbers in an REC table with the option to restrict by strahler order or region
####################################################################################################################
MakePick<-function(MyREC=NULL,MinOrder=5,Region=5){
  if (Region<16){
pick<-MyREC$NZReach[MyREC$ORDER > MinOrder & MyREC$NZReach >= (Region) * 1000000 & MyREC$NZReach < (Region+1) * 1000000]
  }else{ #Otherwise the whole country is selected
    pick<-MyREC$NZReach[MyREC$ORDER > MinOrder]
  }
return(pick)
}
###################################################################################################################
####################################################################################################################
MakePlots<-function(MyREC1=NULL,MyREC2=NULL,MyREC3=NULL,Indicator=1,PlotType=1,Scenarios=NULL){
  if (PlotType==1){ #Create a map 
   MultiMap(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn=varn[Indicator],varname=labeln[Indicator],legname=legname[Indicator],Scenarios=Scenarios,Ind=Indicator)
   
 }else{
  if (PlotType==2){ #Create a histogram
    ScenarioHistograms(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn=varn[Indicator],labeln=labeln[Indicator],Limit=NULL,Scenarios=Scenarios)
    }else{
    ScenarioDensityPlots_by(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn=varn[Indicator],labeln=labeln[Indicator],Scenarios=Scenarios)
        
    }
 }
  
}
###################################################################################################################
####################################################################################################################

ScenarioDensityPlots_by<-function(MyREC1=0,MyREC2=MyREC2,MyREC3=MyREC3,varn="MeanFils",by_var="SimpleClass",labeln="Mean Cover by Filaments",Scenarios){
BasicResultsStack<-GetBasicResults(MyREC1,MyREC2,MyREC3,varn,by_var,Scenarios) 
BasicResultsStack<-BasicResultsStack[which(is.na(BasicResultsStack$values)==F),]
x11(); densityplot(~values | SimpleClass, groups = ind, data = BasicResultsStack,  auto.key = T, scales = list(relation = "free"),
                   xlab = labeln, ylab = "Proportion of locations", main=paste("Predicted ",labeln,"for NZ rivers",sep=" "), plot.points = F)  #  !
#addLine(v=mean(MyREC2[,var], na.rm=T), col="red", lty=2) 
#addLine(v=mean(MyREC3[,var], na.rm=T), col="blue", lty=2) 
#if (MyREC1!=0){addLine(v=mean(MyREC1[,var], na.rm=T), col="green", lty=2)}

}

###################################################################################################################
####################################################################################################################

ScenarioHistograms<-function(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn="MeanFils",by_var="SimpleClass",labeln="Mean Cover by Filaments",Limit=NULL,Scenarios=NULL){
#This function creates a plto with 3 histograms representing each of the low, medium and high scenarios for a specified indicator  
BasicResultsStack<-GetBasicResults(MyREC1,MyREC2,MyREC3,varn,by_var,Scenarios)
nplot<-length(unique(BasicResultsStack$ind))
x11(); histogram(~values | ScenariosFac, data = BasicResultsStack,  auto.key = T, scales = list(relation = "free"),layout=c(1,nplot),xlab = labeln, ylab = "No. of locations", main=labeln, plot.points = F,panel=function(x,...){panel.histogram(x, ...) 
                                                                                                                                                                                                                                panel.abline(v=median(x,na.rm=T),col="red")}) #  !
#if (!is.null(Limit)){  #Add a vertical line to represent the limit if specified
#  addLine(v=Limit, col="red", lty=2, lwd=3) 
 # }
}
###################################################################################################################
####################################################################################################################

GetBasicResults<-function(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn="MeanFils",by_var=SimpleClass,Scenarios=NULL){
  if (length(Scenarios)==3){
    BasicResults <- cbind.data.frame(MyREC1[,varn], MyREC2[,varn], MyREC3[,varn])   
    names(BasicResults) <-c("Low", "Medium", "High")
  }else{
    if (length(Scenarios)==2){
      if (sum(Scenarios)==3){
        BasicResults <- cbind.data.frame(MyREC1[,varn], MyREC2[,varn])   
        names(BasicResults) <-c("Low", "Medium")
      }else{
        if (sum(Scenarios)==4){
          BasicResults <- cbind.data.frame(MyREC1[,varn], MyREC3[,varn])   
          names(BasicResults) <-c("Low", "High")
        }else{
          BasicResults <- cbind.data.frame(MyREC2[,varn], MyREC3[,varn])   
          names(BasicResults) <-c( "Medium", "High")
        }
      }
    }else{
      if (sum(Scenarios)==1){
        BasicResults <- MyREC1[,varn]
        names(BasicResults) <-c("Low")
      }else{
        if (sum(Scenarios)==2){
          BasicResults <-MyREC2[,varn]
          names(BasicResults) <-c("Low")
        }else{
          BasicResults <- MyREC3[,varn]
          names(BasicResults) <-c("High")
        }
      }
    }
  }
  
  BasicResultsStack <- stack(BasicResults)
  BasicResultsStack$ScenariosFac <- factor(BasicResultsStack$ind, levels =c("High", "Medium", "Low"))
  BasicResultsStack$SimpleClass <- MyREC1[,by_var]
  return(BasicResultsStack)
}

###################################################################################################################
####################################################################################################################

MultiMap<-function(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn="MeanFils",varname="Fil Periphyton",legname="Cover",Scenarios=c(1,2,3),Ind=1){
  #This function creates a single plot with 3 maps showing the results from the low medium and high scenarios
  h=20
  if (length(MyREC1[,1])>10000){MinOrder<-3}else{MinOrder<-0}
  if (length(Scenarios)==3){
  x11(height=h, width=60) ; par(mfrow=c(1,3))
  #browser()
  MyBreaks<-unique(as.vector(quantile(c(MyREC1[,varn],MyREC2[,varn],MyREC3[,varn]), probs = seq(0, 1, by = 1/15), na.rm = TRUE)))
  #MyBreaks<-seq(min(c(MyREC1[,varn],MyREC2[,varn],MyREC3[,varn]),na.rm=T),max(c(MyREC1[,varn],MyREC2[,varn],MyREC3[,varn]),na.rm=T),length.out=10)
  MapRivers(RECvar=MyREC1[,varn], REC=MyREC1,MinOrder=MinOrder, pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 1",sep=" "), MyBreaks=MyBreaks,col.ord=coldir[Ind]) 
  MapRivers(RECvar=MyREC2[,varn], REC=MyREC2, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 2",sep=" "),MyBreaks=MyBreaks,col.ord=coldir[Ind])
  MapRivers(RECvar=MyREC3[,varn], REC=MyREC3,MinOrder=MinOrder, pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 3",sep=" "),  MyBreaks=MyBreaks,col.ord=coldir[Ind]) 
  }else{
    if (length(Scenarios)==2){
      x11(height=h, width=40);par(mfrow=c(1,2))
      if (sum(Scenarios)==3){
        #MyBreaks<-seq(min(c(MyREC1[,varn],MyREC2[,varn]),na.rm=T),max(c(MyREC1[,varn],MyREC2[,varn]),na.rm=T),length.out=10)
        MyBreaks<-unique(as.vector(quantile(c(MyREC1[,varn],MyREC2[,varn]), probs = seq(0, 1, by = 1/15), na.rm = TRUE)))
        MapRivers(RECvar=MyREC1[,varn], REC=MyREC1, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 1",sep=" "), MyBreaks=MyBreaks, col.ord=coldir[Ind]) 
        MapRivers(RECvar=MyREC2[,varn], REC=MyREC2, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 2",sep=" "),  MyBreaks=MyBreaks,col.ord=coldir[Ind]) 
        
      }else{
        if (sum(Scenarios)==4){
          #MyBreaks<-seq(min(c(MyREC1[,varn],MyREC3[,varn]),na.rm=T),max(c(MyREC1[,varn],MyREC3[,varn]),na.rm=T),length.out=10)
          MyBreaks<-unique(as.vector(quantile(c(MyREC1[,varn],MyREC3[,varn]), probs = seq(0, 1, by = 1/15), na.rm = TRUE)))
          MapRivers(RECvar=MyREC1[,varn], REC=MyREC1, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 1",sep=" "), MyBreaks=MyBreaks, col.ord=coldir[Ind]) 
          MapRivers(RECvar=MyREC3[,varn], REC=MyREC3,MinOrder=MinOrder, pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 3",sep=" "),  MyBreaks=MyBreaks,col.ord=coldir[Ind]) 
        }else{
          #MyBreaks<-seq(min(c(MyREC3[,varn],MyREC2[,varn]),na.rm=T),max(c(MyREC3[,varn],MyREC2[,varn]),na.rm=T),length.out=10)
          MyBreaks<-unique(as.vector(quantile(c(MyREC2[,varn],MyREC3[,varn]), probs = seq(0, 1, by = 1/15), na.rm = TRUE)))
          MapRivers(RECvar=MyREC2[,varn], REC=MyREC2, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 2",sep=" "), MyBreaks=MyBreaks, col.ord=coldir[Ind]) 
          MapRivers(RECvar=MyREC3[,varn], REC=MyREC3, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 3",sep=" "),  MyBreaks=MyBreaks,col.ord=coldir[Ind])                     
        }
      }
      
    }else{
      x11(height=h, width=20) 
      if(sum(Scenarios)==1){
        MapRivers(RECvar=MyREC1[,varn], REC=MyREC1, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 1",sep=" "), n.breaks=10, col.ord=coldir[Ind])
      }else{
        if (sum(Scenarios)==2){
          MapRivers(RECvar=MyREC2[,varn], REC=MyREC2, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 2",sep=" "), n.breaks=10, col.ord=coldir[Ind])
        }else{
          MapRivers(RECvar=MyREC3[,varn], REC=MyREC3, MinOrder=MinOrder,pos="topleft", myCol = "lightblue", myScale = 0.5, name = legname, main=paste(varname,"Scenario 3",sep=" "), n.breaks=10, col.ord=coldir[Ind])
        }
      } 
    }      
  }
  

}

###################################################################################################################
####################################################################################################################

addLine<- function(a=NULL, b=NULL, v = NULL, h = NULL, ..., once=F) {
  #This function adds a line to all trellis plots
  tcL <- trellis.currentLayout()
  k<-0
  for(i in 1:nrow(tcL))
    for(j in 1:ncol(tcL))
      if (tcL[i,j] > 0) {
         k<-k+1
                trellis.focus("panel", j, i, highlight = FALSE)
         if (once) panel.abline(a=a[k], b=b[k], v=v[k], h=h[k], ...) else
           panel.abline(a=a, b=b, v=v, h=h, ...)
                trellis.unfocus()
                }
   }

###################################################################################################################
####################################################################################################################
ChooseCatch<-function(MyREC=NULL,Terms=1,LgCatch=1){
  if (length(MyREC1[,1])>10000){MinOrder<-3}else{MinOrder<-min(MyREC1$ORDER)}
  x11();ScLeg<-MapRivers(RECvar=MyREC$MeanFlow*1000, REC=MyREC,MinOrder=MinOrder, pos="topleft", name="Mean Flow (l/s)", n.breaks=10, col.ord=F) 
if (Terms==1){
TerminalReach<-match(setdiff(MyREC$tnode,MyREC$fnode),MyREC$tnode)
TerminalREC <- MyREC[TerminalReach, ] #[MyREC$ORDER > 6]only pick out the largest catchments
points(TerminalREC[TerminalREC$ORDER >5, c("segXcentroid", "segYcentroid")], pch=15) #     TerminalReach
print("Click with pointer on map to locate catchment outlet of interest...")
  if (LgCatch==1){TerminalReach<-TerminalReach[TerminalREC$ORDER>5] } #restricts selection so that only the large catchments can be selected
}else{
  TerminalReach<-seq(1,length(pick),1)
  TerminalREC <- MyREC
  points(TerminalREC[TerminalREC$ORDER >5, c("segXcentroid", "segYcentroid")], pch=15) #     TerminalReach
  print("Click with pointer on map to locate site of interest...")
}


Point <- locator(1)

names(Point) <- c("segXcentroid", "segYcentroid") 
TheDists <- as.matrix(dist(rbind(Point, MyREC[TerminalReach, c("segXcentroid", "segYcentroid")])))[-1, 1]     
EndR <- as.numeric(names(TheDists[match(min(TheDists), TheDists)]))# the likely NZReach located
#setwd(dir$fun)
#source('FindCatchment.R')

#Call the function to find the catchment reaches
CatchSel<-c(EndR,FindCatchment(EndR,MyREC))
#Create a new REC list containing only upstream reaches
MyREC_Catch <- MyREC[match(CatchSel, MyREC$NZReach), ]
match(c("fnode", "tnode"), names(MyREC_Catch))

#Create a new REC list containing only upstream reaches
MyREC_Catch <- MyREC[match(CatchSel, MyREC$NZReach), ]

MyREC_Catch$tnode <- MyREC$tnode[match(CatchSel, MyREC$NZReach)]
MyREC_Catch$fnode <- MyREC$fnode[match(CatchSel, MyREC$NZReach)]
x11()
MapRivers(RECvar=MyREC_Catch$MeanFlow*1000, REC=MyREC_Catch,  pos="topright", myCol = "lightblue", MyLegend = ScLeg, name="Mean Flow (l/s)", n.breaks=10, col.ord=F,PlotTitle="Selected catchment")

return(Point)

}

###################################################################################################################
####################################################################################################################

MultiWheel<-function(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,Parlims=MyParlims,INDICATORsel=INDICATORsel,Point=MyPoint,MinQ=MinQ,AllocQ=AllocQ,LandManagement=LandManagement,UPSTind=c(1,0.05),WheelLayout=c(1,4,0)){
  #setwd(dir$fun)
  source(paste(dir$fun,"DrawWaterWheel.R",sep="/"))
    
  TerminalReach<-match(setdiff(MyREC1$tnode,MyREC1$fnode),MyREC1$tnode)
  TheDists <- as.matrix(dist(rbind(Point, MyREC1[TerminalReach, c("segXcentroid", "segYcentroid")])))[-1, 1]     
  MyPointNZReach <- names(TheDists[match(min(TheDists), TheDists)])# the likely NZReach located
  if (UPSTind[1]==1){
    MyOutComes<-Integrate.Upst.Outcomes(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn=varn,p=UPSTind[2],MyPointNZReach=MyPointNZReach)
    sub2<-"Catchment WQ estimates"
  }else{
  MyOutComes <- rbind(MyREC1[MyPointNZReach,],MyREC2[MyPointNZReach,],MyREC3[MyPointNZReach,])
  sub2<-"Segment WQ estimates"
  }
  MyOutComes<-MyOutComes[,varn[INDICATORsel]]
  x11(height=WheelLayout[1]*15, width=WheelLayout[2]*20);par(mfrow=WheelLayout[1:2], mar=c(1,1,3,1)) ## drw the wheels  2,1,2,1
  for (i in 1:3){
    PlotTitle <- paste("Water Management Scenario #", i)
    Parval<-MyOutComes[i,]
    if (WheelLayout[3]==1){
      names(Parval)<-labeln[INDICATORsel]
    }else{
      names(Parval)<-seq(1,length(Parval),1)
    }
    DrawWaterWheel(Parval,Parlim=MyParlims[INDICATORsel],PlotTitle,MinQ=NULL,AllocQ=AllocQ[i],LandManagement=LandManagement[i],sub2=sub2)
  } 
  
  if (WheelLayout[3]==0){
    par(mar=c(4,4,4,4))
    plot(c(0,0,1,1,0),c(0,1,1,0,0),type='l',col='grey75',axes=FALSE,xlab="",ylab="",lwd=2)
  key<-"Key:\n"
  for (i in 1:length(Parval)){
    key<-paste(key,i,". ",labeln[INDICATORsel[i]],"\n")
  }
  text(0.1,0.5,key,adj=c(0,0.5))
  }
  return(MyOutComes)
}

Integrate.Upst.Outcomes<-function(MyREC1=MyREC1,MyREC2=MyREC2,MyREC3=MyREC3,varn=varn,p=0.05,MyPointNZReach=MyPointNZReach){
   #Start with the local estimates
  MyOutComes <- rbind(MyREC1[MyPointNZReach,varn],MyREC2[MyPointNZReach,varn],MyREC3[MyPointNZReach,varn])
  #Then replace those values that should be integrated
  sel<-c("CLAR","MCI","RedWidth","MeanFilsDeltaQ","MaxFilsDeltaQ","MeanMatsDeltaQ","MaxMatsDeltaQ","deltaEEL","deltaShortEEL","deltaTrout","deltaBully")
  ptl<-c(p,p,p,1-p,1-p,1-p,1-p,p,p,p,p) #Some indicators go in opposite direcitons so set limit accordingly
  
  #sel<-c("CLAR","MCI","RedWidth","MeanFilsDeltaQ","MaxFilsDeltaQ","MeanMatsDeltaQ","MaxMatsDeltaQ","deltaEEL","deltaShortEEL","deltaTrout")
  #ptl<-c(p,p,p,1-p,1-p,1-p,1-p,p,p)
  names(ptl)<-sel
  MyOutComes[1,sel]<-(unlist(lapply(sel,function(x) quantile(MyREC1[,x],ptl[x],na.rm = TRUE))))
  MyOutComes[2,sel]<-(unlist(lapply(sel,function(x) quantile(MyREC2[,x],ptl[x],na.rm = TRUE))))
  MyOutComes[3,sel]<-(unlist(lapply(sel,function(x) quantile(MyREC3[,x],ptl[x],na.rm = TRUE))))
  return(MyOutComes)
  
}

Integrate.Upst.Outcomes.K<-function(REC_OUT=REC_OUT,sc_sel=sc_sel,varn=varn,p=0.05,MyPointNZReach=MyPointNZReach,MinQ=MinQ,AllocQ=AllocQ,Swit=0,CatchSel=0){
  #Select only the rows that apply
  if (length(CatchSel)>1){
    ind<-match(CatchSel,REC_OUT[[1]]$NZReach)
    for (i in 1:length(REC_OUT))
    REC_OUT[[i]]<-REC_OUT[[i]][ind,]
  }
  
  #Start with the local estimates
  d_sel<-list()
  
  for (i in 1:length(sc_sel)){
  if (Swit==1){
    d_sel[[i]]<-which(REC_OUT[[sc_sel[i]]][,"TheTake"]>0)
  }else{
    d_sel[[i]]<-seq(1:length(REC_OUT[[1]][,1]))
  }}
  
  #MyOutComes<-colMeans(REC_OUT[[sc_sel[1]]][d_sel[[1]],varn],na.rm=T)
  #for (i in 2:length(sc_sel)){
  #MyOutComes<-rbind(MyOutComes,colMeans(REC_OUT[[sc_sel[i]]][d_sel[[i]],varn],na.rm=T))
  #}
  temp<-sapply(1:length(varn),function(x) median(REC_OUT[[sc_sel[1]]][d_sel[[1]],varn[x]],na.rm=T))
  dim(temp)<-c(1,length(varn))
  MyOutComes<- temp  
  
  for (i in 2:length(sc_sel)){
  temp<-sapply(1:length(varn),function(x) median(REC_OUT[[sc_sel[i]]][d_sel[[i]],varn[x]],na.rm=T))
  dim(temp)<-c(1,length(varn))  
  MyOutComes<-rbind(MyOutComes,temp)
  }
  MyOutComes<-as.data.frame(MyOutComes)
  names(MyOutComes)<-varn
  
 TerminalReach<-match(setdiff(REC_OUT[[1]]$tnode,REC_OUT[[1]]$fnode),REC_OUT[[1]]$tnode)
 #TerminalReach<-c(TerminalReach,which(REC_OUT[[1]][,"MALF"]==max(REC_OUT[[1]][,"MALF"])))
   #for (i in 1:length(sc_sel)){
  #MyOutComes[i,"IrriAreaIrrigated"]<-sum(REC_OUT[[sc_sel[i]]][TerminalReach,"usIrriArea"]*REC_OUT[[sc_sel[i]]][TerminalReach,"IrriAreaIrrigated"])/10000
  #}
#   qqq=0
#   if (qqq==1){
#     setwd(homedir)
#     #browser()
#   RelData<-read.csv("Kaituna_Production.csv",as.is = T, na.strings = "NaN",header=T)
#   for (i in 1:length(sc_sel)){
#     id<-which(RelData$MinQ==MinQ[sc_sel[i]]&RelData$AllocQ==AllocQ[sc_sel[i]])
#     MyOutComes[i,"freq.reliability"]<-RelData$ProdLoss[id]*100
#   }
#   }
  #Then replace those values that should be integrated
  sel<-c("CLAR","MCI","RedWidth","MeanFilsDeltaQ","MaxFilsDeltaQ","MeanMatsDeltaQ","MaxMatsDeltaQ","deltaEEL","deltaShortEEL","deltaTrout","deltaBully","deltaInanga","deltaTorrent","deltaKokopu")
  ptl<-c(p,p,p,1-p,1-p,1-p,1-p,p,p,p,p,p,p,p) #Some indicators go in opposite direcitons so set limit accordingly
  names(ptl)<-sel
  for (i in 1:length(sc_sel)){ 
  MyOutComes[i,sel]<-(unlist(lapply(sel,function(x) quantile(REC_OUT[[sc_sel[i]]][d_sel[[i]],x],ptl[x],na.rm = TRUE))))
  }
  return(MyOutComes)
  
}

MultiWheel_Kaituna<-function(REC_OUT=REC_OUT,sc_sel=sc_sel,Parlims=MyParlims,INDICATORsel=INDICATORsel,Point=MyPoint,MinQ=MinQ,AllocQ=AllocQ,LandManagement=LandManagement,UPSTind=c(1,0.05),WheelLayout=c(1,4,0),CatchSel=0,p_on=1){
 # setwd(dir$fun)
  source(paste(dir$fun,"DrawWaterWheel.R",sep="/"))

  if (UPSTind[1]==1){
    MyOutComes<-Integrate.Upst.Outcomes.K(REC_OUT=REC_OUT,sc_sel=sc_sel,varn=varn,p=UPSTind[2],MyPointNZReach=MyPointNZReach,MinQ=MinQ,AllocQ=AllocQ,Swit=UPSTind[3],CatchSel=CatchSel)
  }else{
 
    MyPointNZReach<-match(Point,REC_OUT[[1]]$NZReach)
    MyOutComes <- Doug.rbind.list(lapply(1:length(sc_sel),function(x) REC_OUT[[x]][MyPointNZReach,]))
  }
 
  if (is.na(match(1,INDICATORsel))==F){
    MyOutComes[,"TheTake"]<-MyOutComes[,"TheTake"]*1000
    labeln[1]<<-"Irrigation Take (l/s)"
  }
  
  MyOutComes<-MyOutComes[,varn[INDICATORsel]]
  if (p_on==1){
  x11(height=WheelLayout[1]*15, width=WheelLayout[2]*20);par(mfrow=WheelLayout[1:2], mar=c(1,1,3,1)) ## drw the wheels  2,1,2,1
  for (i in 1:length(sc_sel)){
    PlotTitle <- paste("Water Management Scenario #", sc_sel[i])
    Parval<-MyOutComes[sc_sel[i],]
    if (WheelLayout[3]==1){
      names(Parval)<-labeln[INDICATORsel]
    }else{
      names(Parval)<-seq(1,length(Parval),1)
    }
    DrawWaterWheel(Parval,Parlim=MyParlims[INDICATORsel],PlotTitle,MinQ=MinQ[sc_sel[i]],AllocQ=AllocQ[sc_sel[i]],LandManagement=NULL,sub2=NULL)
  } 
  
  
  
  if (WheelLayout[3]==0&WheelLayout[4]==0){
    par(mar=c(4,8,4,8))
    plot(c(0,0,1,1,0),c(0,1,1,0,0),type='l',col='grey75',axes=FALSE,xlab="",ylab="",lwd=2)
    key<-"Key:\n"
    for (i in 1:length(Parval)){
      key<-paste(key,i,". ",labeln[INDICATORsel[i]],"\n")
    }
    text(0.1,0.5,key,adj=c(0,0.5))
  }
  }
  return(MyOutComes)
  
}

################################################################################
# HEAT PLOTS  of the outcomes for all scenarios
################################################################################

NiceImage <- function(myMat, xlab=NULL,ylab=NULL,Reverse=F...) {
  par(mar=c(5,5,2,2), xpd = F)
  MyCols <- (heat.colors(12));  if(Reverse == T) MyCols <- rev(MyCols)
  image(myMat, x = 0:nrow(myMat), y = 0:ncol(myMat), axes = F, xlab=xlab,ylab=ylab,col = MyCols)
  abline(v = seq(0,nrow(myMat), length = nrow(myMat) + 1) )
  abline(h = seq(0,ncol(myMat), length = ncol(myMat) + 1) )
  par(xpd = T)
  myX <- seq(0.5, nrow(myMat)-0.5, by = 1)
  myY <- seq(0.5, ncol(myMat)-0.5, by = 1)
  text(x = myX, y = 0, rownames(myMat), pos = 1)
  text(x = 0, y = myY, colnames(myMat), pos = 2)# , srt = 90
  for(n in 1:nrow(myMat)) {
    for(m in 1:ncol(myMat)) {
      text(n-0.5, m-0.4, round(myMat[n,m], digits = 0))
      #text(n-0.5, m-0.6, paste("(", MyOtherMat[n,m], ")", sep=""))
    }
  }
  return(NULL)
}


MultiImage<-function(INDICATORsel,MyOutComes,ManageScenarios){
  #browser()
  r<-ceiling(sqrt(length(INDICATORsel)))
  x11();par(mfrow=c(ceiling(length(INDICATORsel)/r),r), mar=c(1,1,3,1)) 
  for (i in 1:length(INDICATORsel)){
    heat_sel<-INDICATORsel[i]
    MyMat<-MyOutComes[,varn[heat_sel]]
    dim(MyMat)<-c(5,6)
    rownames(MyMat)<-unique(ManageScenarios[,1])
    colnames(MyMat)<-unique(ManageScenarios[,2])
    NiceImage(MyMat,ylab="Minimum Flow (multiple of MALF)",xlab="Allocation (multiple  of MALF)",Reverse=coldir[INDICATORsel[i]]);title(main=labeln[heat_sel])
    #NiceImage(MyMat,ylab="Minimum Flow (proportion of MALF)",xlab="Allocation (proportion of MALF)");title(main="Irrigated Area (ha)")
  }
}

MultiContour<-function(INDICATORsel,MyOutComes,ManageScenarios){
  #browser()
  r<-ceiling(sqrt(length(INDICATORsel)))
  x11();par(mfrow=c(ceiling(length(INDICATORsel)/r),r), mar=c(1,1,3,1)) 
  for (i in 1:length(INDICATORsel)){
    heat_sel<-INDICATORsel[i]
    MyMat<-MyOutComes[,varn[heat_sel]]
    dim(MyMat)<-c(5,6)
    rownames(MyMat)<-unique(ManageScenarios[,1])
    colnames(MyMat)<-unique(ManageScenarios[,2])
    NiceContour(MyMat,ManageScenarios,ylab="Minimum Flow (proportion of MALF)",xlab="Allocation (proportion of MALF)",Reverse=coldir[INDICATORsel[i]]);title(main=labeln[heat_sel])
    #NiceImage(MyMat,ylab="Minimum Flow (proportion of MALF)",xlab="Allocation (proportion of MALF)");title(main="Irrigated Area (ha)")
  }
}

NiceContour <- function(myMat, ManageScenarios,xlab=NULL,ylab=NULL,Reverse=F...) {
  par(mar=c(5,5,2,2), xpd = F)
  #MyCols <- (heat.colors(12));  if(Reverse == T) MyCols <- rev(MyCols)
 browser()
  contour(z=myMat, x = ManageScenarios[,1], y = ManageScenarios[,2], axes = F, xlab=xlab,ylab=ylab)
  myX <- seq(0.5, nrow(myMat)-0.5, by = 1)
  myY <- seq(0.5, ncol(myMat)-0.5, by = 1)
  text(x = myX, y = 0, rownames(myMat), pos = 1)
  text(x = 0, y = myY, colnames(myMat), pos = 2)# , srt = 90
  return(NULL)
}
