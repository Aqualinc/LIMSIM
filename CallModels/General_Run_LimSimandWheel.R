###########################################################################################
#                                                                                         #
#        File to run analysis of change impacts                                           #
#                                                                                         #
###########################################################################################
rm(list=ls()) # clear memory
require(lattice)
require(latticeExtra)
require(randomForest)
##############################################################
#                   Define directories                       #
##############################################################

rootdir<-"H:/GitRepositories/LimSim"
rootdir<-"D:\\Projects\\Aqualinc\\projects\\MPI_WaterStorage\\LimSim" #Tim's computer
setwd(rootdir)
source("SetDirectories.r")
dir<-SetDirectories(rootdir)
#dir$proj<-"H:/GitRepositories/LimSim/CallModel/Datasets"
##############################################################
#                    Load functions                          #
##############################################################
source(file.path(dir$fun,"GeneralFunctions.r"))
source(file.path(dir$fun,"map_REC7.r"))
source(file.path(dir$fun,"LIMSIM Functions.R"))
source(file.path(dir$fun,"Flow&AbstractionFunctions.R"))
source(file.path(dir$fun,"WaterQ&PeriphytonFunctions.R"))
source(file.path(dir$fun,"FindCatchment.R"))
source(file.path(dir$fun,"HabitatFunctions.R"))

##############################################################
#            assemble data for a domain                      #
##############################################################
load(file.path(dir$data,"RECdata_9.RData"))

#Open the FWENZ data from the data directory. If it doesn't exist, stop and explain the problem.
FWENZDataFileName <- "RECandFWENZ.RData"
if (!file.exists(file.path(dir$data,FWENZDataFileName))){stop("File ",file.path(dir$data,FWENZDataFileName), " is needed and doesn't exist")}
load(file.path(dir$data,FWENZDataFileName))

#Open the Groundwater data from the data directory
GWDataFileName <- "RECGroundwater.RData"
if (!file.exists(file.path(dir$data,GWDataFileName))){stop("File ",file.path(dir$data,GWDataFileName), " is needed and doesn't exist")}
load(file.path(dir$data,GWDataFileName))

MyREC <- merge(MyREC,RECGroundwater)

##############################################################
#            select subset of data                           #
##############################################################
minord<-1             #Select minimum order to be included in the analysis
SelectSwitch<-1       #Select whethe to define by region (0) or a specified pt (1)
if (SelectSwitch==0){
      # EITHER......................................................
      # Choose region from:
      # 1.  Northland               9.  Greater Wellington
      # 2.  Auckland                10. Tasman
      # 3.  Waikato                 11. Marlborough 
      # 4.  Bay of Plenty           12. West Coast
      # 5.  East Cape               13. Canterbury
      # 6.  Taranaki                14. Otag0
      # 7.  Manawatu                15. Southland
      # 8.  Hawke's Bay             16. ALL NZ
      
      #And seclect a minimum stream order (from 1 to 8) (Recommend not having smaller than 4 for whole country)
      ##############################################################
      pick <- MakePick(MyREC=MyREC,MinOrder=3,Region=16)
      }else{
      # OR.........................................................
      # select by a downstream locaiton
      
      Pt<-11027203  #This is the outlet of Flaxbourne at Corrie Downs in Marlborough District
      pick<-c(Pt,FindCatchment(Pt,MyREC,minord))}

##############################################################
MyREC<-MyREC[match(pick, MyREC$NZReach),]

##############################################################
#    Set Scenario Characteristics                            #
##############################################################

#_____________________________________________________________
#DEFINING ALLOCATION SCENARIOS

#Defined as absolute values
MinQ_siteRECID  <- 11027203
MinQ_siteMALF   <- MyREC$MALF[MyREC$NZReach == MinQ_siteRECID]
MinQ_absolute   <- 0.025    #'A' block minimum
AllocQ_absolute <- 0.00226  #'A' block allocation

#Transform back into values relative to MALF
MinQ_sel        <- MinQ_absolute / MinQ_siteMALF
AllocQ_sel      <- AllocQ_absolute / MinQ_siteMALF

#MinQ_sel   <-as.double(c(0.7))      #Define minimum flows (as a proportion of MALF)
#AllocQ_sel <-as.double(c(0.5))      #Define Allocation volumes (as a proportion of MALF)
#GWAlloc    <-as.double(c(0.0,0.2,0.4,0.6,0.8,1.0))      #Define groundwater allocation as a fraction of mean annual recharge
GWAlloc    <-as.double(c(0))
#____________________________________________________________________________________________
#USER DEFINED SETTINGS FOR ALL SCENARIOS
TakeAll             <-0                 #Define whether to take ALL allocated water (1) or just that for irrigable area (0)
IrrigableAreaTarget <-100               #Proportion of irrigable area we will aim to irrigate (0:100)
PropIrri            <-1                 #Binary switch to indicate whether to scale land management factors by  proportion of upstream pasture that is irrigated
SysCap              <-0.58              #System capacity (l/s) default to 0.58 l/s (5 mm/d)
EffIrriArea         <-0.8               #Effective irrigable area - scale irrigable area to get actual irrigated area (0:100)
WQModel             <-"Unwin"           #Select which WaterQuality model to use Unwin or McDowell (default is Unwin)
QWModel             <-"Jowett"          #Select which FlowWidth model to use: Jowett or Booker (default is Jowett)
#############################################################
#  POPULATE SCENARIOS

#Function to run "LIMSIM" for one scenario
RunOne=function(MyREC=MyREC,AllocQ=AllocQ,MinQ=MinQ,GWAlloc=GWAlloc,pick=pick,LandManagement="Fair",IrrigableAreaTarget=100,PropIrri=1,
                TakeAll=TakeAll,SysCap=SysCap,EffIrriArea=EffIrriArea,WQModel=WQModel,QWModel=QWModel){
  #browser()
  MyREC1<-Run_FlowandAbstraction(MyREC=MyREC,AllocQ=AllocQ,MinQ=MinQ,GWAlloc=GWAlloc,pick=pick,TakeAll=TakeAll,SysCap=SysCap,
                                 EffIrriArea=EffIrriArea)
  #The water quality and periphyton models, are affected by BFI, nNeg and FRE3.Count which are altered by the flow and abstraction routines above.
  MyREC1<-Run_WQandPeri(MyREC=MyREC1,LandManagement=LandManagement,IrrigableAreaTarget=IrrigableAreaTarget,
                        AllocQ=AllocQ,GWAlloc=GWAlloc,PropIrri=PropIrri,WQModel=WQModel) 
  MyREC1<-Run_HabitatFunctions(MyREC=MyREC1,pick=pick,MinQ=MinQ,method=QWModel)
  return(MyREC1)
}

ManageScenarios<-expand.grid(AllocQ_sel,MinQ_sel,GWAlloc)
names(ManageScenarios)<-c("AllocQ","MinQ","GWAlloc")

LandManagement<-c("Fair")               #Select a management type of: "Good", "Fair" or "Poor"
REC_OUT<-list()
#for (i in 1){      #for testing only
for (i in 1:length(ManageScenarios[,1])){
 #i<-1 
  print(paste("########## SCENARIO ",i," ############",sep=""))
  REC_OUT[[i]]<-RunOne(MyREC=MyREC,AllocQ=ManageScenarios[i,1],MinQ=ManageScenarios[i,2],GWAlloc=ManageScenarios[i,3],
                       pick=pick,LandManagement=LandManagement,IrrigableAreaTarget=IrrigableAreaTarget,
                       PropIrri=PropIrri,TakeAll=TakeAll,SysCap=SysCap,EffIrriArea=EffIrriArea,
                       WQModel=WQModel,QWModel=QWModel)
}

#OR Run just one scenario:
#REC_OUT<-RunOne(MyREC=MyREC,AllocQ=1.5,MinQ=0.5,GWAlloc=0,pick=pick,LandManagement=LandManagement,IrrigableAreaTarget=100,PropIrri=1,TakeAll=TakeAll)

#Or run for a single site
SingleSiteOutput<-RunOne(MyREC=MyREC[1,],AllocQ=0.7,MinQ=0.5,GWAlloc=0,
                       pick=pick[1],LandManagement=LandManagement,IrrigableAreaTarget=IrrigableAreaTarget,
                       PropIrri=PropIrri,TakeAll=TakeAll,SysCap=SysCap,EffIrriArea=EffIrriArea,
                       WQModel=WQModel,QWModel=QWModel)

#Or single site for just flow reliability
SingleSiteOutputReliability<-Run_FlowandAbstraction(MyREC=MyREC[1,],AllocQ=0.7,MinQ=0.5,GWAlloc=0,pick=pick[1],TakeAll=TakeAll,SysCap=SysCap,
                               EffIrriArea=EffIrriArea)

freq.restrict.multiband()

#OPTION TO SAVE FILES FOR REUSE LATER....
#setwd(dir$proj)
#save(REC_OUT,ManageScenarios,file="WairauScenarios_4.RData")
#load("WairauScenarios_4.RData")
#############################################################################
# Create plots of baseline information:
i<-5 #SElect a scenario
pos_sel="bottomright"
#Plot up stream irrigable area 
x11();MapRivers(RECvar=REC_OUT[[i]][,"usIrriArea"]/10000, REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "usIrriArea (ha)", main="Upstream Irrigable Area", MyBreaks=c(0,100,200,300,400,500,1000,5000,50000), col.ord=T)                           
#Plot upstream pasture
x11();MapRivers(RECvar=REC_OUT[[i]][,"usPastoral"]*REC_OUT[[i]][,"usArea"]/10000, REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "usPasture (ha)", main="Upstream Pasture", MyBreaks=c(0,100,200,300,400,500,1000,5000), col.ord=T)                            
#Plot Mean Annual Low Flow
x11();MapRivers(RECvar=REC_OUT[[i]][,"MALF"]*1000, REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "MALF l/s", main="Mean Annual Low Flow", n.breaks=10, col.ord=T)               
#Plot mean flow
x11();MapRivers(RECvar=REC_OUT[[i]][,"MeanFlow"], REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "Mean Flow m3/s", main="Mean Flow",n.breaks=10 , col.ord=T)               
#Plot T95
x11();MapRivers(RECvar=MyREC[,"MALF"], REC=MyREC, pos="topleft", myCol = "lightblue", myScale = 0.5, name = "T95", main="T95",n.breaks=10 , col.ord=T)               
#x11();MapRivers(RECvar=MyREC[,"T95"], REC=MyREC, pos="topleft", myCol = "lightblue", myScale = 0.5, name = "T95", main="T95",MyBreaks=c(0,12,14,16,18,20,22,24,26), col.ord=T,linevar=F)               


TerminalReach<-match(setdiff(REC_OUT[[1]]$tnode,REC_OUT[[1]]$fnode),REC_OUT[[1]]$tnode)#MyBreaks=c(1,10,100,500,1000,5000,10000,50000,100000)/1000
TerminalReach<-c(TerminalReach,which(REC_OUT[[1]][,"MALF"]==max(REC_OUT[[1]][,"MALF"])))
points(REC_OUT[[1]][TerminalReach,"segXcentroid"],REC_OUT[[1]][TerminalReach,"segYcentroid"])


pick2<-match(MyREC$NZReach,row.names(FWENZ))
MyFWENZ<-FWENZ[pick2,];rm(NewREC,FWENZ)
x11();MapRivers(RECvar=MyFWENZ$segAveTWarm/10, REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "segAvTWarm", main="Segment av. warm temperature",n.breaks=10 , col.ord=T)                           
x11()
MapRivers(RECvar=MyFWENZ$segAveTCold/10, REC=REC_OUT[[i]], pos=pos_sel, myCol = "lightblue", myScale = 0.5, name = "segAvTCold", main="Segment av. cold temperature",n.breaks=10 , col.ord=T)                           


##############################################################################
#            CREATE REGIONAL/NATIONAL PLOTS AND/OR MAPS                      #
##############################################################################

# Select indicators from this list (by number):
#   1.  The total abstraction for irrigation (m3/s)
#   2.  The percentage of time abstraction is stopped due to minimum flows (%)
#   3.  The percentage of time abstraciton is restricted due to management flows (%)
#   4.  The reliability of the take (% of total demand taken)
#   5.  The percentage of the upstream irrigable area irrigated\
#   6.  The fraction of the allocated water used
#   7.  Clarity
#   8.  MCI
#   9.  Reduction in River Width
#   10. Periphyton mean filaments
#   11. Periphyton Max annual Filaments
#   12. Periphyton mean mats
#   13. Periphyton max annual mats
#   14. Change in Long fin Eel Habitat
#   15. Change in Short Fin Eel Habitat
#   16. Change in Adult Brown trout Habitat
#   17. CHange in  blue gilled bullies
#   18. CHange in inanaga
#   19. CHange in torrent fish
#   20. CHange in Kokopu

#Seclect Plot type:
#   1.  Map (showing upto 3 maps for one indicator) (Specify which scenarios in Scenarios)
#   2.  Histogram of values for each scenario
#   3.  Density plots for each scenario by class (North/South, Upland/Lowland)###########NA FOR KAITUNA#############

###############################################################################
sel_sc<-c(1,2,3)

MakePlots(MyREC1=REC_OUT[[sel_sc[1]]],MyREC2=REC_OUT[[sel_sc[2]]],MyREC3=REC_OUT[[sel_sc[3]]],Indicator=4,PlotType=1,Scenarios=sel_sc)

##############################################################################
#                   PERFORM CATCHMENT ANALYSIS                               #
##############################################################################
#CHoose limits for the same 20 indicators listed eariler
MyParlims<-list()
MyParlims[[1]]<-c(0,1,2,6,8)            #   1.   The total abstraction for irrigation (m3/s)
MyParlims[[2]]<-c(100,30,20,10,0)       #   2.   The percentage of time abstraction is stopped due to minimum flows (%)
MyParlims[[3]]<-c(100,35,15,10,0)       #   3.   The percentage of time abstraciton is restricted due to management flows (%)
MyParlims[[4]]<-c(0,85,90,95,100)       #   4.   The reliability of the take (% of total demand taken)
MyParlims[[5]]<-c(-1,0,0.5,0.8,1)*100   #   5.   The fraction of the upstream irrigable area irrigated\
MyParlims[[6]]<-c(0,0.5,0.7,0.8,1)      #   6.   The fraction of the allocated water used
MyParlims[[7]]<-c(0.2,1.0,1.6,2.5,4.5)  #   7.   Clarity
MyParlims[[8]]<-c(70,86,105,106,110)    #   8.   MCI
MyParlims[[9]]<-c(-50, -40, -20, -10,0) #   9.   Reduction in River Width
MyParlims[[10]]<-c(50, 30, 15, 10, 0)   #   10.  Periphyton mean filaments
MyParlims[[11]]<-c(50, 30, 20, 10, 0)   #   11.  Periphyton max filaments
MyParlims[[12]]<-c(30, 15, 10, 5, 0)    #   12.  Periphyton mean mats
MyParlims[[13]]<-c(60, 40, 30, 10, 0)   #   13.  Periphyton max mats
MyParlims[[14]]<-c(0, 80, 85, 95, 250)  #   14.  Percentage of baseline habitat (Long fin eel) maintained
MyParlims[[15]]<-c(0, 80, 85, 95, 250)  #   15.  Percentage of baseline habitat (short fine eel) maintaine
MyParlims[[16]]<-c(0, 80, 85, 95, 250)  #   16.  Percentage of baseline habitat (brown trout) maintained
MyParlims[[17]]<-c(0, 90, 95, 120, 250) #   17.  Percentage of baseline habitat (BlueGillBully) maintained
MyParlims[[18]]<-c(0, 90, 95, 120, 250) #   18.  Percentage of baseline habitat (inanaga) maintained
MyParlims[[19]]<-c(0, 90, 95, 120, 250) #   19.  Percentage of baseline habitat (torrent fish) maintained
MyParlims[[20]]<-c(0, 90, 95, 120, 250) #   20.  Percentage of baseline habitat (Kokopu) maintained
#Select the indices of the indicators to include in the water wheel
#INDICATORsel<-c(4,5,9,10,14,16) #USE THIS FOR THE TABLE OF ALL OF THE INDICATORS FOR ALL OF THE SCENARIOS 
INDICATORsel<-c(1,2,3,4,5,6)    #Just the reliability indicators

WheelLayout<-c(2,2,0,1)       #number of rows (1) and columns (2).  NOTE: must be at least (no. scenarios +1) positions
# and (3) is whether to include full text labels [1] or number and legend [0]
sc_sel<-c(1:8)   #USE THIS TO GET A TABLE OF ALL OF THE INDICATORS FOR ALL OF THE SCENARIOS AND TO MAKE INPUTS FOR DECISION SPACE DIAGRAMS

#Determine the subcatchments

  
UPSTind<-c(1,0.95,1)          #1st value is switch to decide whether to do local [0] or catchment analysis [1].  
                            #2nd value is only used for catchment analysis and is the accpetable fraction of exceedance
                            #Third value is a switch, to either calculate over whole set (0) or subset affected by change (1)
MyPoint<-Pt 
MyOutComes<-MultiWheel_Kaituna(REC_OUT=REC_OUT,sc_sel=sc_sel,Parlims=MyParlims,INDICATORsel=INDICATORsel,Point=MyPoint,MinQ=ManageScenarios$MinQ,AllocQ=ManageScenarios$AllocQ,LandManagement="Fair",UPSTind=UPSTind,WheelLayout=WheelLayout,CatchSel=pick,p_on=0)


##############################################################################
#MAKE DECISION SPACE DIAGRAMS
#  MultiImage(INDICATORsel,MyOutComes,ManageScenarios)


##############################################################################
#MAKE Summary wheels


#sc_sel <-c(5,6,7)
sc_sel<-c(1,4,5)
nm<-c("A","B","C")

MyOutComesSummary<-MyOutComes[,varn[INDICATORsel]]
WheelLayout<-WheelLayout<-c(1,3,0,0) 
x11(height=WheelLayout[1]*15, width=WheelLayout[2]*20);par(mfrow=WheelLayout[1:2], mar=c(1,1,3,1)) ## drw the wheels  2,1,2,1
for (i in 1:length(sc_sel)){
  PlotTitle <- paste("Water Management Scenario #", nm[i])
  Parval<-MyOutComesSummary[sc_sel[i],]
  print(Parval)
  if (WheelLayout[3]==1){
    names(Parval)<-labeln[INDICATORsel]
  }else{
    names(Parval)<-seq(1,length(Parval),1)
  }
  DrawWaterWheel(Parval,Parlim=MyParlims[INDICATORsel],PlotTitle,MinQ=NULL,AllocQ=AllocQ[sc_sel[i]],LandManagement=NULL,sub2=NULL)
}

lims<-sapply(INDICATORsel, function(x) MyParlims[[x]])
