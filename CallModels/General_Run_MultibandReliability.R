##############################################################
#        File to run analysis of change impacts              #
##############################################################
rm(list=ls()) # clear memory

#Define directories
{
  ProjectDirectory    <- "D:\\Projects\\Aqualinc\\projects\\MPI_WaterStorage\\LimSim"  #This is for Rainfall.NZ. Edit to match local directory structure
  DataDirectory       <- file.path(ProjectDirectory,"Datasets")                        #Assumes all data is in "Datasets" sub directory of the project directory
  RFunctionsDirectory <- file.path(ProjectDirectory, "Functions")                      #Assumes all additional R files are in a "Functions" sub directory of the project directory
  FlowAllocationFile  <- file.path(ProjectDirectory,"..\\FlaxbourneExample\\WaterTakeLimits_FlaxbourneExample.csv")
}

# Load functions
{
  source(file.path(RFunctionsDirectory,"Flow&AbstractionFunctions.R"))   #This file has the "freq.restrict.multiband" function in it, needed to calculate take reliability
  source(file.path(RFunctionsDirectory,"HabitatFunctions.R"))            #This file has the "IntegratedWidthV2" function in it, which calculates the change in width for two different FDC's
  
  }

# Load general data
{
  load(file.path(DataDirectory,"RECdata_9.RData")) #This is the "master" RECV1 data frame (called MyREC) complete with additional attributes including the flow duration curve parameters
  
  #load the Groundwater data from the data directory. This loads an data frame of Aquifer parameters for each REC reach.
  load(file.path(DataDirectory,"RECGroundwater.RData"))
  
  #Combine all the attributes together
  MyREC <- merge(MyREC,RECGroundwater)
}

# # Set site-specific data
# 
# SiteRECID           <- 11027203                           #The RECV1 reach ID of the site of interest
# MinQ                <- c(0.025,0.045,0.25,0.4,0.6)        #A vector of band minimum flows, in cumecs
# AllocQ              <- c(0.002,0.005,0,0.05,1)            #A vector of band allocation flows, in cumecs
# AllocShare          <- c(50,50,50,50,50)                  #A vector of allocation sharing percentages when flow is less than the band's minimum flow + allocation flow
# 
# TestData <- list(RECV1ReachID=SiteRECID,MinQ=MinQ_absolute,AllocQ=AllocQ_absolute,AllocShare=AllocationSharePct)
# 
# 
# AllSites <- list(TestData)
# names(AllSites) <- "Flaxbourne at Corrie Downs"

# Read in allocation file and convert to a list of relevant parameters for each site
{
  AllocationData <- read.csv(FlowAllocationFile)
  
  # Cut down to just the columns of interest
  AllocationData <- AllocationData[,c("SiteName","RECV1ReachID","MinQ","AllocQ","AllocShare")]
  
  #Split on SiteName into a list of data frames
  SplitDataFrames <- split(AllocationData[2:5],AllocationData$SiteName)
  
  #Convert each site's data frame into a list and remove the duplicate REC SiteID's
  AllSites <- lapply(SplitDataFrames, function(x) {
    x <- as.list(x)
    x$RECV1ReachID <- unique(x$RECV1ReachID)
    return(x)
  })
}

# Find the reliability for each monitoring site. 
Reliabilities <- lapply(seq_along(AllSites), function(SiteIndex){
  SiteName <- names(AllSites)[[SiteIndex]]
  EachSite <- AllSites[[SiteIndex]]
  SiteREC  <- MyREC[match(EachSite[["RECV1ReachID"]], MyREC$NZReach),]
  
  Reliability <- freq.restrict.multiband(ThisNZReach = EachSite[["RECV1ReachID"]],
                                         minFlow=EachSite[["MinQ"]],
                                         allocation=EachSite[["AllocQ"]],
                                         allocation_share=EachSite[["AllocShare"]],
                                         Data = SiteREC, FDCPlot = TRUE,SiteName = SiteName)
  return(Reliability)
})

#Convert the Reliabilities list to a data frame and add a column with the site name
ReliabilitiesData.Frame <- as.data.frame(do.call(rbind,Reliabilities))
ReliabilitiesData.Frame <- cbind("SiteName"=row.names(ReliabilitiesData.Frame),ReliabilitiesData.Frame)

# Save the output to a file
write.table(ReliabilitiesData.Frame, file = file.path(DataDirectory, "Reliabilities.csv"),sep=",",row.names = FALSE,quote = FALSE)

#Undertake change in river width and habitat calculations - a demonstration example. The resut is a one row data frame with all the extra attributes attached
ReachAttributesWithWidthandDeltaHabitats <- Run_HabitatFunctions(MyREC = MyREC,pick=11027203,MinQ=0.1,GWAlloc=0)
