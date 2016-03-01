#########
# function to rbind all the elements in a list of data frames
# works even when all NA's
Doug.rbind.list <- function (dfs)
{
  # dfs <- list(...)
  if (length(dfs) == 0)
    return(list())
  all.names <- unique(unlist(lapply(dfs, names)))
  data.frame(do.call("rbind", lapply(dfs, function(df) df)))
}
#########

#########
DougSplit <- function(mytext, myPart = 1, ...) {
  myOut <- strsplit(mytext, ...)[[1]][myPart]
  return(myOut)
}
#########

# function to cbind all the elements in a list of data frames
# works even when all NA's
Doug.cbind.list <- function (dfs)
{
  data.frame(do.call("cbind", lapply(dfs, function(df) df)))
}
#########
#########

# REMOVE NA values in colunms - replace with means!!
rep.na <- function(df) {
  for (i in 1:ncol(df)) {
    if (sum(as.numeric(is.na(df[,i])))>0) {
      df[,i][is.na(df[,i])] <- mean(df[,i], na.rm = TRUE)
      cat("inserting mean into col", i, "\n")
    } else {
      cat("no NA values in column", i, "\n")
    }
  }
  df
} # end function


###########################################################################
#Does a Principal comonemtns analysis and a nice plot
PC_Plot<-function(acp=NULL,groups=NULL,ptype="text",PlotTitle="PCA"){
  N = length(acp$scores[,1]) 
  if (N>200){psc<-0.5}else{psc<-2}
  if (is.null(groups)==T){
    groups = matrix(1,nrow=N,ncol=1)
  }
  
  # The original biplot for comparison 
  #biplot(acp) 
  # Compute de min/max of new coordinates 
  xmin = min(acp$scores[,1])
  xmax = max(acp$scores[,1])  
  ymin = min(acp$scores[,2]) 
  ymax = max(acp$scores[,2])
  x_extra<-(xmax-xmin)*0.3
  y_extra<-(ymax-ymin)*0.3
  plot(c(xmin-x_extra,xmax+x_extra),c(ymin-y_extra,ymax+y_extra),col="white", xlab="Comp 1", ylab="Comp 2") 
  # Plot the points with colors 
  if (ptype=="text"){
    text(acp$scores[,1],acp$scores[,2],1:N,col= groups,cex=0.5)
  }else{
    points(acp$scores[,1],acp$scores[,2],col= groups,pch=20,cex=psc)
    legend("topleft", legend=sort(unique(groups)), title="Class", col = sort(unique(groups)),
           pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)   
  }  
  
  title(PlotTitle) 
  abline(v=0,lty=2) 
  abline(h=0, lty=2) 
  # Compute and apply a re-scaling for the arrows of old components 
  # xl.min for min of xloadings,... 
  xl.min = min(0,min(acp$loadings[,1])) 
  xl.max = max(0,max(acp$loadings[,1])) 
  yl.min = min(0,min(acp$loadings[,2])) 
  yl.max = max(0,max(acp$loadings[,2])) 
  xl.scale = max(abs(xmax),abs(xmin))/max(abs(xl.max),abs(xl.min))*0.75 
  # 0.75 factor is just for leave some place for the text 
  yl.scale = max(abs(ymax),abs(ymin))/max(abs(yl.max),abs(yl.min))*0.75 
  #Draw old components 
  arrows(rep(0,100),rep(0,100),acp$loadings[,1]*xl.scale, 
         acp$loadings[,2]*yl.scale,col="red") 
  # Names of old components 
  text(acp$loadings[,1]*xl.scale*1.15, acp$loadings[,2]*yl.scale*1.15 
       ,rownames(acp$loadings),col="red",cex=0.7) 
}
###########################################################################

###########################################################################
#GENERAL FUNCITON FOR DOING CLUSTER ANALYSIS
doclust=function(OBS=NULL,Pred=NULL,plotpca=F,doenchain=F,levplot=3,MyREC_ALL=NULL,MyREC_ALL2=NULL,minord=1,rst=F){
  OUTPUTS<-list()
  PCA_obs<-princomp(x=OBS,cor=TRUE)
  if (plotpca==T){x11();PC_Plot(PCA_obs)}
  OUTPUTS[["pVarexpl_obs"]]<-PCA_obs$sdev^2/sum(PCA_obs$sdev^2)    
  
  PCA_mod<-princomp(x=Pred,cor=TRUE)
  if (plotpca==T){x11();PC_Plot(acp=PCA_mod,ptype="points")}
  OUTPUTS[["pVarexpl_mod"]]<-PCA_mod$sdev^2/sum(PCA_mod$sdev^2)    
  ###############################################################################################
  if (rst==F){
    WQClass_mod<-classify.domain(x=PCA_mod$scores, MyMethod="manhattan", clusters="ward", level=seq(2,20,2)) 
    if (doenchain==T){
      sel<-sample(length(Pred[,1]),500)
      test_mod<-enchain(PCA_mod$scores[sel,],classes=WQClass_mod[sel,],method="euclidean",
                        NumClass=seq(2,20,2),stat=anosift,mins=10,permutations=F,nboot=10)
      x11();plot.enchain(test_mod)
    }
    #browser()
    x11(height=40, width=40) ; par(mfrow=c(2,2),mar=c(5,5,5,3))
    MapRivers_cat(classes=WQClass_mod,level=levplot, MinOrder = minord, REC = MyREC_ALL[match(row.names(MyREC_ALL2),MyREC_ALL$NZReach),],
                  VarTable=PCA_mod$scores,VarSel=1:3,PlotTitle="Class",col.rand=T)
    MapRivers_cat(classes=WQClass_mod,level=levplot, MinOrder = minord, REC = MyREC_ALL[match(row.names(MyREC_ALL2),MyREC_ALL$NZReach),],
                  VarTable=PCA_mod$scores,VarSel=1:3,PlotTitle="Class")
  }else{
    PCA_moda<-decostand(PCA_mod$scores,method="range")
    WQClass_mod<-classify.domain(x=PCA_moda, MyMethod="manhattan", clusters="ward", level=seq(2,20,2)) 
    if (doenchain==T){
      sel<-sample(length(Pred[,1]),500)
      test_moda<-enchain(PCA_moda[sel,],classes=WQClass_moda[sel,],method="euclidean",
                         NumClass=seq(2,20,2),stat=anosift,mins=10,permutations=F,nboot=10)
      x11();plot.enchain(test_moda)
    }
    
    x11(height=40, width=40) ; par(mfrow=c(2,2),mar=c(5,5,5,3))
    MapRivers_cat(classes=WQClass_mod,level=levplot, MinOrder = minord, REC = MyREC_ALL[match(row.names(MyREC_ALL2),MyREC_ALL$NZReach),],
                  VarTable=PCA_moda,VarSel=1:3,PlotTitle="Class",col.rand=T)
    MapRivers_cat(classes=WQClass_mod,level=levplot, MinOrder = minord, REC = MyREC_ALL[match(row.names(MyREC_ALL2),MyREC_ALL$NZReach),],
                  VarTable=PCA_moda,VarSel=1:3,PlotTitle="Class")
  }
  par(mar=c(5,5,5,3))
  PC_Plot(acp=PCA_mod,ptype="points",groups=WQClass_mod[,levplot])
  PC_Plot(acp=PCA_obs,ptype="points",groups=WQClass_mod[match(Data$NZReach,row.names(Pred)),levplot])    
  
  
  GP_summary_mod<-apply(Pred,2,function(x) by(x,WQClass_mod[,3],median))
  GP_summary_mod[,-agrep("SS",names(Pred))] <-10^GP_summary_mod[,-agrep("SS",names(Pred))]  
  OUTPUTS[["ModSummary"]]<-GP_summary_mod
  return(OUTPUTS)
}
###############################################################################################
