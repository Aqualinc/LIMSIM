setwd(dir$data)
NZmap <- read.csv("NZmap.csv")  # add coastline
NZmap <- NZmap[NZmap$Section == "North Island" | NZmap$Section == "South Island", ]
setwd(dir$fun)
###############################################################################
#        Map any continuous variable for REC reaches  in n.breaks colours            #
#        AND drawing representaion of river lines 
#        the legend can be saved and reused on new maps!!                                                                    #
###############################################################################    


MapRivers_cat <- function(classes=NULL,level=1, MinOrder = 0, REC = REC.dat, pos="topleft",  
                      myScale = 0.5, name="Variable",PlotTitle=NULL,VarTable=NULL,VarSel=NULL,col.rand=F,...) {
  
  #browser()
  tclass <-  as.factor(classes[,level])
  n.class <- length(which(as.vector((table(tclass) >0))))
  Tclass <-  table(tclass) 
  OrdClass <- names(sort(Tclass, decreasing=T))  # put in decreasing order by size of class
  PropClass <- round(100*(Tclass/length(tclass)),1) # proportion of domain asigned to each class
  class.names <-  names(Tclass)
  n.cat<-length(class.names)
  
  RECvar<- classes[, level]
  TheseSites <- which(REC$ORDER > MinOrder&!is.na(RECvar)==T) 
  
  RECsample <- REC[TheseSites, ]
  REC_plot<-REC[which(REC$ORDER>MinOrder),]
  xy <-  cbind(REC$segXcentroid, REC$segYcentroid) 
  par(mar = c(0.1,0.1,2, 0.1))  
  plot(xy, type = "n", xlab=" ", ylab=" ", axes=FALSE, ...)
  polygon(x = NZmap$NZE, y = NZmap$NZN, col = "grey95", border = "grey80")

  RECvar <-  as.numeric(RECvar[TheseSites])
 # RECvar.col<-rgb(runif(n.cat),runif(n.cat),runif(n.cat))  #Raisample(1:n.cat, n.cat, replace=F)/n.cat
  if (is.null(VarTable)==T|col.rand==T){
  RECvar.col<-1:n.cat
  }else{
    #browser()
    rgb.control <- Doug.rbind.list(by(VarTable[TheseSites,VarSel], as.factor(RECvar), function(x) colMeans(x)))
    ClassRGBLoads <- decostand(rgb.control, "range") # range 0:1
    names(ClassRGBLoads) <- c("red", "green", "blue")
    RECvar.col <- rgb(ClassRGBLoads, max=1)
    names(RECvar.col) <- row.names(ClassRGBLoads)
  }
  
  
if(is.null(REC_plot$fnode)) {
 print("For faster plotting next time run function 'IDfNode' and save ouput to REC.dat as '$fnode'") 
  fnode<-REC_plot$segYcentroid*NA    # Identify the node
 tnode<-REC_plot$segYcentroid*NA    # Identify the node
  for (n in 1:length(fnode)){
    val<-which(REC_plot$Nztnode==REC_plot$Nzfnode[n])
    if (length(val)>=1) {
      if (length(val)>=1) {
          ind=which(REC_plot$ORDER[val]==max(REC_plot$ORDER[val]))
          fnode[n]<-REC_plot$NZReach[val[ind[1]]]
      }else{
      fnode[n]<-REC_plot$NZReach[val[1]]}}       
  val2<-which(REC_plot$Nztnode==REC_plot$Nztnode[n])
   if (length(val2)>=1) {
      if (length(val2)>=1) {
          ind=which(REC_plot$ORDER[val2]==max(REC_plot$ORDER[val2]))
          tnode[n]<-REC_plot$NZReach[val2[ind[1]]]
    }else{
      tnode[n]<-REC_plot$NZReach[val2[1]]}}
  }
 
} else {
fnode <- REC_plot$fnode
tnode<-REC_plot$tnode
}
  
  fnodes<-fnode
  tnodes<-tnode
  fnode<-fnode[match(RECsample$NZReach,REC_plot$NZReach)]
  tnode<-tnode[match(RECsample$NZReach,REC_plot$NZReach)]
  
  #Convert fnode and tnode back into row indexes
  fnode<-unlist(sapply(fnode,function(x) ifelse(is.na(x)==F,which(x==RECsample$NZReach),NA)))
  tnode<-unlist(sapply(tnode,function(x) ifelse(is.na(x)==F,which(x==RECsample$NZReach),NA)))
  
  #Create a 3xn matrix of [UPST,DST,NA] for both x an y (to be reshaped later)
  x=t(cbind(RECsample$segXcentroid[fnode],RECsample$segXcentroid,RECsample$segXcentroid*NA,RECsample$segXcentroid,RECsample$segXcentroid[tnode],RECsample$segXcentroid*NA))
  y=t(cbind(RECsample$segYcentroid[fnode],RECsample$segYcentroid,RECsample$segXcentroid*NA,RECsample$segYcentroid,RECsample$segYcentroid[tnode],RECsample$segXcentroid*NA))
  
  #Create a 3xn matrix of [UPST,DST,NA] for both x an y (to be reshaped later)
  #Convert fnode and tnode back into row indexes
  fnodes<-unlist(sapply(fnodes,function(x) ifelse(is.na(x)==F,which(x==REC_plot$NZReach),NA)))
  tnodes<-unlist(sapply(tnodes,function(x) ifelse(is.na(x)==F,which(x==REC_plot$NZReach),NA)))
  x_s=t(cbind(REC_plot$segXcentroid[fnodes],REC_plot$segXcentroid,REC_plot$segXcentroid*NA,REC_plot$segXcentroid,REC_plot$segXcentroid[tnodes],REC_plot$segXcentroid*NA))
  y_s=t(cbind(REC_plot$segYcentroid[fnodes],REC_plot$segYcentroid,REC_plot$segXcentroid*NA,REC_plot$segYcentroid,REC_plot$segYcentroid[tnodes],REC_plot$segXcentroid*NA))
  
  
  #browser()
  linet<-seq(1:8)*myScale
  test=0
  if (length(RECvar>0)){    
    
  for (i in  min(REC_plot$ORDER):max(REC_plot$ORDER)) {   #Loop through the different order 
    inds<-which(REC_plot$ORDER!=i) 
    x2s<-x_s; x2s[ ,inds]<-NA; dim(x2s)=c(length(x2s),1)    #Reshape the matrix into a vector
    y2s<-y_s; y2s[ ,inds]<-NA; dim(y2s)=c(length(y2s),1)
    lines(x2s, y2s, lwd = linet[i], col = "grey80")
    rm(x2s,y2s) 
    
    for (n in 1:n.cat){
      ind<-which(RECsample$ORDER!=i|RECvar!=n) #identify those rows that are NOT to be plotted
      test=test+(length(RECsample$ORDER)-length(ind))
      x2<-x; x2[ ,ind]<-NA; dim(x2)=c(length(x2),1)    #Reshape the matrix into a vector
      y2<-y;  y2[ ,ind]<-NA; dim(y2)=c(length(y2),1)
      lines(x2, y2, lwd = linet[i], col = RECvar.col[n])
      rm(x2,y2)
    }
  }
  }

  #browser()
  legend(x=pos, legend=paste(class.names, " (", PropClass, ")", sep=""), title="Class (% domain)", col = RECvar.col,
          pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)   
  if(!is.null(PlotTitle)==T) title(main=PlotTitle, col.main="red")  

  
}

IDfNode <- function(REC = REC.dat) { # Identify the fnode and tnode nodea
  fnode<-REC$segYcentroid*NA    # Identify the node
  tnode<-REC$segYcentroid*NA    # Identify the node
  for (n in 1:length(fnode)){
    val<-which(REC$Nztnode==REC$Nzfnode[n])
    if (length(val)>=1) {
      if (length(val)>=1) {
        ind=which(REC$ORDER[val]==max(REC$ORDER[val]))
        fnode[n]<-REC$NZReach[val[ind[1]]]
      }else{
        fnode[n]<-REC$NZReach[val[1]]}}       
    val2<-which(REC$Nztnode==REC$Nztnode[n])
    if (length(val2)>=1) {
      if (length(val2)>=1) {
        ind=which(REC$ORDER[val2]==max(REC$ORDER[val2]))
        tnode[n]<-REC$NZReach[val2[ind[1]]]
      }else{
        tnode[n]<-REC$NZReach[val2[1]]}}
  }
    
    return(cbind(fnode,tnode))
  }

###############################################################################
#        Map any continuous variable for REC reaches  in n.breaks colours            #
#                                                                             #
###############################################################################    

map.REC.var <- function(RECvar, REC = REC.dat, n.breaks=10, col.ord = F, name="Variable", pos="topleft", MyRound = 1, ...) {
 
xy <-  cbind(REC$segXcentroid, REC$segYcentroid) 
par(mar = c(0.1,0.1,4, 0.1))

plot(xy, type = "n", xlab=" ", ylab=" ", axes=FALSE, ...)
polygon(x = NZmap$NZE, y = NZmap$NZN, col = "grey95", border = "grey80")

divs <- unique(as.vector(quantile(RECvar, probs = seq(0, 1, by = 1/n.breaks), na.rm = TRUE)))

RECvar.col <-  rainbow(n=n.breaks, end=4/6)
if(col.ord==T) RECvar.col <-  rev(rainbow(n=n.breaks, end=4/6)) 

div.RECvar.n <- cut(RECvar, breaks = divs, labels = F, include.lowest = T) 
xycol <- sapply(div.RECvar.n, function(n, cols = RECvar.col) {return(RECvar.col[n])})

points(xy,  col = xycol, pch = 15, ...)

# points(mySites, cex = 0.75, pch = 19, col = "black")
legend(x=pos, legend = pretty.cut.labs(divs, p=name, rounding = MyRound), col = RECvar.col,
	text.col = "black", pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)   
      
}

# create legend labels for map.REC
pretty.cut.labs <- function(cuts, RECvar, p, s1 = "<", s2 = "<=", rounding = MyRound) {
  if(any(cuts == 0)) {
    if(any(cuts<0)){
      cuts <- paste("-",(formatC(-cuts, digits = min(1, 1 - floor(log10(median(-cuts)))), format = "f") ))
    }else{
  cuts <- formatC(cuts, digits = max(1, 1 - floor(log10(median(cuts)))), format = "f") 
    }
  } else {
  cuts <- round(cuts, rounding)
  }
	n <- length(cuts)
	labs <- c(paste(p, s2, cuts[2]), paste(cuts[2:(n - 2)], s1, p, s2, cuts[3:(n - 1)]), paste(cuts[n - 1], s1, p))		
	return(labs)
}



###############################################################################
#        Map any CATEGORICAL variable for REC reaches                         #
###############################################################################    

map.REC.cat <- function(RECcat, REC = REC.dat, col.ord = T, col.rev = T, name="Classes", pos="topleft", props=F, ...) {
 #col.ord  = FALSE will scramble the colour schememe and may help differentiation
 # props=T will write the proportion of NZReach in each cat in the legend
NZmap <- read.csv("NZmap.csv")  # add coastline
NZmap <- NZmap[NZmap$Section == "North Island" | NZmap$Section == "South Island", ]
 
xy <-  cbind(REC$segXcentroid, REC$segYcentroid) 
par(mar = c(0.1,0.1,4, 0.1))

plot(xy, type = "n", xlab=" ", ylab=" ", axes=FALSE, main=name)
polygon(x = NZmap$NZE, y = NZmap$NZN, col = "grey95", border = "grey80")

cats <- table(RECcat)
n.cats <- length(cats) 
RECcat.col <- rainbow(n.cats)

if (props==T) { 
 prop.class <- round((cats/length(RECcat))*100,1) 
cat.labs <- paste(names(cats), " ; ", prop.class, "%") 
  } else {
   cat.labs <- names(cats)
   } 

      if (col.ord==FALSE) {
      rands <- sample(1:n.cats)
      RECcat.col <- RECcat.col[rands]
      }
      
if(col.rev==T) RECcat.col <- rev(RECcat.col) 
            
xycol <- RECcat.col[match(RECcat,names(cats))] 
points(xy, col = xycol, pch = 15, ...)
  
if(!is.null(pos))  
legend(pos, title=name, legend = cat.labs, col = RECcat.col,
text.col = "black", pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)         
}



