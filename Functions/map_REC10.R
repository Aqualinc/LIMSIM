###############################################################################    
MapRivers <- function(RECvar=REC.dat$usArea, 
                      REC = REC.dat, 
                      BlockColor = F, 
                      MinOrder = 0, 
                      LogValues = F, 
                      pos="topleft", 
                      n.breaks=10, 
                      MyLegend = NULL,linevar=TRUE, 
                      MyRound = 1, myCol = "lightblue", 
                      myScale = 0.5, name="Variable", 
                      col.ord=F,
                      MyBreaks = NULL,MySeq = NULL,PlotTitle=NULL,
                      p_leg=TRUE,is.cat=FALSE,props=TRUE,col.rnd=FALSE,col.lookup=FALSE,
                      do.sea.shade = T,                                 # apply coastal sea shading?
                      do.hill.shade = F,
                      do.lake = T,
                      min.lake.area = 10^6,                      # smallest lake areas to be plotted
                      do.RC = T,                              # include regional council boundaries?
                      ...) {
  
  #BlockColor = T will produce a grey river lines map.
  
  #browser()
if(do.sea.shade | do.hill.shade | do.RC | do.lake) { 
                        if(!exists("nzhires")) stop("Load NZmap_hires.RData from C:/TonsData/MyTools/Tons R SCRIPTS/MartinsMapping/NZmap_hires.RData")
}
  
  if(LogValues) {  # better breaks if the data are highly skewed  
    LogRECvar <- log10(RECvar)
    bb <- hist(log10(RECvar), breaks=n.breaks)
    LogBreaks <- bb$breaks
    CutBreaks <- cut(LogRECvar, breaks= LogBreaks) 
    UnLogBreaks <- as.vector(by(LogRECvar, CutBreaks, min))
    myBreaks <- RECvar[match(UnLogBreaks, LogRECvar)]
    if(any(is.infinite(LogRECvar)))  myBreaks <- c(0, myBreaks)
    MyBreaks <- myBreaks 
  }
  
  TheseSites <- which(REC$ORDER > MinOrder & !is.na(RECvar)==T) 
  
  RECsample <- REC[TheseSites, ]
  REC_plot<-REC[which(REC$ORDER>MinOrder),]
  xy <-  cbind(REC$segXcentroid, REC$segYcentroid) 
  par(mar = c(0.1,0.1,2, 0.1))  
  plot(xy, type = "n", xlab=" ", ylab=" ", axes=FALSE, asp=1, ...)
  
  if (do.sea.shade) {points(x = nzhires$nzmge, y = nzhires$nzmgn,
                            col = "aliceblue", pch = 16, cex = 10)}
  
  polygon(x = NZmap$NZE, y = NZmap$NZN, col = "grey95", border = "grey80")
  
  if (do.hill.shade) {
    ynorth <- as.numeric(names(nz.hill.shade.xy))             # for hill shading
    cellsize <- 100
    shade.y.resolution <- 5
    myRamp <- colorRamp(c("grey65", "grey100"))       # specify colours to taste
    cols.alt <- rgb(myRamp(seq(0, 1, length = 100)), max = 255)  # colour vector
    theserows <- seq(1, length(ynorth), by = shade.y.resolution)
    junk <- lapply(ynorth[theserows], function(k) {
      xeast <- nz.hill.shade.xy[[as.character(k)]]
      if (length(xeast) > 0) {
        altcols <- sapply(round(100 * xeast$shade / 256, 0),
                          function(x) {cols.alt[x + 1]})
        points(xeast$east, rep(k, nrow(xeast)), pch = ".", col = altcols)
      }
    }
    )
  }
  
  
  
  RECvar <-  RECvar[TheseSites]
  
  if (is.cat==F){
    if(is.null(MyLegend)) {
      if(is.null(MyBreaks)&is.null(MySeq)){
        divs <- unique(as.vector(quantile(RECvar, probs = seq(0, 1, by = 1/n.breaks), na.rm = TRUE)))
      }else{
        if(is.null(MyBreaks)){
          divs<-as.vector(seq(MySeq[1],MySeq[2],MySeq[3]))  
        }else{
          if(max(MyBreaks)<max(RECvar)){
            divs<-c(MyBreaks,max(RECvar))
          }else{
            divs<-MyBreaks
          }
          
        }
      }
      n.breaks<-length(divs)-1
      RECvar.col <-  rainbow(n=n.breaks, end=4/6)
      if(col.ord==T) RECvar.col <-  rev(rainbow(n=n.breaks, end=4/6)) 
      # browser()
      if (length(RECvar>0)){
        div.RECvar.n <- cut(RECvar, breaks = divs, labels = F, include.lowest = T) 
        xycol <- sapply(div.RECvar.n, function(n, cols = RECvar.col) {return(RECvar.col[n])})
      }
    } else {
      newBreaks <- MyLegend$divs
      # ensure the upper and lower limits comply
      if (newBreaks[1]>min(RECvar,na.rm=T)) newBreaks[1]<-min(RECvar,na.rm=T)
      if (newBreaks[length(newBreaks)]<max(RECvar,na.rm=T)) newBreaks[length(newBreaks)]<-max(RECvar,na.rm=T)
      
      RECvar.col <-  rainbow(n=length(MyLegend$divs)-1, end=4/6)
      if(MyLegend$col.ord==T) RECvar.col <-  rev(RECvar.col) 
      if (length(RECvar>0)){
        div.RECvar.n <- cut(RECvar, breaks = newBreaks, labels = F, include.lowest = T) 
        xycol <- as.character(sapply(div.RECvar.n, function(n) {return(RECvar.col[n])}))
      }
      divs <- newBreaks
    }
  }else{
    cats <- levels(RECvar)
    n.cats <- length(cats) 
    
    if (col.lookup[1]==FALSE){
      RECvar.col <- rainbow(n.cats)
      if (col.rnd==TRUE) {
        rands <- sample(1:n.cats)
        RECvar.col <- RECvar.col[rands]
      }
      if (col.ord==FALSE){RECvar.col <- rev(RECvar.col)}
    }else{
      #browser()
      RECvar.col<-col.lookup[match(cats,col.lookup[,1]),2]
    }
    
    if (props==T) { 
      prop.class <-round(sapply(cats,function(x) length(which(RECvar==x))/length(RECvar)*100),1)
      #prop.class <-round((cats/length(RECvar))*100,1) 
      cat.labs <- paste(as.character(cats), " ; ", prop.class, "%") 
    } else {
      cat.labs <- as.character(cats)
    } 
    xycol <- RECvar.col[match(RECvar,cats)] 
    divs<-cats
    
  }

#browser()
  
  if(BlockColor) { # this is a hack to swithch off the colour off (switch off legend too) 
    RECvar.col <- rep("blue", length(RECvar.col))
    p_leg <- F
  }
  
  if(is.null(REC_plot$fnode)) {
    print("For faster plotting next time run function 'IDfNode' and save ouput to REC.dat as '$fnode' & 'Stnode' ") 
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
  
  #print("Convert fnode and tnode back into row indexes")
  fnodes<-match(fnodes,REC_plot$NZReach)
  tnodes<-match(tnodes,REC_plot$NZReach)
  x_s=t(cbind(REC_plot$segXcentroid[fnodes],REC_plot$segXcentroid,REC_plot$segXcentroid*NA,REC_plot$segXcentroid,REC_plot$segXcentroid[tnodes],REC_plot$segXcentroid*NA))
  y_s=t(cbind(REC_plot$segYcentroid[fnodes],REC_plot$segYcentroid,REC_plot$segXcentroid*NA,REC_plot$segYcentroid,REC_plot$segYcentroid[tnodes],REC_plot$segXcentroid*NA))
  
  fnode<-match(fnode,RECsample$NZReach)
  tnode<-match(tnode,RECsample$NZReach)
  x=t(cbind(RECsample$segXcentroid[fnode],RECsample$segXcentroid,RECsample$segXcentroid*NA,RECsample$segXcentroid,RECsample$segXcentroid[tnode],RECsample$segXcentroid*NA))
  y=t(cbind(RECsample$segYcentroid[fnode],RECsample$segYcentroid,RECsample$segXcentroid*NA,RECsample$segYcentroid,RECsample$segYcentroid[tnode],RECsample$segXcentroid*NA))
  

#browser()
  #MAKE MAP
  if(linevar==T){linet<-seq(1:8)*myScale}else{linet<-matrix(myScale,nrow=8,ncol=1)}
  test=0
  if (length(RECvar>0)){   
    
    for (i in  min(REC_plot$ORDER):max(REC_plot$ORDER)) {   #Loop through the different order 
      inds<-which(REC_plot$ORDER!=i) 
      x2s<-x_s; x2s[ ,inds]<-NA; dim(x2s)=c(length(x2s),1)    #Reshape the matrix into a vector
      y2s<-y_s; y2s[ ,inds]<-NA; dim(y2s)=c(length(y2s),1)
      lines(x2s, y2s, lwd = linet[i], col = "grey80")
      rm(x2s,y2s) 
     
      for (n in 1:length(RECvar.col)){
        ind<-which(RECsample$ORDER!=i|xycol!=RECvar.col[n]) #identify those rows that are NOT to be plotted
        test=test+(length(RECsample$ORDER)-length(ind))
        x2<-x; x2[ ,ind]<-NA; dim(x2)=c(length(x2),1)    #Reshape the matrix into a vector
        y2<-y;  y2[ ,ind]<-NA; dim(y2)=c(length(y2),1)
        lines(x2, y2, lwd = linet[i], col = RECvar.col[n])
        rm(x2,y2)
      }
    }
  }
  
  # add lake polygons
if(do.lake) {
  myLakes <- allLakes[allLakes$Area >= min.lake.area, ]
  if (nrow(myLakes) > 0) {
    Nskip <- 10       # number of boundary points to skip between plotted points
    polysToPlot <- merge(lakes, myLakes[ , c("LakeID", "Area")], by = "LakeID")
    polysToPlot <- polysToPlot[order(polysToPlot$lpID), ]
    polysToPlot <- polysToPlot[polysToPlot$lpID%%Nskip == 0 |
                                 is.na(polysToPlot$east), ]
    polygon(x = polysToPlot$east, y = polysToPlot$north,
            col = "lightblue3", border = "lightblue3")}
}
  # add RC boundaries
  if (do.RC) points(North ~ East, data = RCbounds, type = "l", col = "grey40")
  # add scale and north arrow
#   if (do.scale) {
#     lines(x = c(2600000, 2800000), y = rep(5500000, 2))
#     lines(x = rep(2600000, 2), y = c(5495000, 5505000))
#     lines(x = rep(2800000, 2), y = c(5495000, 5505000))
#     text(x = 2700000, y = 5530000, "200 km", font = 2, adj = 0.5)
#     lines(x = rep(2900000, 2), y = c(5400000, 5600000))
#     polygon(x = c(2900000, 2890000, 2910000, 2900000), y = c(5600000, 5580000,
#                                                              5580000, 5600000), col = "black", border = "black")
#   }
  
  #browser()
  if (p_leg==T){
    if(is.cat==F){
      legend(x=pos, legend = pretty.cut.labs(divs, p=name, rounding = MyRound), col = RECvar.col,
             pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)   
    }else{
      legend(pos, title=name, legend = cat.labs, col = RECvar.col,
             text.col = "black", pch = c(15), bg = 'white', inset = c(0.05, 0.15), cex = 0.8)    
    }
  }else{
    
  }
  if(!is.null(PlotTitle)==T) title(main=PlotTitle, col.main="red")
  
  TheLegend <- list(divs=divs, col.ord=col.ord) 
  return(TheLegend)
  
}  # end MapRivers


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