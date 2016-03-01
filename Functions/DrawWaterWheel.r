DrawWaterWheel <- function(Parval,Parlim,PlotTitle,MinQ=NULL,AllocQ=NULL,LandManagement=NULL,sub2=""){
nPar<-length(Parval)
ang=360/nPar #this defines the angle between spokes on the wheel

if (length(PlotTitle)==0) PlotTitle<-"Wheel of Water"

#Define the thickness of the background and foreground bars for the chart
th1<-0.05
th2<- 0.03
innerC<-th2/sin(ang/2*pi/180) #This defines the size of the centre dot

# Prepare the basline circle information
x1<-seq(0,1+innerC,0.01)
y1<-(1-x1^2)^0.5 
x<-c(-rev(x1),x1,rev(x1),-x1)
y<-c(rev(y1),y1,-rev(y1),-y1)

#Define the colour values for the wheel
 Cval<-c(554,90,652,257)

#create figure window
#dev.new()
#x11()

# browser()
#Plot the background circles
plot(x,y,type='l',col='grey75',lty=1,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1,axes=FALSE,xlab="",ylab="",lwd=2)
lines(x*((0.75+innerC)/(innerC+1)),y*((0.75+innerC)/(innerC+1)),type='l',col='green',lty=2,lwd=2)
lines(x*((0.5+innerC)/(innerC+1)),y*((0.5+innerC)/(innerC+1)),type='l',col='goldenrod1',lty=2,lwd=2)
lines(x*((0.25+innerC)/(innerC+1)),y*((0.25+innerC)/(innerC+1)),type='l',col='red',lty=2,lwd=2)

if(is.null(MinQ)==T){
title(main=PlotTitle, col.main="red")
}else{
  #title(main=paste(PlotTitle,"\n","Min flow=",MinQ,"xMALF, Alloc=",AllocQ,"xMALF, Land Management=",LandManagement,"\n",sub2), col.main="red")
  title(main=PlotTitle, col.main="red")
  if (is.null(LandManagement)){
  text(0,1.57,paste("Min flow=",MinQ,"xMALF, Alloc=",AllocQ,"xMALF"),col="blue",cex=1)
  }else{
    text(0,1.57,paste("Min flow=",MinQ,"xMALF, Alloc=",AllocQ,"xMALF, Land Management=",LandManagement),col="blue",cex=1)
  }
  if (!is.null(sub2)) text(0,1.4,sub2,col="blue",cex=1.2)
}
#______________________________________________________________________________
#Plot the background grey bars
for (i in 1:nPar) {
    xn<-c(0,0,cos(((i-1)*ang+90)*pi/180),cos(((i-1)*ang+90)*pi/180))
    yn<-c(0,0,sin(((i-1)*ang+90)*pi/180),sin(((i-1)*ang+90)*pi/180))

    xm=(cos(((i-1)*ang+180)*pi/180)*th1);
    ym=(sin(((i-1)*ang+180)*pi/180)*th1);

    xm1<-c(-xm,xm,xm,-xm)
    ym1<-c(-ym,ym,ym,-ym)

    polygon(xn+xm1,yn+ym1,col='grey85',border=NA)
    }
#______________________________________________________________________________
#Plot the coloured bars for the indicator values
#browser()
for (i in 1:nPar) {

    xn<-c(0,0,cos(((i-1)*ang+90)*pi/180),cos(((i-1)*ang+90)*pi/180))
    yn<-c(0,0,sin(((i-1)*ang+90)*pi/180),sin(((i-1)*ang+90)*pi/180))

    xm<-(cos(((i-1)*ang+180)*pi/180)*th1)
    ym<-(sin(((i-1)*ang+180)*pi/180)*th1)

    xm2<-c(-xm,xm,xm,-xm)
    ym2<-c(-ym,ym,ym,-ym)
    
    #browser()
    PARsc<-approx(Parlim[[i]],c(0.1,0.25,0.5,0.75,1),Parval[[i]])

    Lev<-ceiling(PARsc$y*4)

    if (!is.na(PARsc$y)==T){
    polygon(xn*(PARsc$y+innerC)/(1+innerC)+xm2*th2/th1,yn*(PARsc$y+innerC)/(1+innerC)+ym2*th2/th1,col=colors()[Cval[Lev]],border=NA)  
    text(1.2*cos(((i-1)*ang+90)*pi/180),1.2*sin(((i-1)*ang+90)*pi/180),names(Parval[i]))
    }else{
      polygon(xn+xm2,yn+ym2,col='grey97',border=NA)
    }
    
    }
#______________________________________________________________________________
# Plot a black circle on the middle
polygon(x*innerC,y*innerC,col='grey80',border=NA)
}