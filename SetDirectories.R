SetDirectories<-function(rootdir){
  direct<-list()
  direct$fun<-paste(rootdir,"/Functions",sep="")
  direct$mod<-paste(rootdir,"/CallModels",sep="") 
  direct$data<-paste(rootdir,"/Datasets",sep="")
  direct$root<-rootdir
  
  return(direct)
}