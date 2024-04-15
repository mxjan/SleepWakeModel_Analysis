GetCol<-function(){
  library(RColorBrewer)
  col_vector<-c(brewer.pal(9, "Set1"),brewer.pal(8, "Dark2"),brewer.pal(8, "Accent"))
  return(col_vector)
}

ShapeCond<-function(){
  return(c("NSD"=19,"SD" = 17))
}


AddRectDarkggplot2<-function(ggplotobj,maxv=96,ymax=Inf,color=rgb(0,0,0,0),AddLightRect=F){
  for (i in seq(12,maxv,by=24)){
    ggplotobj <- ggplotobj + annotate("rect", xmin = i, xmax = i+12, ymin = -Inf, ymax = ymax, fill="black",alpha = 0.2,color=color)
    if (AddLightRect == T){
      ggplotobj <- ggplotobj + annotate("rect", xmin = i-12, xmax = i, ymin = -Inf, ymax = ymax, fill="white",alpha = 0.2,color=color)
    }
  }
  
  return(ggplotobj)
}
  

AddRectSDggplot2<-function(ggplotobj,ymax="scaled"){
  rangev<-ggplot_build(ggplotobj)$layout$panel_params[[1]]$y.range
  if (ymax == "scaled"){
    ggplotobj <- ggplotobj + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = rangev[[1]]+(rangev[[2]]-rangev[[1]])*0.1, fill="red",alpha = 0.3)
  }else{
    ggplotobj <- ggplotobj + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = ymax, fill="red",alpha = 0.3)
  }
  
  return(ggplotobj)
}


ColorCode<-function(){
  
  colcode<-c()
  
  col<-rgb(154/255, 164/255, 157/255)#rgb(194/255, 204/255, 197/255)
  names(col)<-"Wake"
  colcode<-c(colcode,col)
  
  col<-rgb(20/255,46/255, 84/255)
  names(col)<-"Sleep"
  colcode<-c(colcode,col)
  
  #col<-"firebrick4"
  col<-"grey10"
  names(col)<-"Modelfit"
  colcode<-c(colcode,col)
  
  #col<-rgb(20 /255,61/255, 89/255)
  col<-"#5D5BA4"
  names(col)<-"SWforce"
  colcode<-c(colcode,col)
  
  col<-rgb(224/255,160/255, 0/255)#rgb(244/255,180/255, 26/255)
  names(col)<-"SCNforce"
  colcode<-c(colcode,col)
  
  col<-rgb(40/255,169/255, 225/255)
  names(col)<-"Cortex"
  colcode<-c(colcode,col)
  
  # col<-rgb(190/255,22/255, 33/255)
  # names(col)<-"Liver"
  # colcode<-c(colcode,col)
  
  col<-rgb(103/255,76/255, 71/255)
  names(col)<-"Liver"
  colcode<-c(colcode,col)
  
  col<-rgb(138/255,3/255, 3/255)
  names(col)<-"Blood"
  colcode<-c(colcode,col)
  
  return(colcode)
  
}

GetContZTCols<-function(ZT,colorres=100){
  require(pals)
  library(colorspace)
  #colorband<-kovesi.cyclic_mrybm_35_75_c68_s25(colorres)
  #colorband<- ocean.phase(colorres)
  
  # col selection
  # colsel<-colorRampPalette(c("red","green","blue","purple","orange","red"))
  # colsel<-colorRampPalette(brewer.set1(5)[c(1,3,2,4,5,1)],interpolate ="linear")
  
  colsel<-colorRampPalette(c("darkred",darken("gold",.2),"darkorange1","#7F7100","darkgreen","#00327F","royalblue4","#45007F","darkred"),interpolate ="linear")
  #plot(seq(1,100)/100,col=colsel(100),pch=19,cex=2)
  
  colorband<-colsel(100)
  
  #z=matrix(1:100,nrow=1)
  #x=1
  #y=seq(3,2345,len=100) # supposing 3 and 2345 are the range of your data
  #image(x,y,z,col=colsel(100),axes=FALSE,xlab="",ylab="")
  
  return(colorband[as.integer((ZT %% 24)/24*colorres)+1])
}


ZTcolor<-function(Tissue){
  library(grDevices)
  library(colorspace)
  color<-GetCol()
  #BaseColor<-color[[Tissue]]
  
  # Generate cyclic color function
  #colfunc <- colorRampPalette(c(rgb(104/255,13/255,31/255), BaseColor,rgb(244/255,167/255,129/255), rgb(214/255,99/255,77/255),rgb(104/255,13/255,31/255)))
  
  cols <- c("0" =color[1] ,"3" = darken(color[2],.2), "6" = color[3], "12" = color[4], "18" = color[5])
  return(cols)
}
