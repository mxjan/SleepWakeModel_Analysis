######################## 
# SOME FUNCTIONS TO MANIPULATE GGPLOT
########################


RemoveXaxis<-function(gg){
  gg<-gg + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  gg <- gg + xlab(NULL) #added
  return(gg)
  
}

RemoveYaxis<-function(gg){
  
  gg<-gg + theme(axis.title.y=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y.left=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.line.y = element_blank())
  gg <- gg + ylab(NULL)
  return(gg)
  
}

RemoveXlab<-function(gg){
  gg <- gg + xlab(NULL)
  return(gg)
}

RemoveYlab<-function(gg){
  gg <- gg + ylab(NULL)
  return(gg)
}

GetTitleH<-function(Text,angle=0,bg="white",fontface="plain",fontCol="black",size=10,boxcolor=NULL,margin=c(0,0,0,0)){
  gg<-ggdraw() +draw_label(
    Text,
    color = fontCol,
    size = size,
    fontface=fontface,
    angle=angle
  ) +
    theme(plot.margin = unit(margin,"cm"),
          plot.background = element_rect(fill = bg,color=bg))
  if (! is.null(boxcolor)){
    gg<-gg+theme(plot.background = element_rect(fill = bg,color=boxcolor))
  }
  return(gg)
}



squish_trans <- function(from, to, factor) {
  require(scales)
  trans <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    if (any(is.na(x))) return(x)
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}

AddDayLine<-function(p){
  p<-p+geom_vline(xintercept = seq(0,252,by=24),linetype="dashed",col="grey50",size=.25)
  return(p)
}

# Rudolf Cardinal, March 2011
# Simple extensions to ggplot2 (v0.8.7); see http://pobox.com/~rudolf/statistics/R
# Modified 5 Jan 2013 for ggplot2 0.9.3 (NB: use sessionInfo() to find current package versions)
# - fetch ggplot2 source with: git clone https://github.com/hadley/ggplot2.git
# Changes, because ggplot2 has changed its internal calling mechanisms:
# - opts() deprecated in favour of theme()
# - "Element panel.border must be an element_rect object" (error from validate_element() in theme-elements.r)
#   ... so change all class = "theme" to class = c("element_rect", "element")
# - "cannot coerce type 'closure' to vector of type 'list'"
#   ... a closure is a function (see ?typeof)
#   ... change class to be of class c("MYCLASS", "element_rect", "element")
# - then element_grob.MYCLASS not called by element_render()/element_grob()/UseMethod()... environment/namespace problem
#   tried setMethod("element_grob", "theme_border", function(STUFF) { STUFF} , where = as.environment("package:ggplot2")
#   but the environment is locked
#   ggplot2's theme-elements.r defines e.g. element_rect (exported) and element_grob.element_rect (not exported, does the work)
#   However, we can't override an internal function:
#       ... e.g. rewrite "validate_element" to crash
#           set environment(validate_element) <- as.environment("package:ggplot2") -- doesn't break the plotting.
# - Upshot: now impossible to hack through this way (locked environment).
# - http://obeautifulcode.com/R/How-R-Searches-And-Finds-Stuff/
# - http://stackoverflow.com/questions/8204008/redirect-intercept-function-calls-within-a-package-function
# - These don't fix it:
#   library(proto)
#   theme <- with(proto(environment(ggplot2::theme), theme = ggplot2::theme, element_grob.theme_border = my.element_grob.theme_border), theme) --- doesn't work
#   ggplot <- with(proto(environment(ggplot2::ggplot), ggplot = ggplot2::ggplot, element_grob.theme_border = my.element_grob.theme_border), ggplot) --- breaks!
# - Fix by Baptiste Auguie 8/1/2013: inherit from element_blank instead; then it works fine.

#-------------------------------------------------------------------------------
# Requirements
#-------------------------------------------------------------------------------

library(grid) # for gpar

#-------------------------------------------------------------------------------
# Code duplicated from ggplot2 source (not exposed to wider namespace) for convenience
#-------------------------------------------------------------------------------

.pt <- 1 / 0.352777778
len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}

#-------------------------------------------------------------------------------
# Generic panel border (can set any combination of left/right/top/bottom)
#-------------------------------------------------------------------------------

theme_border <- function(
  type = c("left", "right", "bottom", "top", "none"),
  colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + opts( panel.border=theme_border(type=c("bottom","left")) ) + ...
  type <- match.arg(type, several.ok=TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}
element_grob.theme_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (length(type)==0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: "L" (left + bottom) border
#-------------------------------------------------------------------------------

theme_L_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_L_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_L_border", "element_blank", "element")
  )
}
element_grob.theme_L_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x+width, x, x), y = c(y,y,y+height), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: bottom border only
#-------------------------------------------------------------------------------

theme_bottom_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_bottom_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_bottom_border", "element_blank", "element")
  )
}
element_grob.theme_bottom_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x, x+width), y = c(y,y), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}

#-------------------------------------------------------------------------------
# For convenience: left border only
#-------------------------------------------------------------------------------

theme_left_border <- function(colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_left_border() ) + ...
  structure(
    list(colour = colour, size = size, linetype = linetype),
    class = c("theme_left_border", "element_blank", "element")
  )
}
element_grob.theme_left_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = c(x, x), y = c(y, y+height), ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}



#-------------------------------------------------------------------------------
# Border selection by number
#-------------------------------------------------------------------------------

theme_border_numerictype <- function(type, colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + theme( panel.border=theme_border(type=9) ) + ...
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border_numerictype", "element_blank", "element")
  )
}
element_grob.theme_border_numerictype <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  # numerical types from: library(gridExtra); example(borderGrob)
  # 1=none, 2=bottom, 3=right, 4=top, 5=left, 6=B+R, 7=T+R, 8=T+L, 9=B+L, 10=T+B, 11=L+R, 12=T+B+R, 13=T+L+R, 14=T+B+L, 15=B+L+R, 16=T+B+L+R
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if (type==2 || type==6 || type==9 || type==10 || type==12 || type==14 || type==15 || type==16) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if (type==4 || type==7 || type==8 || type==10 || type==12 || type==13 || type==14 || type==16) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if (type==5 || type==8 || type==9 || type==11 || type==13 || type==14 || type==15 || type==16) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if (type==3 || type==6 || type==7 || type==11 || type==12 || type==13 || type==15 || type==16) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (type==1) { # blank; can't pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}




FittingPlot_xscaled<-function(object,params,xlim=NULL,noFit=F,noy=F,
                              othertitle=NA,bsldashed=F,colexpr="forestgreen",bsldashedcol=rgb(70/255,130/255,180/255,1),
                              fitcol="red",ylab="log2 CPM",ExpTransform=F,lineSize=2,lineSizeBSL=1,PointSize=1,MeanPointSize=2,BootstrapResults=NULL){
  
  
  library(ggplot2)
  
  FittedData<-SWDMrFit(object = object,params = params)
  
  if (ExpTransform == T){
    FittedData$y1<-2^(FittedData$y1)-1
  }
  
  Gene<-object@VarExp
  
  if (ExpTransform == T){
    object@Gexp[,Gene]<-2^object@Gexp[,Gene]-1
  }
  dfd<-as.data.frame(cbind(object@Gexp$Time,object@Gexp[,object@VarExp]))
  
  colnames(dfd)<-c("Time","Gene")
  
  p_Gene <- ggplot(dfd,aes(x=Time,y=Gene))
  # LD rect
  for (i in seq(12,96,by=24)){
    p_Gene <- p_Gene + annotate("rect", xmin = i, xmax = i+12, ymin = -Inf, ymax = Inf, fill="black",alpha = 0.2)
  }
  i<-204
  p_Gene <- p_Gene + annotate("rect", xmin = i, xmax = i+12, ymin = -Inf, ymax = Inf, fill="black",alpha = 0.2)
  
  
  # Boxplot + Data points
  p_Gene <- p_Gene + geom_point(color = "black",size=PointSize)
  #p_Gene <- p_Gene + geom_boxplot(fill="forestgreen",aes(group = cut_width(Time, 2))) #+ geom_jitter(shape=16, position=position_jitter(0.2))
  # Mean point and sd per time points
  meansd<-SWDMr:::StatsPerTimePoint(object)
  
  p_Gene <- p_Gene + ylab(ylab)+xlab("Time")
  #p_Gene <- p_Gene + scale_x_continuous(breaks=seq(0,max(xlim),12),limits = xlim)
  p_Gene <- p_Gene + scale_x_continuous(trans = squish_trans(108, 203, 30),
                                        breaks = c(seq(0, 96, by = 12),204,216,228),limits=xlim)
  p_Gene <- p_Gene + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA)) 
  
  
  p_Gene <- p_Gene + annotate(x=unique(object@Gexp$Time),"errorbar",ymin=meansd[,"mean"]-meansd[,"se"],ymax=meansd[,"mean"]+meansd[,"se"],colour="grey25", width=1,size=1)
  rangev<-ggplot_build(p_Gene)$layout$panel_params[[1]]$y.range
  p_Gene <- p_Gene + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = rangev[[1]]+(rangev[[2]]-rangev[[1]])*0.1, fill="red",alpha = 0.3)
  p_Gene <- p_Gene + annotate("point",x=unique(object@Gexp$Time),y=meansd[,"mean"],size=MeanPointSize,col="black",shape=21,fill=colexpr)
  
  if (bsldashed ==T){
    bslF<-rep(FittedData$y1[FittedData$time > 0.0 & FittedData$time <=24.0],8)
    p_Gene <- p_Gene + annotate("line",x=seq(48.1,240,by=.1),y=bslF,color=bsldashedcol,size=lineSizeBSL,alpha=1,linetype = "dashed") 
  }
  
  # SD rect
  if (noFit == F){
    p_Gene <- p_Gene + annotate("line",x=FittedData$time[FittedData$time>=0],y=FittedData$y1[FittedData$time>=0],color=fitcol,size=lineSize,alpha=.8)
  }
  
  if (noy ==T ){
    p_Gene <- p_Gene + theme(axis.title.y=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank())
  }
  if (! is.na(othertitle)){
    p_Gene <- p_Gene + ggtitle(othertitle)
  }else{
    p_Gene <- p_Gene + ggtitle(Gene)
  }
  
  if (! is.null(BootstrapResults)){
    p_Gene <-p_Gene + annotate("ribbon",x=BootstrapResults$df.mc$Time, ymin=BootstrapResults$df.mc$lwr.pred, ymax=BootstrapResults$df.mc$upr.pred, alpha=0.2, fill="blue")
    p_Gene <-p_Gene + annotate("ribbon",x=BootstrapResults$df.mc$Time, ymin=BootstrapResults$df.mc$lwr.conf, ymax=BootstrapResults$df.mc$upr.conf, alpha=0.4, fill="#339900")
  }
  
  
  p_Gene<-p_Gene + annotate("rect",xmin=105,xmax=204, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
  
  return(p_Gene)
}


PlotMeGene<-function(Gene,Tissue){
  library(SWDMr)
  source("../../Analysis/RFunctions/C57BL6_ModelBuilding.R",chdir=T)
  source("../../Analysis/RFunctions/PlotFitMouse.R")
  source("../../Analysis/RFunctions/color.R")
  source("../../Analysis/RFunctions/SWdrivenPlot.R")
  
  if (Tissue == "Liver"){
    load("../../Data/C57BL6J_Liver_TimeCourse_MJCHN2020/B6Liver_FilePrepSWDMr.Rdata")
    fits<-read.table("../../Results&Figures/CHN2019_Fitting/C57BL6_LiverGenesFits_Full_131020.txt",header=T)
    fits$dampratio<-fits$dampratio/(2*fits$omega)^2
  }else{
    load("../../Data/C57BL6J_Cortex_TimeCourse_CHN2019/B6Cortex_FilePrepSWDMr.Rdata")
    fits<-read.table("../../Results&Figures/CHN2019_Fitting/C57BL6_CortexGenesFits_Full_131020.txt",header=T)
    fits$dampratio<-fits$dampratio/(2*fits$omega)^2
  }

  cols<-ColorCode()
  swdmr <- SWDMr(SWdist=SWdf, Gexp=rna_expr)
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  
  pFUll1<-FittingPlot_xscaled(model,fits[fits$Gene == Gene,],xlim=c(0,228),
                              noFit = F,bsldashed = T,colexpr = cols[[Tissue]],
                              bsldashedcol = cols[[Tissue]],fitcol = cols[["Modelfit"]],
                              ExpTransform=F,PointSize = 2,lineSize = 1,lineSizeBSL = 1,
                              MeanPointSize=3)
  pFUll1<-pFUll1+theme(plot.margin = margin(1, 1, 1, 1),axis.text.x = element_text(angle = 45, hjust = 1),title =element_text(size=10, face='bold'))
  return(pFUll1)
  
}


PlotMeProbe<-function(ProbeID){
  library(SWDMr)
  source("../../Analysis/RFunctions/Human_ModelBuilding.R")
  source("../../Analysis/RFunctions/Plot_HT_SleepData_OneProbe.R")
  source("../../Analysis/RFunctions/color.R")
  source("../../Analysis/RFunctions/SWdrivenPlot.R")
  load("../../Data/HumanDataset/HS_BloodTranscriptomeAndSleepWake.Rdata")
  cols<-ColorCode()
  
  fits<-read.table("../../Results&Figures/HumanFits/Human_BloodProbesFits_Full_111120.txt",header=T)
  fits$dampratio<-fits$dampratio/(2*fits$omega)^2
  pFull<-HT_PlotProbe(ProbeID,fits,cols)
  pFull<-pFull+theme(plot.margin = margin(1, 1, 1, 1),axis.text.x = element_text(angle = 45, hjust = 1),title =element_text(size=10, face='bold'))
  return(pFull)
  
}
