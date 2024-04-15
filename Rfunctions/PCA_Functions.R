PCAplot<-function(Expr,TimeExpr,fits,fitted_full,fitted_circcomp,fitted_swcomp,
                  linedashedHex="11",xlim=c(-120,100),ylim=c(-60,60),AddULeq=F,title){
  
  AddTheme<-function(p){
    p<-p+theme(panel.background=element_rect(fill = 'white', colour = 'black'),panel.grid.major=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank())
    return(p)
  }
  
  LineSize<-.5
  TextSize<-3
  AlphaEllipse<-.2
  
  library(ggrepel)
  library(FactoMineR)
  library(factoextra)
  
  idx_BF<-fits$Gene[fits$BF>exp(1)]
  idx_BF<-idx_BF[!duplicated(idx_BF)]
  
  TimeFactor<-as.factor(paste("T",TimeExpr,sep=""))
  res.pca <- PCA(cbind.data.frame(Expr[,idx_BF],TIME=TimeFactor),
                 scale.unit=T,quali.sup = dim(Expr[,idx_BF])[2]+1,graph=FALSE,ncp = 10)
  
  message("PC1 top contributors:",paste0(rev(names(tail(sort(res.pca$var$contrib[,1])))),collapse = " "))
  message("PC2 top contributors:",paste0(rev(names(tail(sort(res.pca$var$contrib[,2])))),collapse = " "))
  
  ProjFitC<-predict(res.pca,newdata = t(fitted_full[idx_BF,])) # Projection of fits
  ProjFitcirc<-predict(res.pca,newdata = t(fitted_circcomp[idx_BF,]))
  ProjFitsw<-predict(res.pca,newdata = t(fitted_swcomp[idx_BF,]))
  ellipse.coord = coord.ellipse(cbind.data.frame(TimeFactor, res.pca$ind$coord),bary=T,level.conf = .95) # ellipses per time-points
  
  # Plot data under baseline condition (24h sleep-wake cycle)
  data<-cbind.data.frame(PC1=ProjFitC$coord[1:241,1],PC2=ProjFitC$coord[1:241,2])
  p_pca<-ggplot(aes(PC1,PC2),data=data)
  # Add ellipse
  for (i in unique(TimeFactor)){
    Time<-as.numeric(gsub("T","",i))
    p_pca<-p_pca + annotate("polygon",x=ellipse.coord$res$Dim.1[ellipse.coord$res$TimeFactor == i],
                                  y=ellipse.coord$res$Dim.2[ellipse.coord$res$TimeFactor == i],fill=alpha(GetContZTCols(Time),AlphaEllipse),col=rgb(0,0,0,AlphaEllipse))
  }
  
  # Add data sync (28h sleep-wake cycle)
  p_pca<-p_pca+annotate("path",x=ProjFitC$coord[481:840,1],y=ProjFitC$coord[481:840,2],size=LineSize,col=colorCode[["Modelfit"]],linetype="solid")
  # Plot fitted lines
  p_pca<-p_pca+geom_path(col=colorCode[["Modelfit"]],size=LineSize,linetype=linedashedHex)
  
  # # Add text
  # TimeToWrite<-c("T24","T27","T30","T36","T42","T51","T54","T60","T66")
  # idx_TTW<-ellipse.coord$res$TimeFactor %in% TimeToWrite
  # Agg<-aggregate(ellipse.coord$res[idx_TTW,c(2,3)],list(ellipse.coord$res$TimeFactor[idx_TTW]),mean)
  # Agg$Label<-paste("ZT",as.numeric(gsub("T","",Agg$Group.1)) %% 24,sep="")
  # p_pca<-p_pca+ annotate("text",x=Agg$Dim.1,y=Agg$Dim.2,label=Agg$Label,size=TextSize,color=GetContZTCols(as.numeric(gsub("T","",Agg$Group.1))),fontface=2)
  

  # Add text
  TimeToWrite<-list("T24"=expression(T24[ZT0]),
                    "T27"=expression(T27[ZT3]),
                    "T30"=expression(T30[ZT6]),
                    "T36"=expression(T36[ZT12]),
                    "T42"=expression(T42[ZT18]),
                    "T51"=expression(T51[ZT3]),
                    "T54"=expression(T54[ZT6]),
                    "T60"=expression(T60[ZT12]),
                    "T66"=expression(T66[ZT18]))
  idx_TTW<-ellipse.coord$res$TimeFactor %in% names(TimeToWrite)
  Agg<-aggregate(ellipse.coord$res[idx_TTW,c(2,3)],list(ellipse.coord$res$TimeFactor[idx_TTW]),mean)
  for (i in names(TimeToWrite)){
    p_pca<-p_pca+ annotate("text_repel",x=Agg[Agg$Group.1 == i,"Dim.1"],y=Agg[Agg$Group.1 == i,"Dim.2"],
                           label=TimeToWrite[[i]],
                           size=TextSize,
                           color=GetContZTCols(as.numeric(gsub("T","",i))),
                           fontface=2,
                           box.padding = 0)
  }

  # Add arrows
  p_pca<-p_pca+annotate("segment",x=ProjFitC$coord[560,1],y=ProjFitC$coord[560,2],
                 xend=ProjFitC$coord[565,1],yend=ProjFitC$coord[565,2],
                 col=colorCode[["Modelfit"]], #
                 arrow=arrow(type = "closed",length = unit(.25,"cm")))
  
  p_pca<-p_pca+annotate("segment",x=ProjFitC$coord[390,1],y=ProjFitC$coord[390,2],
                        xend=ProjFitC$coord[395,1],yend=ProjFitC$coord[395,2],
                        col=colorCode[["Modelfit"]], #
                        arrow=arrow(type = "closed",length = unit(.25,"cm")))
  
  # Limits
  p_pca<-p_pca+ylab(paste("PC2 [",round(res.pca$eig[2,2]),"%]",sep=""))+
    coord_cartesian(xlim=c(xlim[1],xlim[2]),ylim=c(ylim[1],ylim[2]))#+xlim(xlim[1],xlim[2])+ylim(ylim[1],ylim[2])
  
  
  p_pca<-AddTheme(p_pca)
  
  ########### ADD SIZE PLOTS
  
  library(patchwork)
  
  dd<-cbind.data.frame(Time=seq(.1,24,by=.1),PC2=ProjFitcirc$coord[1:240,2])
  
  p_pca_right<-ggplot(aes(x=Time,y=PC2),data=dd)+
    geom_path(size=LineSize,col=colorCode[["SCNforce"]],linetype="solid")+
    scale_y_continuous(limits=c(ylim[1],ylim[2]),position = "right")
  p_pca_right<-p_pca_right+
    annotate("path",x=seq(.1,24,by=.1),y=ProjFitsw$coord[1:240,2],col=colorCode[["SWforce"]],linetype=linedashedHex,size=LineSize)
  p_pca_right<-p_pca_right+
    annotate("path",x=seq(.1,24,by=.1),y=ProjFitsw$coord[481:720,2],col=colorCode[["SWforce"]],linetype="solid",size=LineSize)+
    scale_x_continuous(breaks=seq(0,24,by=6),labels=paste(seq(24,48,by=6),seq(48,72,by=6),sep="|"),position = "top")+
    theme(axis.text.x.top = element_text(angle=90,vjust = 0.5))
  
  p_pca+ylim(ylim[1],ylim[2])|p_pca_right
  
  dd<-cbind.data.frame(Time=seq(.1,24,by=.1),PC1=ProjFitsw$coord[1:240,1])
  p_pca_bott<-ggplot(aes(x=PC1,y=Time),data=dd)+geom_path(col=colorCode[["SWforce"]],linetype=linedashedHex,size=LineSize)+
    annotate("path",x=ProjFitsw$coord[481:720,1],y=seq(.1,24,by=.1),col=colorCode[["SWforce"]],linetype="solid",size=LineSize)+
    xlim(xlim[1],xlim[2])+annotate("path",x=ProjFitcirc$coord[481:720,1],y=seq(.1,24,by=.1),col=colorCode[["SCNforce"]],size=LineSize)+
    scale_y_continuous(breaks=seq(0,24,by=6),labels=seq(0,24,by=6))
  (p_pca_bott+xlim(xlim[1],xlim[2]))/(p_pca_bott)
  
  p_pca_bott<-p_pca_bott+ylab("Time [h]")+xlab(paste("PC1 [",round(res.pca$eig[1,2]),"%]",sep=""))+
    scale_y_reverse(breaks=seq(0,24,by=6),labels=paste(seq(24,48,by=6),seq(48,72,by=6),sep="|"))
    #scale_y_reverse(breaks=seq(0,24,by=6),labels=seq(0,24,by=6))
  p_pca_right<-p_pca_right+xlab("Time [h]")
  
  
  # Lower and upper equilibrium
  if (AddULeq ==T){
    UpLo<-GetWakeSleepEquilibrium(idx_BF,fits,Expr)
    Upl<-predict(res.pca,newdata = t(cbind(UpLo$WakeEq[idx_BF],UpLo$SleepEq[idx_BF])))
    p_pca<-p_pca+geom_vline(xintercept = Upl$coord[,1],color="red",linetype="22")
    p_pca<-p_pca+geom_hline(yintercept = Upl$coord[,2],color="red",linetype="22")
    p_pca_bott<-p_pca_bott + geom_vline(xintercept = Upl$coord[,1],color="red",linetype="22")
    p_pca_right<-p_pca_right + geom_hline(yintercept = Upl$coord[,2],color="red",linetype="22")
    #p_pca<-p_pca+annotate("point",x=Upl$coord[,1],y=Upl$coord[,2],color="red")
  }
  # 
  p_pca<-p_pca+theme(plot.margin = margin(0,0,0,0))
  p_pca_bott<-p_pca_bott+theme(plot.margin = margin(0,0,0,0))
  p_pca_right<-p_pca_right+theme(plot.margin = margin(0,0,0,0))
  # 
  # 
  p_pca_right<-AddTheme(p_pca_right)
  p_pca_bott<-AddTheme(p_pca_bott)
  p_pca_space<-plot_spacer()+theme(plot.margin = unit(c(0,0,0,0), "cm")) #,plot.background = element_rect(fill = "yellow")
  
  p_pca<-p_pca + scale_x_continuous(position = "top") # Move axis to stick panel together
  p_pca<-p_pca + ggtitle(title)
  
  
  layout <- '
  AB
  CD
  '
  gf<-wrap_plots(A = RemoveXaxis(p_pca), C = p_pca_bott, B = RemoveXlab(RemoveYaxis(p_pca_right)) , D = p_pca_space , design = layout,widths = c(1,.3),heights = c(1,.3))

  
  return(gf)
}


PCAplotTimeCourse<-function(Expr,TimeExpr,fits,fitted_full,fitted_circcomp,fitted_swcomp,PC=c(1,2),linedashedHex="11",xlim=c(-120,100),ylim=c(-60,60)){
  LineSize<-.5
  TextSize<-3
  AlphaEllipse<-.2
  
  
  idx_BF<-fits$Gene[fits$BF>exp(1)]
  idx_BF<-idx_BF[!duplicated(idx_BF)]
  
  TimeFactor<-as.factor(paste("T",TimeExpr,sep=""))
  res.pca <- PCA(cbind.data.frame(Expr[,idx_BF],TIME=TimeFactor),
                 scale.unit=T,quali.sup = dim(Expr[,idx_BF])[2]+1,graph=FALSE,ncp = 10)
  
  
  ProjFitC<-predict(res.pca,newdata = t(fitted_full[colnames(Expr[,idx_BF]),])) # Projection of fits
  ProjFitcirc<-predict(res.pca,newdata = t(fitted_circcomp[colnames(Expr[,idx_BF]),]))
  ProjFitsw<-predict(res.pca,newdata = t(fitted_swcomp[colnames(Expr[,idx_BF]),]))
  
  data<-cbind.data.frame(PC1=ProjFitC$coord[240:960,1],PC2=ProjFitC$coord[240:960,2],Time=seq(24,96,by=.1),
                         PC1_circ=ProjFitcirc$coord[240:960,1],PC2_circ=ProjFitcirc$coord[240:960,2],
                         PC1_sw=ProjFitsw$coord[240:960,1],PC2_sw=ProjFitsw$coord[240:960,2])
  
  
  # Print SWrc
  AmpCir1<-max(ProjFitcirc$coord[1:240,1])-min(ProjFitcirc$coord[1:240,1])
  AmpCir2<-max(ProjFitcirc$coord[1:240,2])-min(ProjFitcirc$coord[1:240,2])
  AmpSWr1<-max(ProjFitsw$coord[1:240,1])-min(ProjFitsw$coord[1:240,1])
  AmpSWr2<-max(ProjFitsw$coord[1:240,2])-min(ProjFitsw$coord[1:240,2])
  print(paste("SWrc PC1:",AmpSWr1/(AmpSWr1+AmpCir1)))
  print(paste("SWrc PC2:",AmpSWr2/(AmpSWr2+AmpCir2)))
  
  # PC1
  idx<-data$Time <=48
  gf<-ggplot(aes(x=Time,y=PC1),data=data[idx,])
  
  rectgrey<-"grey80"
  for (i in seq(36,96,by=24)){
    gf<-gf+annotate("rect",xmin=c(i),xmax=c(i+12),ymin=-Inf,ymax=Inf,fill=rectgrey)
  }
  gf<-AddRectSDggplot2(gf,ymax=Inf)
  
  
  gf<-gf+geom_line(size=LineSize,color=colorCode[["Modelfit"]],linetype=linedashedHex)+scale_x_continuous(breaks=seq(24,96,by=6),labels = seq(24,96,by=6) %% 24)
  gf<-gf + annotate("path",x=data$Time[idx],y=data$PC1_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype=linedashedHex)
  gf<-gf + annotate("path",x=data$Time[idx],y=data$PC1_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype=linedashedHex)
  
  idx<-data$Time > 48
  gf<-gf + annotate("path",x=data$Time[idx],y=data$PC1[idx],color=colorCode[["Modelfit"]],size=LineSize,linetype="solid")
  gf<-gf + annotate("path",x=data$Time[idx],y=data$PC1_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype="solid")
  gf<-gf + annotate("path",x=data$Time[idx],y=data$PC1_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype="solid")
  
  gf<-gf+theme_classic()
  
  # PC2
  idx<-data$Time <=48
  gf2<-ggplot(aes(x=Time,y=PC2),data=data[idx,])
  
  rectgrey<-"grey80"
  for (i in seq(36,96,by=24)){
    gf2<-gf2+annotate("rect",xmin=c(i),xmax=c(i+12),ymin=-Inf,ymax=Inf,fill=rectgrey)
  }
  gf2<-AddRectSDggplot2(gf2,ymax=Inf)
  
  
  gf2<-gf2+geom_line(size=LineSize,color=colorCode[["Modelfit"]],linetype=linedashedHex)+
    scale_x_continuous(breaks=seq(24,96,by=6),labels = seq(24,96,by=6))
  gf2<-gf2 + annotate("path",x=data$Time[idx],y=data$PC2_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype=linedashedHex)
  gf2<-gf2 + annotate("path",x=data$Time[idx],y=data$PC2_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype=linedashedHex)
  
  idx<-data$Time > 48
  gf2<-gf2 + annotate("path",x=data$Time[idx],y=data$PC2[idx],color=colorCode[["Modelfit"]],size=LineSize,linetype="solid")
  gf2<-gf2 + annotate("path",x=data$Time[idx],y=data$PC2_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype="solid")
  gf2<-gf2 + annotate("path",x=data$Time[idx],y=data$PC2_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype="solid")
  
  gf2<-gf2+theme_classic()
  gf2<-gf2+xlab("Time [h]")
  
  gf2<-gf2+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  gf<-gf+scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  
  
  return(RemoveXaxis(gf)/gf2)
  
}


GetPCA<-function(Expr,TimeExpr,fits){
  
  library(FactoMineR)
  library(factoextra)
  
  if ("ProbeName" %in% colnames(fits)){
    idx_BF<-fits$ProbeName[fits$BF>exp(1)]
    AddGenes<-T
  }else{
    idx_BF<-fits$Gene[fits$BF>exp(1)]
    AddGenes<-F
  }
  
  idx_BF<-idx_BF[!duplicated(idx_BF)]
  
  TimeFactor<-as.factor(paste("T",TimeExpr,sep=""))
  res.pca <- PCA(cbind.data.frame(Expr[,idx_BF],TIME=TimeFactor),
                 scale.unit=T,quali.sup = dim(Expr[,idx_BF])[2]+1,graph=FALSE,ncp = 10)
  
  if (AddGenes == T){
    res.pca$var$Genes<-fits[rownames(res.pca$var$contrib),"Gene"]
    names(res.pca$var$Genes)<-rownames(res.pca$var$contrib)
  }
  
  return(res.pca)  
}


#########################################################

PCAplot_Human<-function(Expr,TimeExpr,
                        Dataset="Sync",fits,
                        fitted_full,fitted_circcomp,fitted_swcomp,
                        linedashedHex="11",xlim=c(-120,100),ylim=c(-60,60),SWdistr,title){
  
  AddTheme<-function(p){
    p<-p+theme(panel.background=element_rect(fill = 'white', colour = 'black'),panel.grid.major=element_blank(),panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank())
    return(p)
  }
  
  if (Dataset=="Sync"){
    TimeToPlot<-c(44,48,52,56,60,64,68)
    FittedLim<-c(441,720)
    arrowbsl<-c(240)
    arrowexp<-c(640)
    FittedSidePlotLim<-c(481,720)
  }
  
  if (Dataset=="Desync"){
    TimeToPlot<-c(116,120,124,128,132,136,140)
    FittedLim<-c(1161,1440)
    arrowbsl<-c(240)
    arrowexp<-c(1260)
    FittedSidePlotLim<-c(1201,1440)
  }
  if (Dataset == "Ext"){
    TimeToPlot<-c(208,211,214,217,220,223,226,229,232,235)
    FittedLim<-c(2060,2350)
    arrowexp<-c(2220)
    arrowbsl<-c(240)
    FittedSidePlotLim<-c(2161,2400)
  }
  if (Dataset == "Res"){
    TimeToPlot<-c(209,212,215,218,221,227,230,233,236,206,224)
    FittedLim<-c(2060,2350)
    arrowexp<-c(2220)
    arrowbsl<-c(240)
    FittedSidePlotLim<-c(2161,2400)
  }
  
  LineSize<-.5
  TextSize<-3
  AlphaEllipse<-.2
  
  library(ggrepel)
  library(FactoMineR)
  library(factoextra)
  
  idx_BF<-fits$ProbeName[fits$BF>exp(1)]
  idx_BF<-idx_BF[!duplicated(idx_BF)]
  idx_BF<-idx_BF[idx_BF %in% rownames(fitted_full)]
  
  TimeFactor<-as.factor(paste("T",TimeExpr,sep=""))
  res.pca <- PCA(cbind.data.frame(Expr[,idx_BF],TIME=TimeFactor),
                 scale.unit=T,quali.sup = dim(Expr[,idx_BF])[2]+1,graph=FALSE,ncp = 10)
  
  
  ProjFitC<-predict(res.pca,newdata = t(fitted_full[idx_BF,])) # Projection of fits
  ProjFitcirc<-predict(res.pca,newdata = t(fitted_circcomp[idx_BF,]))
  ProjFitsw<-predict(res.pca,newdata = t(fitted_swcomp[idx_BF,]))
  ellipse.coord = coord.ellipse(cbind.data.frame(TimeFactor, res.pca$ind$coord),bary=T,level.conf = .95) # ellipses per time-points
  
  # Print SWrc
  AmpCir1<-max(ProjFitcirc$coord[1:240,1])-min(ProjFitcirc$coord[1:240,1])
  AmpCir2<-max(ProjFitcirc$coord[1:240,2])-min(ProjFitcirc$coord[1:240,2])
  AmpSWr1<-max(ProjFitsw$coord[1:240,1])-min(ProjFitsw$coord[1:240,1])
  AmpSWr2<-max(ProjFitsw$coord[1:240,2])-min(ProjFitsw$coord[1:240,2])
  print(paste("SWrc PC1 [bsl]:",AmpSWr1/(AmpSWr1+AmpCir1)))
  print(paste("SWrc PC2 [bsl]:",AmpSWr2/(AmpSWr2+AmpCir2)))
  
  FittedLim1<-1
  FittedLim2<-nrow(ProjFitcirc$coord)
  AmpCir1<-max(ProjFitcirc$coord[FittedLim1:FittedLim2,1])-min(ProjFitcirc$coord[FittedLim1:FittedLim2,1])
  AmpCir2<-max(ProjFitcirc$coord[FittedLim1:FittedLim2,2])-min(ProjFitcirc$coord[FittedLim1:FittedLim2,2])
  AmpSWr1<-max(ProjFitsw$coord[FittedLim1:FittedLim2,1])-min(ProjFitsw$coord[FittedLim1:FittedLim2,1])
  AmpSWr2<-max(ProjFitsw$coord[FittedLim1:FittedLim2,2])-min(ProjFitsw$coord[FittedLim1:FittedLim2,2])
  print(paste("SWrc PC1 [full time-course]:",AmpSWr1/(AmpSWr1+AmpCir1)))
  print(paste("SWrc PC2 [full time-course]:",AmpSWr2/(AmpSWr2+AmpCir2)))

  
  
  # Plot data under baseline condition (24h sleep-wake cycle)
  data<-cbind.data.frame(PC1=ProjFitC$coord[1:241,1],PC2=ProjFitC$coord[1:241,2])
  p_pca<-ggplot(aes(PC1,PC2),data=data)
  # Add ellipse
  for (i in unique(TimeFactor)){
    Time<-as.numeric(gsub("T","",i))
    if (Time %in% TimeToPlot)
    p_pca<-p_pca + annotate("polygon",x=ellipse.coord$res$Dim.1[ellipse.coord$res$TimeFactor == i],
                            y=ellipse.coord$res$Dim.2[ellipse.coord$res$TimeFactor == i],fill=alpha(GetContZTCols(Time),AlphaEllipse),col=rgb(0,0,0,AlphaEllipse))
  }
  
  # Add data sync (28h sleep-wake cycle)
  p_pca<-p_pca+annotate("path",x=ProjFitC$coord[FittedLim[1]:FittedLim[2],1],
                        y=ProjFitC$coord[FittedLim[1]:FittedLim[2],2],size=LineSize,col=colorCode[["Modelfit"]],linetype="solid")
  # Plot fitted lines
  p_pca<-p_pca+geom_path(col=colorCode[["Modelfit"]],size=LineSize,linetype=linedashedHex)
  
  # # Add text
  # TimeToWrite<-c("T24","T27","T30","T36","T42","T51","T54","T60","T66")
  # idx_TTW<-ellipse.coord$res$TimeFactor %in% TimeToWrite
  # Agg<-aggregate(ellipse.coord$res[idx_TTW,c(2,3)],list(ellipse.coord$res$TimeFactor[idx_TTW]),mean)
  # Agg$Label<-paste("ZT",as.numeric(gsub("T","",Agg$Group.1)) %% 24,sep="")
  # p_pca<-p_pca+ annotate("text",x=Agg$Dim.1,y=Agg$Dim.2,label=Agg$Label,size=TextSize,color=GetContZTCols(as.numeric(gsub("T","",Agg$Group.1))),fontface=2)
  
  
  # Add text
  idx<-ellipse.coord$res$TimeFactor %in% TimeFactor[as.numeric(gsub("T","",TimeFactor)) %in% TimeToPlot]
  Agg<-aggregate(ellipse.coord$res[idx,c(2,3)],list(ellipse.coord$res$TimeFactor[idx]),mean)
  TimeToWrite<-list("T44"=expression(T44[ZT20]),
                    "T48"=expression(T48[ZT0]),
                    "T52"=expression(T52[ZT4]),
                    "T56"=expression(T56[ZT8]),
                    "T60"=expression(T60[ZT12]),
                    "T64"=expression(T64[ZT16]),
                    "T68"=expression(T68[ZT20]),
                    "T116"=expression(T116[ZT20]),
                    "T120"=expression(T120[ZT0]),
                    "T124"=expression(T124[ZT4]),
                    "T128"=expression(T128[ZT8]),
                    "T132"=expression(T132[ZT12]),
                    "T136"=expression(T136[ZT16]),
                    "T140"=expression(T140[ZT20]),
                    "T208"=expression(T208[ZT16]),
                    "T211"=expression(T211[ZT19]),
                    "T214"=expression(T214[ZT22]),
                    "T217"=expression(T217[ZT1]),
                    "T220"=expression(T220[ZT4]),
                    "T223"=expression(T223[ZT7]),
                    "T226"=expression(T226[ZT10]),
                    "T229"=expression(T229[ZT13]),
                    "T232"=expression(T232[ZT16]),
                    "T235"=expression(T235[ZT19]),
                    "T206"=expression(T206[ZT14]),
                    "T209"=expression(T209[ZT17]),
                    "T212"=expression(T212[ZT20]),
                    "T215"=expression(T215[ZT23]),
                    "T218"=expression(T218[ZT2]),
                    "T221"=expression(T221[ZT5]),
                    "T224"=expression(T224[ZT8]),
                    "T227"=expression(T227[ZT11]),
                    "T230"=expression(T230[ZT14]),
                    "T233"=expression(T233[ZT17]),
                    "T236"=expression(T236[ZT20]))
  
  for (i in as.character(Agg$Group.1)){
    p_pca<-p_pca+ annotate("text_repel",x=Agg[Agg$Group.1 == i,"Dim.1"],y=Agg[Agg$Group.1 == i,"Dim.2"],
                           label=TimeToWrite[[i]],
                           size=TextSize,
                           color=GetContZTCols(as.numeric(gsub("T","",i))),
                           fontface=2,
                           box.padding = 0)
  }
  
  
  
  # Add arrows
  p_pca<-p_pca+annotate("segment",x=ProjFitC$coord[arrowbsl,1],y=ProjFitC$coord[arrowbsl,2],
                        xend=ProjFitC$coord[arrowbsl+5,1],yend=ProjFitC$coord[arrowbsl+5,2],
                        col=colorCode[["Modelfit"]], #
                        arrow=arrow(type = "closed",length = unit(.25,"cm")))
  
  p_pca<-p_pca+annotate("segment",x=ProjFitC$coord[arrowexp,1],
                        y=ProjFitC$coord[arrowexp,2],
                        xend=ProjFitC$coord[arrowexp+5,1],
                        yend=ProjFitC$coord[arrowexp+5,2],
                        col=colorCode[["Modelfit"]], #
                        arrow=arrow(type = "closed",length = unit(.25,"cm")))
  
  # Limits
  p_pca<-p_pca+ylab(paste("PC2 [",round(res.pca$eig[2,2]),"%]",sep=""))+
    coord_cartesian(xlim=c(xlim[1],xlim[2]),ylim=c(ylim[1],ylim[2]))#+xlim(xlim[1],xlim[2])+ylim(ylim[1],ylim[2])
  
  
  p_pca<-AddTheme(p_pca)
  
  ########### ADD SIZE PLOTS
  
  library(patchwork)
  
  dd<-cbind.data.frame(Time=seq(.1,24,by=.1),PC2=ProjFitcirc$coord[1:240,2])
  
  p_pca_right<-ggplot(aes(x=Time,y=PC2),data=dd)+
    geom_path(size=LineSize,col=colorCode[["SCNforce"]],linetype="solid")+
    scale_y_continuous(limits=c(ylim[1],ylim[2]),position = "right")
  
  p_pca_right<-p_pca_right+
    annotate("path",x=seq(.1,24,by=.1),y=ProjFitsw$coord[1:240,2],col=colorCode[["SWforce"]],linetype=linedashedHex,size=LineSize)
  
  p_pca_right<-p_pca_right+
    annotate("path",x=seq(.1,24,by=.1),y=ProjFitsw$coord[FittedSidePlotLim[1]:FittedSidePlotLim[2],2],col=colorCode[["SWforce"]],linetype="solid",size=LineSize)
  if (Dataset %in% c("Sync","Desync")){
    p_pca_right<-p_pca_right+scale_x_continuous(breaks=seq(0,24,by=6),
                       labels=paste(seq(0,24,by=6),seq(TimeToPlot[-1][1],TimeToPlot[-1][1]+24,by=6),sep="|"),
                       position = "top")+
    theme(axis.text.x.top = element_text(angle=90,vjust=0.5))
  }else{
    p_pca_right<-p_pca_right+scale_x_continuous(breaks=seq(0,24,by=6),
                                                labels=paste(seq(0,24,by=6),seq(216,240,by=6),sep="|"),
                                                position = "top")+
      theme(axis.text.x.top = element_text(angle=90,vjust=0.5))
  }
  #p_pca+ylim(ylim[1],ylim[2])|p_pca_right
  
  dd<-cbind.data.frame(Time=seq(.1,24,by=.1),PC1=ProjFitsw$coord[1:240,1])
  p_pca_bott<-ggplot(aes(x=PC1,y=Time),data=dd)+
    geom_path(col=colorCode[["SWforce"]],linetype=linedashedHex,size=LineSize)+
    annotate("path",x=ProjFitsw$coord[FittedSidePlotLim[1]:FittedSidePlotLim[2],1],y=seq(.1,24,by=.1),col=colorCode[["SWforce"]],linetype="solid",size=LineSize)+
    xlim(xlim[1],xlim[2])+
    annotate("path",x=ProjFitcirc$coord[FittedSidePlotLim[1]:FittedSidePlotLim[2],1],y=seq(.1,24,by=.1),col=colorCode[["SCNforce"]],size=LineSize)
  print(paste(seq(0,24,by=6),seq(TimeToPlot[-1][1],TimeToPlot[-1][1]+24,by=6),sep="|"))
  #(p_pca_bott+xlim(xlim[1],xlim[2]))/(p_pca_bott)
  
  p_pca_bott<-p_pca_bott+
    ylab("Time [h]")+
    xlab(paste("PC1 [",round(res.pca$eig[1,2]),"%]",sep=""))
  if (Dataset %in% c("Sync","Desync")){
    p_pca_bott<-p_pca_bott+scale_y_reverse(breaks=seq(0,24,by=6),
                    labels=paste(seq(0,24,by=6),seq(TimeToPlot[-1][1],TimeToPlot[-1][1]+24,by=6),sep="|"))
  }else{
    p_pca_bott<-p_pca_bott+scale_y_reverse(breaks=seq(0,24,by=6),
                                           labels=paste(seq(0,24,by=6),seq(216,240,by=6),sep="|"))   
  }
    
  p_pca_right<-p_pca_right+xlab("Time [h]")
  

  
  p_pca<-p_pca+theme(plot.margin = margin(0,0,0,0))
  p_pca_bott<-p_pca_bott+theme(plot.margin = margin(0,0,0,0))
  p_pca_right<-p_pca_right+theme(plot.margin = margin(0,0,0,0))
  
  p_pca_right<-AddTheme(p_pca_right)
  p_pca_bott<-AddTheme(p_pca_bott)
  
  p_pca_space<-plot_spacer()+theme(plot.margin = unit(c(0,0,0,0), "cm")) #,plot.background = element_rect(fill = "yellow")
  p_pca<-p_pca + 
    scale_x_continuous(position = "top") # Move axis to stick panel together  
  p_pca<-p_pca + ggtitle(title)
  
  
  
  layout <- '
  AB
  CD
  '
  gf<-wrap_plots(A = RemoveXaxis(p_pca), C = p_pca_bott, B = RemoveXlab(RemoveYaxis(p_pca_right)),D=p_pca_space , design = layout,widths = c(1,.3),heights = c(1,.3))
  
  # Time course plot
  dd<-cbind.data.frame(Time=seq(0,(length(ProjFitC$coord[,1])-1)/10,by=.1),PC1=ProjFitC$coord[,1],
                       PC1_sw=ProjFitsw$coord[,1],PC1_circ=ProjFitcirc$coord[,1],
                       PC2=ProjFitC$coord[,2],
                       PC2_sw=ProjFitsw$coord[,2],PC2_circ=ProjFitcirc$coord[,2])
  idx<-dd$Time<44
  ts_plot_1<-ggplot(aes(x=Time,y=PC1),data=dd[idx,])
  startsleep<-SWdistr$Time[SWdistr$Sleep>0][which(diff(c(0,SWdistr$Time[SWdistr$Sleep>0]))>1)]
  endsleep<-c(SWdistr$Time[SWdistr$Sleep>0][which(diff(c(0,SWdistr$Time[SWdistr$Sleep>0]))>1)-1],rev(SWdistr$Time[SWdistr$Sleep>0])[1])
  startsleep<-c(startsleep[1]-24,startsleep)
  endsleep<-c(endsleep[1]-24,endsleep)
  ts_plot_1 <- ts_plot_1 + annotate("rect", xmin = startsleep, xmax = endsleep, 
                      ymin = -Inf, ymax = Inf, 
                      fill="grey75",alpha = 0.5)
  
  ts_plot_1<-ts_plot_1+geom_line(size=LineSize,color=colorCode[["Modelfit"]],linetype=linedashedHex)
  ts_plot_1<-ts_plot_1 + annotate("path",x=dd$Time[idx],y=dd$PC1_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype=linedashedHex)
  ts_plot_1<-ts_plot_1 + annotate("path",x=dd$Time[idx],y=dd$PC1_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype=linedashedHex)
  
  idx<-dd$Time >= 44 
  ts_plot_1<-ts_plot_1 + annotate("path",x=dd$Time[idx],y=dd$PC1[idx],color=colorCode[["Modelfit"]],size=LineSize)
  ts_plot_1<-ts_plot_1 + annotate("path",x=dd$Time[idx],y=dd$PC1_circ[idx],color=colorCode[["SCNforce"]],size=LineSize)
  ts_plot_1<-ts_plot_1 + annotate("path",x=dd$Time[idx],y=dd$PC1_sw[idx],color=colorCode[["SWforce"]],size=LineSize)
  ts_plot_1<-ts_plot_1 + geom_vline(xintercept = seq(0,225,by=24))
  ts_plot_1<-ts_plot_1+theme_classic()
  
  if (Dataset %in% c("Sync","Desync")){
    ts_plot_1<-ts_plot_1+scale_x_continuous(breaks=seq(12,228,by=6),labels = seq(12,228,by=6) %% 24,limits = c(12,228))
    # FD line
    ts_plot_1<-ts_plot_1+coord_cartesian(ylim = xlim,clip = "off",expand = F)
    ts_plot_1<-ts_plot_1+annotate("line",x=c(44,68),y=c(xlim[2]-1,xlim[2]-1),size=.5,color="red",linetype="solid")
    ts_plot_1<-ts_plot_1+annotate("text",label="in-phase",x=56,y=xlim[2],vjust=0,color="red",fontface = "bold")
    ts_plot_1<-ts_plot_1+annotate("line",x=c(116,140),y=c(xlim[2]-1,xlim[2]-1),size=.5,color="red",linetype="solid")
    ts_plot_1<-ts_plot_1+annotate("text",label="anti-phase",x=128,y=xlim[2],vjust=0,color="red",fontface = "bold")
    
    
  }else{
    seqbreak<-c(seq(0,96,by=12),seq(204,258,by=12))
    ts_plot_1<-ts_plot_1+scale_x_continuous(breaks=seqbreak,labels = seqbreak %% 24,trans = squish_trans(90,192, 50),limits = c(12,258))
    # Hide squish
    ts_plot_1<-ts_plot_1 + annotate("rect",xmin=90,xmax=192, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
    
    # CR line
    ts_plot_1<-ts_plot_1+coord_cartesian(ylim = xlim,clip = "off",expand = F)
    ts_plot_1<-ts_plot_1+annotate("line",x=c(206,240),y=c(xlim[2]-1,xlim[2]-1),size=.5,color="red",linetype="solid")
    ts_plot_1<-ts_plot_1+annotate("text",label="C.R.",x=221,y=xlim[2],vjust=0,color="red",fontface = "bold")
  }
  
  
  idx<-dd$Time<44
  ts_plot_2<-ggplot(aes(x=Time,y=PC2),data=dd[idx,])
  ts_plot_2 <- ts_plot_2 + annotate("rect", xmin = startsleep, xmax = endsleep, 
                                    ymin = -Inf, ymax = Inf, 
                                    fill="grey75",alpha = 0.5)
  
  ts_plot_2<-ts_plot_2+geom_line(size=LineSize,color=colorCode[["Modelfit"]],linetype=linedashedHex)
  ts_plot_2<-ts_plot_2 + annotate("path",x=dd$Time[idx],y=dd$PC2_circ[idx],color=colorCode[["SCNforce"]],size=LineSize,linetype=linedashedHex)
  ts_plot_2<-ts_plot_2 + annotate("path",x=dd$Time[idx],y=dd$PC2_sw[idx],color=colorCode[["SWforce"]],size=LineSize,linetype=linedashedHex)
  
  idx<-dd$Time >= 44 
  ts_plot_2<-ts_plot_2 + annotate("path",x=dd$Time[idx],y=dd$PC2[idx],color=colorCode[["Modelfit"]],size=LineSize)
  ts_plot_2<-ts_plot_2 + annotate("path",x=dd$Time[idx],y=dd$PC2_circ[idx],color=colorCode[["SCNforce"]],size=LineSize)
  ts_plot_2<-ts_plot_2 + annotate("path",x=dd$Time[idx],y=dd$PC2_sw[idx],color=colorCode[["SWforce"]],size=LineSize)
  ts_plot_2<-ts_plot_2 + geom_vline(xintercept = seq(0,225,by=24))
  ts_plot_2<-ts_plot_2+theme_classic()
  
  if (Dataset %in% c("Sync","Desync")){
    ts_plot_2<-ts_plot_2+scale_x_continuous(breaks=seq(12,228,by=6),labels = seq(12,228,by=6),limits = c(12,228))
    ts_plot_2<-ts_plot_2+coord_cartesian(ylim = ylim,clip = "off",expand = F)
  }else{
    seqbreak<-c(seq(0,96,by=12),seq(204,258,by=12))
    ts_plot_2<-ts_plot_2+scale_x_continuous(breaks=seqbreak,labels = seqbreak,trans = squish_trans(90,192, 50),limits = c(12,258))
    # Hide squish
    ts_plot_2<-ts_plot_2 + annotate("rect",xmin=90,xmax=192, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
    ts_plot_2<-ts_plot_2+coord_cartesian(ylim = ylim,clip = "off",expand = F)
  }
  
  ts_plot_1<-ts_plot_1+theme(plot.margin = margin(10,0,0,0))
  
  
  ts_plot_1<-ts_plot_1+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
  
  ts_plot_2<-ts_plot_2+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
    xlab("Time [h]")
  
  RemoveXaxis(ts_plot_1)/ts_plot_2
  
  return(list(pca_plot=gf,TimeCoursePlot=RemoveXaxis(ts_plot_1)/ts_plot_2))
  
}





#################################
PlotParams<-function(res.pca,SWrcs,data,fits,PC=1){
  
  # PCA contributors
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  Getomega<-sapply(ss,function(x){
    fits$omega[fits$Gene %in% x][1]
  })
  GetPhiC<-sapply(ss,function(x){
    Amp<-fits$AmpSin[fits$Gene %in% x][1]
    Phi<-fits$PhiSin[fits$Gene %in% x][1]
    if (Amp<0){Phi<-Phi+pi}
    return(Phi %% (2*pi))
  })
  Getloggamma<-sapply(ss,function(x){
    fits$loggamma[fits$Gene %in% x][1]
  })
  dd<-cbind.data.frame(Zeta=log10(GetZeta(ss,data$fits)[ss]),
                       omega=Getomega[ss],
                       upper=GetUpperLowerEquilibrium(ss,fits = data$fits,data$rna_expr,asSD = T)$upper,
                       pl=GetPhaseLag(ss,data$fits,PosAmp=T)[ss],
                       PCA_Contrib=res.pca$var$contrib[ss,PC],
                       SWrc=SWrcs[ss],
                       Phase_Circ=GetPhiC[ss],
                       loggamma=Getloggamma[ss])
  
  dd<-dd[order(dd$PCA_Contrib),]
  dd$Damped<-rep(NA)
  dd$Damped[dd$Zeta>0]<-"Overdamped"
  dd$Damped[dd$Zeta<0]<-"Underdamped"
  
  alpha<-.75
  # gg1<-ggplot(aes(pl,Zeta,fill=SWrc,size=PCA_Contrib),data=dd)+geom_point(shape=21)+scale_fill_gradient2(low=colorCode[["SCNforce"]],midpoint = .5,high=colorCode[["SWforce"]],mid="forestgreen")+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)))+ylab(expression(paste("log10",zeta)))+scale_size_continuous(range = c(0,3))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))
  gg<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Overdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+ggtitle(expression(zeta > 1 ))+labs(size="PC cumulative\ncontribution [%]")+ylim(0,0.55)
  gg2<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Underdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+theme(legend.position = "none")+ggtitle(expression(zeta < 1 ))+ylim(0,0.55)
  return(list(plot=gg2|gg,data=dd))
}

PlotParams2<-function(res.pca,SWrcs,data,fits,PC=1){
  
  require(ggrastr)
  
  # PCA contributors
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  Getomega<-sapply(ss,function(x){
    fits$omega[fits$Gene %in% x][1]
  })
  GetWake<-sapply(ss,function(x){
    fits$Wake[fits$Gene %in% x][1]
  })
  GetPhiC<-sapply(ss,function(x){
    Amp<-fits$AmpSin[fits$Gene %in% x][1]
    Phi<-fits$PhiSin[fits$Gene %in% x][1]
    if (Amp<0){Phi<-Phi+pi}
    return(Phi %% (2*pi))
  })
  Getloggamma<-sapply(ss,function(x){
    fits$loggamma[fits$Gene %in% x][1]
  })
  dd<-cbind.data.frame(Zeta=log10(GetZeta(ss,data$fits)[ss]),
                       omega=Getomega[ss],
                       upper=GetUpperLowerEquilibrium(ss,fits = data$fits,data$rna_expr,asSD = T)$upper,
                       pl=GetPhaseLag(ss,data$fits,PosAmp=T)[ss],
                       PCA_Contrib=res.pca$var$contrib[ss,PC],
                       SWrc=SWrcs[ss],
                       Phase_Circ=GetPhiC[ss],
                       loggamma=Getloggamma[ss],
                       Wake=scale(GetWake[ss]))
  
  dd<-dd[order(dd$PCA_Contrib,decreasing = T),]
  dd$Damped<-rep(NA)
  dd$Damped[dd$Zeta>0]<-"Overdamped"
  dd$Damped[dd$Zeta<0]<-"Underdamped"
  dd$WakeSign<-rep(NA)
  dd$WakeSign[dd$Wake>0]<-"Positive"
  dd$WakeSign[dd$Wake<0]<-"Negative"
  
  loggammaseq<-seq(-4.6,3,by=.001)
  w<-2*pi/24 # freq driving C
  
  alpha<-.75
  # gg1<-ggplot(aes(pl,Zeta,fill=SWrc,size=PCA_Contrib),data=dd)+geom_point(shape=21)+scale_fill_gradient2(low=colorCode[["SCNforce"]],midpoint = .5,high=colorCode[["SWforce"]],mid="forestgreen")+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)))+ylab(expression(paste("log10",zeta)))+scale_size_continuous(range = c(0,3))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))
  gg<-ggplot(aes(loggamma,omega,fill=SWrc,size=PCA_Contrib,shape=WakeSign),data=dd)+#ggplot(aes(loggamma,omega,fill=SWrc,size=PCA_Contrib,shape=Damped),data=dd)+
    geom_point(color=rgb(0,0,0,.55))+
    scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+
    theme_bw()+
    xlab(expression(paste(log,gamma)))+
    ylab(expression(paste(omega[0])))+
    scale_size_continuous(range = c(5,0))+
    scale_shape_manual(values=c("Positive"=21,"Negative"=24))+
    # scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+
    labs(size="PC cumulative\ncontribution [%]")+
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/2,color="black")+ # zeta = 1
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/1,color="grey50")+ # zeta = 0.5
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/.2,color="grey50")+ # zeta = 0.1
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/4,color="grey50")+ # Zeta = 2
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/20,color="grey50")+ # Zeta = 10
    annotate("path",x=loggammaseq,y=rep(w,length(loggammaseq)),linetype="22",col="red")+ # phase-lag = pi/2
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/2.666) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/4) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/8) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/16) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.715)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.5)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.25)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.125)) + w^2)  ,linetype="22",col="red")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    coord_cartesian(xlim = c(-4.6,3),ylim=c(0.05,0.55))
  
  
  dd$Phase_Osc<-(dd$Phase_Circ-dd$pl) %% (2*pi)
  dd$PhaseWake<-((pi-dd$pl)) %% (2*pi)
  dd$PhaseWake[dd$Wake<0]<- (dd$PhaseWake[dd$Wake<0] + pi) %% (2*pi)
  dd$PhaseDiff<-atan2(sin(dd$Phase_Osc-dd$PhaseWake), cos(dd$Phase_Osc-dd$PhaseWake))
  
  gg2<-ggplot(aes(Phase_Osc,PhaseDiff,fill=SWrc,size=PCA_Contrib,shape=WakeSign),data=dd)+#ggplot(aes(Phase_Osc,PhaseDiff,fill=SWrc,size=PCA_Contrib,shape=Damped),data=dd)+
    scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+
    geom_point(color=rgb(0,0,0,.45))+
    scale_shape_manual("Wake effect sign",values=c("Positive"=21,"Negative"=24))+
    #scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+
    coord_polar()+
    theme_bw()+    
    labs(size="PC cumulative\ncontribution [%]")+
    scale_x_continuous(breaks = c(0,pi/2,pi,3/4*2*pi,2*pi),limits = c(0,2*pi),labels=c(0,expression(pi/2),expression(pi),expression(3/2*pi),0))+
    scale_size_continuous(range = c(5,0))+
    scale_y_continuous(breaks = c(-pi,0,pi),labels = c(expression(-pi),0,expression(pi)),limits = c(-5,pi))+
    geom_hline(yintercept = 0)+
    ylab(expression(paste(Delta,phi,"(Cr-SWr)")))+
    xlab(expression(paste(phi,"-Circadian response (Cr)")))
  
  #    ggtitle(expression(zeta > 1 ))+
  #gg2<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Underdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+theme(legend.position = "none")+ggtitle(expression(zeta < 1 ))
  
  gg<-rasterize(gg, layers='Point', dpi=200)
  gg2<-rasterize(gg2, layers='Point', dpi=200)
  
  return(list(plot1=gg,plot2=gg2,data=dd))
}


PlotParamsH2<-function(res.pca,SWrcs,fits,PC=1){
  
  require(ggrastr)
  
  # PCA contributors
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  
  Getomega<-sapply(ss,function(x){
    fits$omega[fits$ProbeName %in% x][1]
  })
  GetWake<-sapply(ss,function(x){
    fits$Wake[fits$ProbeName %in% x][1]
  })
  GetPhiC<-sapply(ss,function(x){
    Amp<-fits$AmpSin[fits$ProbeName %in% x][1]
    Phi<-fits$PhiSin[fits$ProbeName %in% x][1]
    if (Amp<0){Phi<-Phi+pi}
    return(Phi %% (2*pi))
  })
  Getloggamma<-sapply(ss,function(x){
    fits$loggamma[fits$ProbeName %in% x][1]
  })
  dd<-cbind.data.frame(Zeta=log10(GetZeta(ss,fits)[ss]),
                       omega=Getomega[ss],
                       pl=GetPhaseLag(ss,fits,PosAmp=T)[ss],
                       PCA_Contrib=res.pca$var$contrib[ss,PC],
                       SWrc=SWrcs[ss],
                       Phase_Circ=GetPhiC[ss],
                       loggamma=Getloggamma[ss],
                       Wake=scale(GetWake[ss]))
  
  dd<-dd[order(dd$PCA_Contrib,decreasing = T),]
  dd$Damped<-rep(NA)
  dd$Damped[dd$Zeta>0]<-"Overdamped"
  dd$Damped[dd$Zeta<0]<-"Underdamped"
  dd$WakeSign<-rep(NA)
  dd$WakeSign[dd$Wake>0]<-"Positive"
  dd$WakeSign[dd$Wake<0]<-"Negative"
  
  
  loggammaseq<-seq(-4.6,3,by=.001)
  w<-2*pi/24 # freq driving C
  
  alpha<-.75
  
  gg<-ggplot(aes(loggamma,omega,fill=SWrc,size=PCA_Contrib,shape=WakeSign),data=dd)+#ggplot(aes(loggamma,omega,fill=SWrc,size=PCA_Contrib,shape=Damped),data=dd)+
    geom_point(color=rgb(0,0,0,.45))+
    scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+
    theme_bw()+
    xlab(expression(paste(log,gamma)))+
    ylab(expression(paste(omega[0])))+
    scale_size_continuous(range = c(5,0))+
    scale_shape_manual("Wake effect sign",values=c("Positive"=21,"Negative"=24))+
    #scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+
    labs(size="PC cumulative\ncontribution [%]")+
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/2,color="black")+ # zeta = 1
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/1,color="grey50")+ # zeta = 0.5
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/.2,color="grey50")+ # zeta = 0.1
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/4,color="grey50")+ # Zeta = 2
    annotate("path",x=loggammaseq,y=exp(loggammaseq)/20,color="grey50")+ # Zeta = 10
    annotate("path",x=loggammaseq,y=rep(w,length(loggammaseq)),linetype="22",col="red")+ # phase-lag = pi/2
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/2.666) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/4) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/8) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan(pi/16) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.715)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.5)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.25)) + w^2)  ,linetype="22",col="red")+
    annotate("path",x = loggammaseq, y = sqrt((exp(loggammaseq)*w)/tan((pi/1.125)) + w^2)  ,linetype="22",col="red")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    coord_cartesian(xlim = c(-4.6,3.5),ylim=c(0,0.7))
  
  
  dd$Phase_Osc<-(dd$Phase_Circ-dd$pl) %% (2*pi)
  dd$PhaseWake<-((pi-dd$pl)) %% (2*pi)
  dd$PhaseWake[dd$Wake<0]<- (dd$PhaseWake[dd$Wake<0] + pi) %% (2*pi)
  dd$PhaseDiff<-atan2(sin(dd$Phase_Osc-dd$PhaseWake), cos(dd$Phase_Osc-dd$PhaseWake))
  
  gg2<-ggplot(aes(Phase_Osc,PhaseDiff,fill=SWrc,size=PCA_Contrib,shape=WakeSign),data=dd)+#ggplot(aes(Phase_Osc,PhaseDiff,fill=SWrc,size=PCA_Contrib,shape=Damped),data=dd)+
    scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+
    geom_point(color=rgb(0,0,0,.55))+
    scale_shape_manual(values=c("Positive"=21,"Negative"=24))+
    #scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+
    coord_polar()+
    theme_bw()+    
    labs(size="PC cumulative\ncontribution [%]")+
    scale_x_continuous(breaks = c(0,pi/2,pi,3/4*2*pi,2*pi),limits = c(0,2*pi),labels=c(0,expression(pi/2),expression(pi),expression(3/2*pi),0))+
    scale_size_continuous(range = c(5,0))+
    scale_y_continuous(breaks = c(-pi,0,pi),labels = c(expression(-pi),0,expression(pi)),limits = c(-5,pi))+
    geom_hline(yintercept = 0)+
    ylab(expression(paste(Delta,phi,"(Cr-SWr)")))+
    xlab(expression(paste(phi,"-Circadian response (Cr)")))
  
  #    ggtitle(expression(zeta > 1 ))+
  #gg2<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Underdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+theme(legend.position = "none")+ggtitle(expression(zeta < 1 ))
  
  gg<-rasterize(gg, layers='Point', dpi=200)
  gg2<-rasterize(gg2, layers='Point', dpi=200)
  
  return(list(plot1=gg,plot2=gg2,data=dd))  

}




PlotParamsH<-function(res.pca,SWrcs,fits,PC=1){
  
  # PCA contributors
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  
  
  Getomega<-fits[ss,"omega"]
  GetPhiC<-fits[ss,"PhiSin"]
  GetPhiC[fits[ss,"AmpSin"]<0]<-(GetPhiC[fits[ss,"AmpSin"]<0]+pi) %% (2*pi)
  Getloggamma<-fits[ss,"loggamma"]
  dd<-cbind.data.frame(Zeta=log10(GetZeta(ss,fits)[ss]),
                       omega=Getomega,
                       pl=GetPhaseLag(ss,fits,PosAmp=T)[ss],
                       PCA_Contrib=res.pca$var$contrib[ss,PC],
                       SWrc=SWrcs[ss])
  
  dd<-dd[order(dd$PCA_Contrib),]
  dd$Damped<-rep(NA)
  dd$Damped[dd$Zeta>0]<-"Overdamped"
  dd$Damped[dd$Zeta<0]<-"Underdamped"
  
  alpha<-.75
  # gg1<-ggplot(aes(pl,Zeta,fill=SWrc,size=PCA_Contrib),data=dd)+geom_point(shape=21)+scale_fill_gradient2(low=colorCode[["SCNforce"]],midpoint = .5,high=colorCode[["SWforce"]],mid="forestgreen")+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)))+ylab(expression(paste("log10",zeta)))+scale_size_continuous(range = c(0,3))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))
  gg<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Overdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+ggtitle(expression(zeta > 1 ))+labs(size="PC cumulative\ncontribution [%]")
  
  gg2<-ggplot(aes(pl,omega,fill=SWrc,size=PCA_Contrib),data=dd[dd$Damped=="Underdamped",])+geom_point(shape=21)+scale_fill_gradient2(low=alpha(colorCode[["SCNforce"]],alpha),midpoint = .5,high=alpha(colorCode[["SWforce"]],alpha),mid=alpha("forestgreen",alpha),limits=c(0,1))+theme_bw()+xlab(expression(paste(phi,"-lag")))+scale_x_continuous(breaks=c(0,pi/2,pi),labels=c(0,expression(pi/2),expression(pi)),limits = c(0,pi))+ylab(expression(paste(omega[0])))+scale_size_continuous(range = c(5,0))+scale_shape_manual(values=c("Underdamped"=21,"Overdamped"=24))+theme(legend.position = "none")+ggtitle(expression(zeta < 1 ))
  return(list(plot=gg2|gg,data=dd))
}


AmplitdueDiffBetweenPC<-function(res.pca,rna_expr,human=F){
  
  # PC1
  PC<-1
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  TimeThresh<-48
  if (human == T){
    TimeThresh<-68
  }
  
  AmpPC1<-sapply(ss,function(x){

    EE<-rna_expr[rna_expr$Time<=TimeThresh,x]
    Time<-rna_expr[rna_expr$Time<=TimeThresh,"Time"]
    modelcos<-lm(EE~cos(2*pi/24*Time)+sin(2*pi/24*Time))
    amp<-sqrt(coef(modelcos)[2]^2+coef(modelcos)[3]^2)[[1]]
    return(amp)
  })
  
  # PC1
  PC<-2
  ss<-names(sort(res.pca$var$contrib[,PC],decreasing = T))
  res.pca$var$contrib<-res.pca$var$contrib[ss,]
  res.pca$var$contrib[,PC]<-cumsum(res.pca$var$contrib[,PC])
  ss<-ss[1:which.min(abs(res.pca$var$contrib[,PC] - 25))] # top explain 50%
  
  AmpPC2<-sapply(ss,function(x){
    EE<-rna_expr[rna_expr$Time<=TimeThresh,x]
    Time<-rna_expr[rna_expr$Time<=TimeThresh,"Time"]
    modelcos<-lm(EE~cos(2*pi/24*Time)+sin(2*pi/24*Time))
    amp<-sqrt(coef(modelcos)[2]^2+coef(modelcos)[3]^2)[[1]]
    return(amp)
  })
  
  return(list(AmpPC1=AmpPC1,AmpPC2=AmpPC2))
  
}















