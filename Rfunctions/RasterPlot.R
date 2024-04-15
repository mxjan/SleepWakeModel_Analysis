####################################
# SOME FUNCTIONS TO PLOT FIT FOR GENE
# AND PROBES
####################################

################ GENERAL RASTER PLOT FUNCTION ###################
RasterPlot<-function(ProbeGene,Dataset,# Cortex,Liver,Desync,Res,Ctr
                     DaysShows=seq(1,10),ReplicatesSamples=F,PreviousHourShow=6, # Plot option
                     npretty=3,ylim=NULL,xlim_max=24,
                     MeanPointSize=1,AddCISegment=T, # Expression points
                     ModelFitLine=1,FitMethod="Solve",AddASkeldonSol=F,BSL_RepLine=0, # model fits
                     AddFits=NULL,ColsAddFits=NULL,LtyAddFits=NULL,AddLowerUpperEquilibrium=F,AddSleepWake=F){
  
  require(ggplot2)
  require(SWDMr)
  library(scales)
  library(dplyr)
  
  # Load data given dataset
  data<-LoadData(Dataset)
  # Color Codes 
  cols<-ColorCode()
  
  # Get Model
  model<-GetModel(Dataset,data,ProbeGene)
  
  # Get Confidence Interval and mean time point of expression
  CI_Expr_Raster<-GetCIRaster(Dataset,data,ProbeGene,DaysShows,ylim,PreviousHourShow,ReplicatesSamples)
  CI_expr<-CI_Expr_Raster$CI_expr
  rescale_range<-CI_Expr_Raster$rescale_range
  CI_Expr_Raster<-CI_Expr_Raster$CI_expr_raster
  

  # Get Model fitted line(s)
  Model_outputs<-GetModelOutput(model = model,Dataset = Dataset,data = data,
                                ProbeGene =  ProbeGene,FitMethod = FitMethod,
                                AddASkeldonSol = AddASkeldonSol,rescale_range = rescale_range,
                                DaysShows = DaysShows,PreviousHourShow = PreviousHourShow)
  
  # Given other fits ?
  if (! is.null(AddFits)){
    for (i in 1:length(AddFits)){
      dd<-cbind.data.frame(Time=AddFits[[i]]$time,y1=AddFits[[i]]$y1)
      AddFits[[i]]<-ProcessDataRaster(dd,value=c("y1"),DoReplicate=T,Dorescale = T,rescalefrom = rescale_range,replastH =PreviousHourShow)
      AddFits[[i]]<-AddFits[[i]][AddFits[[i]]$Day_axis %in% DaysShows,]
    }
  }
  
  # Data for Plot
  dfd<-as.data.frame(cbind(CI_Expr_Raster$ZT,CI_Expr_Raster[,"Mean"],CI_Expr_Raster$Day_axis))
  colnames(dfd)<-c("ZTime","Gene","Day")
  
  ####### PLOT ########
  p_Gene <- ggplot(dfd,aes(x=ZTime,y=Gene,group=Day))
  p_Gene<-p_Gene+scale_x_continuous(breaks=seq(-6,24,by=6),labels=seq(-6,24,by=6))
  
  # Add LD block given experiments
  p_Gene<-AddLDblock(p_Gene,Dataset,DaysShows,PreviousHourShow)
  
  # Add some lines
  p_Gene<-p_Gene+geom_vline(xintercept=0,color="grey60",size=0.25,linetype="21")
  #p_Gene<-p_Gene+geom_vline(xintercept=12,color="grey60",size=0.25,linetype="21")
  #p_Gene<-p_Gene+geom_vline(xintercept=c(6,18),color="grey85",size=0.25,linetype="23")
  #p_Gene<-p_Gene+geom_vline(xintercept=24,color="grey15",size=0.25,linetype="solid")# line ZT 0
  #p_Gene<-p_Gene+geom_hline(yintercept=unique(c(DaysShows*-1,(DaysShows-1)*-1)),linetype="solid",color="grey10",size=.25) # Separation by day
  
  # theme
  p_Gene<-p_Gene+ theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  
  # Add Expression on Y-axis
  pretty_fun<-scales::breaks_pretty(n=npretty)
  Expr_labels<-pretty_fun(CI_expr$Mean)
  Expr_rescaled<-rescale(Expr_labels,rescale_range,to=c(0.1,0.9))
  
  # if (any(sign(Expr_rescaled) == -1 | Expr_rescaled > 1)){
  #   idx<-which(sign(Expr_rescaled) == -1 | Expr_rescaled > 1)
  #   Expr_rescaled<-Expr_rescaled[-idx]
  #   Expr_labels<-Expr_labels[-idx]
  # }
  Expr_rescaled_raster<-as.vector(sapply(DaysShows,function(x){return(Expr_rescaled-(x))}))
  Expr_labels_raster<-rep(Expr_labels,length(DaysShows))
  
  
  # Squish Y-axis for day without data points
  #p_Gene <- p_Gene + scale_y_continuous(trans = squish_trans(-9, -5, 30),breaks = Expr_rescaled-2,labels=Expr_labels)
  Rdays<-seq(-1*range(DaysShows)[1],-1*range(DaysShows)[2])
  idx_Rdays<-! Rdays %in% (-1*(DaysShows))
  if (any(idx_Rdays)){
    #p_Gene <- p_Gene + scale_y_continuous(trans = squish_day(Rdays,idx_Rdays),breaks = Expr_rescaled-min(DaysShows),labels=Expr_labels) # Add y-axis labels only for day 1
    p_Gene <- p_Gene + scale_y_continuous(trans = squish_day(Rdays,idx_Rdays),breaks = Expr_rescaled_raster,labels=Expr_labels_raster) # Add y-axis labels on all days
    p_Gene<-p_Gene + geom_hline(yintercept = Expr_rescaled_raster, color=rgb(0,0,0,.1),linetype="solid",size=0.25)
  }else{
    p_Gene <- p_Gene + scale_y_continuous(breaks =  Expr_rescaled-min(DaysShows),labels=Expr_labels)
    p_Gene<-p_Gene + geom_hline(yintercept = Expr_rescaled_raster, color=rgb(0,0,0,.1),linetype="solid",size=0.25)
  }
  
  # Add fits if given
  if (! is.null(AddFits)){
    for (i in 1:length(AddFits)){
      if (is.null(LtyAddFits)){lty<-"solid"}else{lty<-LtyAddFits[i]}
      p_Gene<-p_Gene + annotate("path",x=AddFits[[i]]$ZT,y=AddFits[[i]]$y1,
                                group=AddFits[[i]]$Day_axis,color=ColsAddFits[i],size=.75,linetype=LtyAddFits[[i]])
    }
  }
  
  # Given other fits ?
  if (AddSleepWake == T){
    SW<-model@SWdist
    SW<-SW %>% group_by(ceiling(Time)) %>% mutate_at("Wake",mean)
    SW$Wake<-SW$Wake * (0.6 - 0.01) + 0.01
    dd<-cbind.data.frame(Time=SW$Time,y1=SW$Wake)
    dd<-ProcessDataRaster(dd,value=c("y1"),DoReplicate=T,Dorescale = F,rescalefrom = rescale_range,replastH =PreviousHourShow)
    dd<-dd[dd$Day_axis %in% DaysShows,]
    for (i in DaysShows){
      p_Gene<-p_Gene +
        annotate("ribbon",x=dd$ZT[dd$Day_axis==i],ymax=dd$y1[dd$Day_axis==i],ymin=floor(min(dd$y1[dd$Day_axis==i])),
                                fill=alpha(cols[["SWforce"]],.2))+
        annotate("line",x=dd$ZT[dd$Day_axis==i],y=dd$y1[dd$Day_axis==i],color=alpha(cols[["SWforce"]],.7))
    }

  }
  
  # BSL replicated line
  if (BSL_RepLine != 0){
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_bsl$ZT,y=Model_outputs$out_response_bsl$y1,
                              group=Model_outputs$out_response_bsl$Day_axis,color=alpha("black",.65),size=BSL_RepLine,linetype="42")
  }
  
  # Fitted lines
  if (ModelFitLine != 0 & AddASkeldonSol ==F){
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response$ZT,y=Model_outputs$out_response$y1,
                              group=Model_outputs$out_response$Day_axis,color="black",size=ModelFitLine,alpha=1)
  }
  
  # Add line for Anne Skeldon solution
  if (AddASkeldonSol ==T){
    if (Dataset %in% c("Cortex","Liver","Desync")){
      ddint<-as.data.frame(cbind("Time"=12+((DaysShows-1)*24),"Intercept"=rep(model@intercept,length(DaysShows))))
    }else{
      intercept<-data$fits[data$fits$ProbeName == ProbeGene,"intercept"]
      ddint<-as.data.frame(cbind("Time"=12+((DaysShows-1)*24),"Intercept"=rep(intercept,length(DaysShows))))
    }
    
    Intrast<-ProcessDataRaster(ddint,value="Intercept",DoReplicate = T,Dorescale = T,rescalefrom = rescale_range)
    #p_Gene<-p_Gene+ geom_hline(yintercept=Intrast$Intercept,linetype="solid",color=alpha("black",.2),size=1)
    
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_circ_response$ZT,y=Model_outputs$out_response_circ_response$y1,
                              group=Model_outputs$out_response_circ_response$Day_axis,color=alpha(cols[["SCNforce"]],1),
                              size=ModelFitLine,linetype="solid")
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_sw_response$ZT,y=Model_outputs$out_response_sw_response$y1,
                              group=Model_outputs$out_response_sw_response$Day_axis,color=alpha(cols[["SWforce"]],1),
                              size=ModelFitLine,linetype="solid")
    if (ModelFitLine != 0){
      p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response$ZT,y=Model_outputs$out_response$y1,
                                group=Model_outputs$out_response$Day_axis,color="black",size=ModelFitLine,linetype="solid",alpha=1)
    }

    
  }
  
  # Add points and segments
  if (AddCISegment){
    TissueCols<-list("Cortex"="Cortex","Liver"="Liver","Desync"="Blood","Ctr"="Blood","Res"="Blood")
    Tissue<-TissueCols[[Dataset]]
    p_Gene <- p_Gene + annotate(x=CI_Expr_Raster$ZT,xend=CI_Expr_Raster$ZT,"segment",y=CI_Expr_Raster$lw,yend=CI_Expr_Raster$up,colour=cols[[Tissue]],size=1.2)
  }
  if (MeanPointSize > 0){
    p_Gene <- p_Gene + geom_point(shape=21,fill=cols[[Tissue]],col="black",size=MeanPointSize)
  }
  
  # hide some line
  p_Gene<-p_Gene + annotate("rect",xmin=24,xmax=Inf,ymin=-Inf,ymax=Inf,fill="white")
  
  # Add info for days
  p_Gene<-p_Gene + annotate("text",x=rep(-PreviousHourShow+.3,length(DaysShows)),y=(DaysShows-1)*(-1)-.1,label=paste("Day",DaysShows,sep=""),hjust=0,vjust=1,fontface="bold",size=3)
  
  # Add limits
  p_Gene<-p_Gene+coord_cartesian(xlim=c(-PreviousHourShow,xlim_max),ylim=c((range(DaysShows)[2])*-1,(range(DaysShows)[1]-1)*-1),clip = "off",expand = F)
  
  # X-lab
  p_Gene <- p_Gene + xlab("ZT [h]")
  
  
  # Add lower and upper equilibrium
  if (AddLowerUpperEquilibrium == T){
    LU<-AddUpperLowerEquilibrium(model = model,Dataset = Dataset,data = data,
                             ProbeGene =  ProbeGene)
    
    dd<-cbind.data.frame(Time=Model_outputs$out_response$Time,y1=rep(LU$Lower,length(Model_outputs$out_response$Time)))
    ddp<-ProcessDataRaster(dd,value=c("y1"),DoReplicate=T,Dorescale = T,rescalefrom = rescale_range,replastH =PreviousHourShow)
    ddp<-ddp[ddp$Day_axis %in% DaysShows,]
    p_Gene<-p_Gene + geom_hline(yintercept = unique(ddp$y1),col="red",linetype="22")
    
    dd<-cbind.data.frame(Time=Model_outputs$out_response$Time,y1=rep(LU$Upper,length(Model_outputs$out_response$Time)))
    ddp<-ProcessDataRaster(dd,value=c("y1"),DoReplicate=T,Dorescale = T,rescalefrom = rescale_range,replastH =PreviousHourShow)
    ddp<-ddp[ddp$Day_axis %in% DaysShows,]
    p_Gene<-p_Gene + geom_hline(yintercept = unique(ddp$y1),col="red",linetype="22")
    
  }
  
  return(p_Gene)
  
} 


AddUpperLowerEquilibrium<-function(model,Dataset,data,
                         ProbeGene){
  if (Dataset %in% c("Liver","Cortex")){
    WakeEffect<-data$fits[data$fits$Gene == ProbeGene,"Wake"]
    SleepEffect<-data$fits[data$fits$Gene == ProbeGene,"Sleep"]
    SpringConstant<-data$fits[data$fits$Gene == ProbeGene,"omega"]^2
    Upper<-model@intercept + WakeEffect*0.1/SpringConstant
    Lower<-model@intercept + SleepEffect*0.1/SpringConstant
  }else if (Dataset == "Desync"){
    WakeEffect<-data$fits[data$fits$ProbeName == ProbeGene,"Wake"]
    SleepEffect<-data$fits[data$fits$ProbeName == ProbeGene,"Sleep"]
    SpringConstant<-data$fits[data$fits$ProbeName == ProbeGene,"omega"]^2
    Upper<-model@intercept + WakeEffect*0.1/SpringConstant
    Lower<-model@intercept + SleepEffect*0.1/SpringConstant
  }else{
    WakeEffect<-data$fits[data$fits$ProbeName == ProbeGene,"Wake"]
    SleepEffect<-data$fits[data$fits$ProbeName == ProbeGene,"Sleep"]
    SpringConstant<-data$fits[data$fits$ProbeName == ProbeGene,"omega"]^2
    Upper<-data$fits[data$fits$ProbeName == ProbeGene,"Intercept"] + WakeEffect*0.1/SpringConstant
    Lower<-data$fits[data$fits$ProbeName == ProbeGene,"Intercept"] + SleepEffect*0.1/SpringConstant
  }

  message(paste("Lower Eq.:",round(Lower,2)," ; Upper Eq.:",round(Upper,2)))
  return(list(Lower=Lower,Upper=Upper))
}

AddLDblock<-function(p_Gene,Dataset,DaysShows,PreviousHourShow){
  rectgrey<-"grey80"
  if (Dataset %in% c("Liver","Cortex")){
    for (i in DaysShows){
      p_Gene<-p_Gene+annotate("rect",xmin=c(-Inf),xmax=c(0),ymin=c(i*-1),ymax=c(i*-1+1),fill=rectgrey)
      p_Gene<-p_Gene+annotate("rect",xmin=c(12),xmax=c(24),ymin=c(i*-1),ymax=c(i*-1+1),fill=rectgrey)
      if (3 %in% DaysShows){
        p_Gene<-p_Gene + annotate("rect",xmin=0,xmax=6,ymin=-3,ymax=-2,fill=rgb(1,0,0,.05)) # SD block
      }
    }
  }
  if (Dataset == "Desync"){
    p_Gene<-AddRectDesync(p_Gene,DaysShows,rectgrey,PreviousHourShow)
  }else if (Dataset == "Res"){
    p_Gene<-AddRectRest(p_Gene,DaysShows,rectgrey)
  }else if (Dataset == "Ctr"){
    p_Gene<-AddRectExt(p_Gene,DaysShows,rectgrey)
  }
  return(p_Gene)
}


GetModelOutput<-function(model,Dataset,data,ProbeGene,FitMethod,AddASkeldonSol,rescale_range,DaysShows,PreviousHourShow){
  
  if (Dataset %in% c("Liver","Cortex")){
    param_model<-data$fits[data$fits$Gene == ProbeGene,]
    intercept<-model@intercept
  }
  if (Dataset %in% c("Desync","Res","Ctr")){
    param_model<-data$fits[data$fits$ProbeName == ProbeGene,]
    if (Dataset == "Desync"){
      intercept<-model@intercept
    }else{
      intercept<-data$fits[data$fits$ProbeName == ProbeGene,"intercept"]
    }
  }
  out_response<-SWDMrFit(model,params = param_model,method=FitMethod)
  
  # Get unmodified version
  out_response_unrescaled<-out_response
  out_response_unrescaled$intercept<-intercept
  
  
  out_response_circ_response<-NULL
  out_response_sw_response<-NULL
  if (AddASkeldonSol ==T){
    
    
    out_response_circ_response<-cbind.data.frame(Time=out_response$time[-1],y1=out_response$circ_sol+intercept)
    out_response_circ_response<-ProcessDataRaster(out_response_circ_response,value="y1",DoReplicate = T,Dorescale = T,
                                                  rescalefrom = rescale_range,replastH = PreviousHourShow)
    out_response_circ_response<-out_response_circ_response[out_response_circ_response$Day_axis %in% DaysShows,]
    
    out_response_sw_response<-cbind.data.frame(Time=out_response$time[-1],y1=out_response$SW_response+out_response$trans_sol+intercept)
    out_response_sw_response<-ProcessDataRaster(out_response_sw_response,value="y1",DoReplicate = T,Dorescale = T,
                                                rescalefrom = rescale_range,replastH = PreviousHourShow)
    out_response_sw_response<-out_response_sw_response[out_response_sw_response$Day_axis %in% DaysShows,]
  }
  
  out_response<-cbind.data.frame(Time=out_response$time,y1=out_response$y1)
  
  out_response_bsl<-cbind.data.frame(Time=cumsum(rep(.1,2640)),
                                     y1=rep(out_response$y1[out_response$Time > 0.05 & out_response$Time <=24.0],11))
  
  out_response<-ProcessDataRaster(out_response,value="y1",DoReplicate = T,Dorescale = T,
                                  rescalefrom = rescale_range,replastH = PreviousHourShow)
  out_response<-out_response[out_response$Day_axis %in% DaysShows,]
  
  out_response_bsl<-ProcessDataRaster(out_response_bsl,value="y1",DoReplicate = T,Dorescale = T,
                                      rescalefrom = rescale_range,replastH = PreviousHourShow)
  out_response_bsl<-out_response_bsl[out_response_bsl$Day_axis %in% DaysShows & ! out_response_bsl$Day_axis %in% c(1,2),]
  
  return(list(out_response=out_response,out_response_bsl=out_response_bsl,
              out_response_circ_response=out_response_circ_response,
              out_response_sw_response=out_response_sw_response,
              out_response_unrescaled=out_response_unrescaled))
  
}


GetCIRaster<-function(Dataset,data,ProbeGene,DaysShows,ylim,PreviousHourShow,ReplicatesSamples){
  if (Dataset %in% c("Liver","Cortex")){
    CI_expr<-GetCI(data,ProbeGene)
    idx<-c(2,3,4)
  }
  if (Dataset %in% c("Desync","Res","Ctr")){
    # Add mean and SE based on emmeans
    CI_expr<-GetEMMEANS(ProbeGene,data)
    CI_expr$lw<-CI_expr$Mean-1.96*CI_expr$SE
    CI_expr$up<-CI_expr$Mean+1.96*CI_expr$SE
    idx<-c(4,8,9)
  }
  
  # Limits of the plot
  if (is.null(ylim)){
    rescale_range<-c(min(CI_expr[ceiling(CI_expr$Time/24) %in% DaysShows,idx]),max(CI_expr[ceiling(CI_expr$Time/24) %in% DaysShows,idx]))
  }else{
    rescale_range<-ylim
  }
  
  # Rasterize data
  CI_expr_raster<-ProcessDataRaster(CI_expr,value=c("Mean","lw","up"),DoReplicate=ReplicatesSamples,Dorescale = T,rescalefrom = rescale_range,replastH =PreviousHourShow)
  CI_expr_raster<-CI_expr_raster[CI_expr_raster$Day_axis %in% DaysShows,]
  if (Dataset %in% c("Desync","Res","Ctr")){
    CI_expr_raster<-CI_expr_raster[CI_expr_raster$Exp == Dataset,]
  }
  
  return(list(CI_expr_raster=CI_expr_raster,rescale_range=rescale_range,CI_expr=CI_expr))
}


GetModel<-function(Dataset,data,ProbeGene){
  
  if (Dataset %in% c("Liver","Cortex")){
    # Create object and model
    swdmr<- SWDMr(SWdist=data$SWdf, Gexp=data$rna_expr)
    model<-C57BL6_TimeCourse_GetFullModel(swdmr,ProbeGene)
  }
  if (Dataset %in% c("Desync","Res","Ctr")){
    # Get equilibrium position Desync experiment
    Expr<-data$HT_Desync[data$HT_Desync$Time < 72,ProbeGene]
    Time<-data$HT_Desync[data$HT_Desync$Time < 72,"Time"]
    MeanGeneExprInBaseline<-coef(lm(Expr~sin(2*pi/24*Time) + cos(2*pi/24*Time) ))[["(Intercept)"]]
    
    # Intercept difference
    IntDiff<-MeanGeneExprInBaseline-data$fits[data$fits$ProbeName == ProbeGene,"intercept"]
    
    # Build models
    if (Dataset == "Desync"){
      model<-Human_TimeCourse_GetFullModel(data$MeanSWdf,data$HT_Desync,ProbeGene,intercept=MeanGeneExprInBaseline)
    }else if (Dataset == "Ctr"){
      model<-Human_TimeCourse_GetFullModel(data$MeanSWdfExt,data$HT_Ext,ProbeGene)
    } else if (Dataset == "Res"){
      model<-Human_TimeCourse_GetFullModel(data$MeanSWdfRest,data$HT_Res,ProbeGene)
    }
  }

  return(model)
}


LoadData<-function(Dataset){
  
  if (Dataset == "Cortex"){
    data<-LOADMOUSE("Cortex")
  }
  if (Dataset == "Liver"){
    data<-LOADMOUSE("Liver")
  }
  if (Dataset %in% c("Desync","Ctr","Res")){
    data<-LOADHUMAN()
  }
  return(data)
}


AddRectDesync<-function(ggDF2,DaysShows,rectgrey,PreviousHourShow){
  
  for (i in DaysShows){
    if (i <= 2){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=8,ymin=i*-1,ymax=(i-1)*-1,fill=rectgrey)
    }else if (i == 3){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=9.2,ymin=-3,ymax=-2,fill=rectgrey)
    }else if (i == 4){
      ggDF2<-ggDF2+annotate("rect",xmin=0+4,xmax=9.2+4,ymin=-4,ymax=-3,fill=rectgrey)
    }else if (i == 5){
      ggDF2<-ggDF2+annotate("rect",xmin=0+8,xmax=9.2+8,ymin=-5,ymax=-4,fill=rectgrey)
    }else if (i == 6){
      ggDF2<-ggDF2+annotate("rect",xmin=0+12,xmax=9.2+12,ymin=-6,ymax=-5,fill=rectgrey)
    }else if (i == 7){
      ggDF2<-ggDF2+annotate("rect",xmin=0+16,xmax=24,ymin=-7,ymax=-6,fill=rectgrey)
      ggDF2<-ggDF2+annotate("rect",xmin=max(-PreviousHourShow,-2.8-9.2),xmax=-2.8,ymin=-7,ymax=-6,fill=rectgrey)
    }else if (i == 8){
      ggDF2<-ggDF2+annotate("rect",xmin=20,xmax=24,ymin=-8,ymax=-7,fill=rectgrey)
      ggDF2<-ggDF2+annotate("rect",xmin=max(-PreviousHourShow,1.2-9.2),xmax=1.2,ymin=-8,ymax=-7,fill=rectgrey)
    }else if (i == 9){
      ggDF2<-ggDF2+annotate("rect",xmin=max(-PreviousHourShow,5.2-9.2),xmax=5.2,ymin=-9,ymax=-8,fill=rectgrey)
    }else if (i == 10){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=9.2,ymin=-10,ymax=-9,fill=rectgrey)
    }
  }

  return(ggDF2)
}

AddRectExt<-function(ggDF2,DaysShows,rectgrey){
  for (i in DaysShows){
    print(i)
    if (i <= 2){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=8,ymin=i*-1,ymax=Inf,fill=rectgrey)
    }else if (i <= 9){
      ggDF2<-ggDF2+annotate("rect",xmin=-1,xmax=9,ymin=i*-1,ymax=(i-1)*-1,fill=rectgrey)
    }else if (i >= 11){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=12,ymin=i*-1,ymax=(i-1)*-1,fill=rectgrey)
    }
  }
  return(ggDF2)
}

AddRectRest<-function(ggDF2,DaysShows,rectgrey){
  
  for (i in DaysShows){
    print(i)
    if (i <= 2){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=8,ymin=i*-1,ymax=Inf,fill=rectgrey)
    }else if (i <= 9){
      ggDF2<-ggDF2+annotate("rect",xmin=1,xmax=7,ymin=i*-1,ymax=(i-1)*-1,fill=rectgrey)
    }else if (i >= 11){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=12,ymin=i*-1,ymax=(i-1)*-1,fill=rectgrey)
    }
  }
  return(ggDF2)
}


GetCI<-function(data,Gene){
  TF<-as.factor(data$rna_expr$Time)
  model_lm<-lm(data$rna_expr[,Gene]~0+TF)
  MeanCI<-cbind.data.frame(Time=as.numeric(gsub("TF","",names(coef(model_lm)))),Mean=coef(model_lm),lw=confint(model_lm)[,1],up=confint(model_lm)[,2])
  return(MeanCI)
}


# Functions for raster plot
# Data must contain 2 column, Time and other
ProcessDataRaster<-function(data,value="Sleep",rescalefrom=c(0,1),noreplast=F,Dorescale=T,DoReplicate=T,replastH=6){

  # Number of row to draw (n day/2)
  data$Day<-paste("Day",ceiling(data$Time/24))
  Num2Day<-max(floor(data$Time/24))
  
  # Add ZT time
  data$ZT<-(data$Time %% 24)
  data$ZT[data$ZT == 0]<-24
  
  data$Day_axis<-ceiling(data$Time/24)
  
  # Repeat
  dataRep<- RepeatData(data,noreplast,replastH)
  if (DoReplicate == F){
    dataRep<-dataRep[-which(dataRep$ZT >= 24 | dataRep$ZT <0 ),]
  }
  if (Dorescale ==T){
    dataRep[,value]<-apply(dataRep[,value,drop=F],2,function(x){rescale(x,rescalefrom,to=c(0.05,0.95))})
  }
  dataRep[,value]<-apply(dataRep[,value,drop=F],2,function(x){x-dataRep$Day_axis})
  return(dataRep)
}

# Repeat data for double plot
RepeatData<-function(data,noreplast,replastH=6){
  datarep<-data
  if (noreplast ==T){
    newline<-which(data$ZT>=(24-replastH) & data$Day_axis != max(data$Day_axis))
  }else{
    newline<-which(data$ZT>=(24-replastH))
  }
  for (i in newline){
    ToRep<-data[i,]
    ToRep$Day_axis<-ToRep$Day_axis+1
    ToRep$ZT<-ToRep$ZT-24
    datarep<-rbind(datarep,ToRep)
  }
  datarep<-datarep[order(datarep$Day_axis,datarep$Time),]
  return(datarep)
}



############ FUNCTION TO LOAD MOUSE DATA
LOADMOUSE<-function(Tissue="Cortex"){
  
  require(SWDMr)
  
  source(MOUSEMODELS_FUN) # Building functions
  source(COLORFUN) # Color code
  
  if (Tissue == "Liver"){
    load(LIVER_SWDMRDATA_RDATA)
    fits<-read.table(LIVERFITS,header=T)
  }else{
    load(CORTEX_SWDMRDATA_RDATA)
    fits<-read.table(CORTEXFITS,header=T)
  }
  
  return(list(rna_expr=rna_expr,SWdf=SWdf,fits=fits))
}


LOADHUMAN<-function(){
  
  require(SWDMr)
  
  source(HUMANMODELS_FUN) # Building functions
  source(COLORFUN) # Color code
  
  load(HUMANDATASET_RDATA)
  fits<-read.table(HUMANFITS,header=T)
  return(list(HT_Desync=HT_Desync,HT_Ext=HT_Ext,HT_Res=HT_Res,MeanSWdf=MeanSWdf,MeanSWdfExt=MeanSWdfExt,MeanSWdfRest=MeanSWdfRest,fits=fits))
}


##### GET CI human data
GetEMMEANS<-function(ProbeID,data){
  require(lme4)
  load(HUMANTRANSCRIPTOMEEMEANS_RDATA)
  emmeanres<-EmmeansMixModel(ExprData[rownames(ExprData) == ProbeID,],Experiment=Experiment2,TimeCond=TimeCond,Subjs=Subjs)
  
  df1<-cbind.data.frame(GSM=c(rownames(data$HT_Desync),rownames(data$HT_Ext),rownames(data$HT_Res)),
                        Time=c(data$HT_Desync$Time,data$HT_Ext$Time,data$HT_Res$Time),
                        Exp=c(rep("Desync",nrow(data$HT_Desync)),rep("Ctr",nrow(data$HT_Ext)) ,rep("Res",nrow(data$HT_Res)) ))
  rownames(df1)<-df1$GSM
  df2<-cbind.data.frame(GSM=colnames(ExprData),TimeCond=TimeCond)
  df1<-df1[colnames(ExprData),]
  
  dfall<-cbind.data.frame(df1,df2)
  dfall<-unique(dfall[,c("Time","Exp","TimeCond")])
  rownames(dfall)<-dfall$TimeCond
  
  dfall$Mean<-emmeanres$means[rownames(dfall)]
  dfall$SE<-emmeanres$se[rownames(dfall)]
  dfall$lwCI<-emmeanres$lwCI[rownames(dfall)]
  dfall$upCI<-emmeanres$upCI[rownames(dfall)]
  
  return(dfall)  
}

############ FUNCTION USE TO SQUISH X AXIS IN GIVEN INTERVAL
squish_day<-function(Rdays,idx_Rdays){
  require(scales)
  trans<-function(x){
    
    if (any(is.na(x))) return(x)
    
    sf<-rep(1,length(Rdays));sf[idx_Rdays]<-0.01
    transd<-cumsum(-1*rep(1,length(Rdays))*sf)-(abs(Rdays[1])-1)
    apfun<-approxfun(x=Rdays,y=transd,method="linear")
    
    # Max give value
    idxmax<-x>max(Rdays)
    idxmin<-x<min(Rdays)
    idxranged<-x<=max(Rdays) & x>=min(Rdays)
    x[idxmax]<-x[idxmax]
    x[idxmin]<-x[idxmin]-(min(Rdays)-min(transd))
    x[idxranged]<-apfun(x[idxranged])
    
    return(x)
  }
  
  inv<-function(x){

    if (any(is.na(x))) return(x)
    
    sf<-rep(1,length(Rdays));sf[idx_Rdays]<-0.01
    transd<-cumsum(-1*rep(1,length(Rdays))*sf)-(abs(Rdays[1])-1)
    apfun<-approxfun(y=Rdays,x=transd,method="linear")
    
    # Max give value
    idxmax<-x>max(transd)
    idxmin<-x<min(transd)
    idxranged<-x<=max(transd) & x>=min(transd)
    x[idxmax]<-x[idxmax]
    x[idxmin]<-x[idxmin]+((min(Rdays))-(min(transd)))
    x[idxranged]<-apfun(x[idxranged])
    
    return(x)
  }
  
  return(trans_new("squished", trans, inv))
}



########### SOME FIT TO ADD TO THE PLOT

MeanBSLfit<-function(dataexpr,Gene,bsl){
  Expr<-dataexpr[dataexpr$Time <= bsl,Gene]
  Time<-dataexpr[dataexpr$Time<= bsl,"Time"]
  ZT<-as.factor(Time)
  model<-lm(Expr~ZT)
  return(list(time=unique(Time),y1=predict(model,newdata=cbind.data.frame(ZT=unique(ZT)))))
}


MeanBSLfit_Rep<-function(dataexpr,Gene,bsl,EvalTime=NULL){
  Expr<-dataexpr[dataexpr$Time <= bsl,Gene]
  Time<-dataexpr[dataexpr$Time<= bsl,"Time"] %% 24
  ZT<-as.factor(Time)
  model<-lm(Expr~ZT)
  if (is.null(EvalTime)){
    Time_dataset<-sort(unique(dataexpr$Time))
  }else{
    Time_dataset<-EvalTime
  }
  ZTime_dataset<-as.factor(Time_dataset %% 24)
  
  return(list(time=Time_dataset,y1=predict(model,newdata=cbind.data.frame(ZT=ZTime_dataset))))
}


BSLRepfit<-function(dataexpr,Gene,bsl,EvalTime=NULL){
  Expr<-data$rna_expr[data$rna_expr$Time<=tlimits,Gene]
  Time<-data$rna_expr[data$rna_expr$Time<=tlimits,"Time"] %% 24
  ZT<-as.factor(Time)
  model<-lm(Expr~ZT)
  if (is.null(EvalTime)){
    Time<-data$SWdf$Time[data$SWdf$Time>=48]
    Time<-Time[(Time %% 24) %in% c(0,3,6,12,18)]
  }else{
    
  }

  ZT<-as.factor(Time %% 24)
  fitted<-predict(model,newdata = ZT)
  return(list(time=Time,y1=fitted))
}

RecRepfit<-function(Gene,data){
  Expr<-data$rna_expr[data$rna_expr$Time>=tlimits,Gene]
  Time<-data$rna_expr[data$rna_expr$Time>=tlimits,"Time"]
  ZT<-as.factor(Time)
  model<-lm(Expr~ZT)
  return(list(time=Time,y1=fitted(model)))
}

########### LINEAR PLOT


AddLDblocks<-function(gg,starts,lengthblock=12){
  for (i in starts){
    gg<-gg + annotate("rect", xmin = i, xmax = i+lengthblock, ymin = -Inf, ymax = Inf, fill="black",alpha = 0.2)
  }
  return(gg)
}


PlotForces<-function(Gene,Tissue,xlim=c(0,228)){
  
  # Color code 
  cols<-ColorCode()
  
  require(ggplot2)
  
  # LOAD DATA AND FITS
  data<-LOADMOUSE(Tissue=Tissue)
  
  fits<-data$fits[data$fits$Gene == Gene,]
  
  dfd<-as.data.frame(cbind(Time=data$SWdf$Time,SW=data$SWdf$Sleep*fits$Sleep+data$SWdf$Wake*fits$Wake,
                           Circ=fits$AmpSin*sin(2*pi/24*data$SWdf$Time+fits$PhiSin)))
  
  # Initiate plot
  p_Forces <- ggplot(dfd,aes(x=Time,y=SW))
  
  # Add LD blocks
  p_Forces<-AddLDblocks(p_Forces,seq(12,204,by=24))
  p_Forces <- p_Forces + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = Inf, fill="red",alpha = 0.3)
  
  p_Forces<-p_Forces+geom_line(color=cols[["SWforce"]])
  p_Forces<-p_Forces + annotate("line",x=dfd$Time,y=dfd$Circ,color=cols[["SCNforce"]])
  
  
  # Squish time 108 to 203
  p_Forces <- p_Forces + scale_x_continuous(trans = squish_trans(108, 203, 30),
                                            breaks = c(seq(0, 96, by = 12),204,216,228),limits=xlim)
  # Set transparent background
  p_Forces <- p_Forces + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA)) 
  
  # Hide squish
  p_Forces<-p_Forces + annotate("rect",xmin=105,xmax=204, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
  
  p_Forces<-p_Forces+ylab("Forces")
  
  p_Forces
  
}


########### PLOT GENE FITS FOR MOUSE
PlotMouseFit<-function(Gene,Tissue,PointSize=1,xlim=c(0,228),
                       lineSize=1,MeanPointSize=3,bsldashed=F,
                       lineSizeBSL=1,AddFcline=F,AddFsline=F,
                       goldbsl=F,BootstrapResults=NULL,ASContribution=F,AddSleepWake=F,...){
  require(ggplot2)
  require(SWDMr)
  require(dplyr)
  
  # LOAD DATA AND FITS
  data<-LOADMOUSE(Tissue=Tissue)
  
  # Color code 
  cols<-ColorCode()
  
  # Create object and model
  swdmr<- SWDMr(SWdist=data$SWdf, Gexp=data$rna_expr)
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  
  # STANDARD PLOT
  if (ASContribution == F){
    # model simulation
    FittedData<-SWDMrFit(object = model,params = data$fits[data$fits$Gene == Gene,])
    
    # Gene expression
    dfd<-as.data.frame(cbind(model@Gexp$Time,model@Gexp[,model@VarExp]))
    colnames(dfd)<-c("Time","Gene")
    
    # Initiate plot
    p_Gene <- ggplot(dfd,aes(x=Time,y=Gene))
    
    # Add LD blocks
    p_Gene<-AddLDblocks(p_Gene,seq(12,204,by=24))
    
    # Given other fits ?
    if (AddSleepWake == T){
      SW<-model@SWdist
      SW<-SW %>% group_by(ceiling(Time)) %>% mutate_at("Wake",mean)
      # Rescale SW
      minv<-min(dfd$Gene)
      maxv<-max(dfd$Gene)
      diff<-maxv-minv
      SW$Wake<-SW$Wake * (minv+0.3*diff-minv) + minv
      dd<-cbind.data.frame(Time=SW$Time,y1=SW$Wake)
      idxrib<-dd$Time>=min(xlim) & dd$Time<=max(xlim)
      p_Gene<-p_Gene +
        annotate("ribbon",x=dd$Time[idxrib],ymax=dd$y1[idxrib],ymin=-Inf,#-Inf,minv
                 fill=alpha(cols[["SWforce"]],.2))+
        annotate("line",x=dd$Time[idxrib],y=dd$y1[idxrib],color=alpha(cols[["SWforce"]],.7))
    }
    
    
    # Compute Mean, SD and SE for each time points
    #meansd<-SWDMr:::StatsPerTimePoint(model)
    TF<-as.factor(data$rna_expr$Time)
    model<-lm(data$rna_expr[,Gene]~0+TF)
    MeanCI<-cbind.data.frame(Time=as.numeric(gsub("TF","",names(coef(model)))),Mean=coef(model),lw=confint(model)[,1],up=confint(model)[,2])
    # Add error bar and mean 
    #p_Gene <- p_Gene + annotate(x=MeanCI$Time,"errorbar",ymin=MeanCI$lw,ymax=MeanCI$up,colour=cols[[Tissue]], width=1,size=1)
    p_Gene <- p_Gene + annotate(x=MeanCI$Time,xend=MeanCI$Time,"segment",y=MeanCI$lw,yend=MeanCI$up,colour=cols[[Tissue]],size=1)
    p_Gene <- p_Gene + annotate("point",x=MeanCI$Time,y=MeanCI$Mean,size=MeanPointSize,shape=21,fill=cols[[Tissue]],col=cols[[Tissue]]) #
    
    # p_Gene <- p_Gene + annotate(x=unique(model@Gexp$Time),"errorbar",ymin=meansd[,"mean"]-meansd[,"se"],ymax=meansd[,"mean"]+meansd[,"se"],colour=cols[[Tissue]], width=1,size=1)
    # p_Gene <- p_Gene + annotate("point",x=unique(model@Gexp$Time),y=meansd[,"mean"],size=MeanPointSize,shape=21,fill=cols[[Tissue]]) #,col="black"
    
    # add points
    if (PointSize > 0){
      p_Gene <- p_Gene + geom_point(color = cols[[Tissue]],size=PointSize)
    }
    
    
    # Squish time 108 to 203
    p_Gene <- p_Gene + scale_x_continuous(trans = squish_trans(108, 203, 30),
                                          breaks = c(seq(0, 96, by = 12),204,216,228),limits=xlim)
    
    # Set transparent background
    p_Gene <- p_Gene + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA)) 
    
    
    if (bsldashed ==T){
      bslF<-rep(FittedData$y1[FittedData$time > 0.0 & FittedData$time <=24.0],8)
      p_Gene <- p_Gene + annotate("line",x=seq(48.1,240,by=.1),y=bslF,color=rgb(0,0,0,.25),size=lineSizeBSL,linetype = "dashed") 
    }
    
    # SD rect
    rangev<-ggplot_build(p_Gene)$layout$panel_params[[1]]$y.range
    #p_Gene <- p_Gene + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = rangev[[1]]+(rangev[[2]]-rangev[[1]])*0.1, fill="red",alpha = 0.3)
    p_Gene <- p_Gene + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = Inf, fill="red",alpha = 0.3)
    
    # Fitted line
    if (lineSize>0){
      p_Gene <- p_Gene + annotate("line",x=FittedData$time[FittedData$time>=0],y=FittedData$y1[FittedData$time>=0],color=cols[["Modelfit"]],size=lineSize,alpha=1)
    }
    if (goldbsl == T){
      p_Gene <- p_Gene + annotate("line",x=FittedData$time[FittedData$time>=0 & FittedData$time<=24],y=FittedData$y1[FittedData$time>=0 & FittedData$time<=24],
                                  color=cols[["SCNforce"]],size=lineSize,alpha=1)
    }
    
    # Fitted line Fcircadian only
    if (AddFcline == T){
      if (Tissue == "Liver"){
        fit_Fc<-read.table(FILE_LIVERFITS_DropSW,header=T)
      }else{
        fit_Fc<-read.table(FILE_CORTEXFITS_DropSW,header=T)
      }
      fit_Fc$dampratio<-fit_Fc$dampratio/(2*fit_Fc$omega)^2
      modelFc<-C57BL6_TimeCourse_GetCircadianModel(swdmr,Gene)
      FittedDataFc<-SWDMrFit(object = modelFc,params = fit_Fc[fit_Fc$Gene == Gene,])
      p_Gene <- p_Gene + annotate("line",x=FittedDataFc$time[FittedDataFc$time>=0],
                                  y=FittedDataFc$y1[FittedDataFc$time>=0],color=cols[["SCNforce"]],size=lineSize,alpha=.8,linetype = "solid")
    }
    
    # Fitted line Fsleep-wake only
    if (AddFsline == T){
      if (Tissue == "Liver"){
        fit_Fs<-read.table(FILE_LIVERFITS_DropSin,header=T)
      }else{
        fit_Fs<-read.table(FILE_CORTEXFITS_DropSin,header=T)
      } 
      fit_Fs$dampratio<-fit_Fs$dampratio/(2*fit_Fs$omega)^2
      modelFs<-C57BL6_TimeCourse_GetNoSinWaveModel(swdmr,Gene)
      FittedDataFs<-SWDMrFit(object = modelFs,params = fit_Fs[fit_Fs$Gene == Gene,])
      p_Gene <- p_Gene + annotate("line",x=FittedDataFs$time[FittedDataFs$time>=0],
                                  y=FittedDataFs$y1[FittedDataFs$time>=0],color=cols[["SWforce"]],size=lineSize,alpha=.8,linetype = "solid")
      
    }
    
    
    if (! is.null(BootstrapResults)){
      p_Gene <-p_Gene + annotate("ribbon",x=BootstrapResults$df.mc$Time, ymin=BootstrapResults$df.mc$lwr.conf, ymax=BootstrapResults$df.mc$upr.conf, alpha=0.4, fill=cols[[Tissue]])
    }
    
    # if (ASContribution == T){
    #   model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
    #   out_response<-SWDMrFit(model,data$fits[data$fits$Gene == Gene,],method="Solve")
    #   p_Gene<-p_Gene + annotate("line",x=out_response$time[-1],y=out_response$circ_sol,col=cols[["SCNforce"]])
    #   p_Gene<-p_Gene + annotate("line",x=out_response$time[-1],y=out_response$trans_sol + out_response$SW_response,col=cols[["SWforce"]])
    # }
    
    # Hide squish
    p_Gene<-p_Gene + annotate("rect",xmin=105,xmax=204, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
    
    
    
    return(p_Gene)
    
    # PLOT using contribution of sleep-wake and circadian 
  }else{
    coef_sleep<-data$fits[data$fits$Gene == Gene,"Sleep"]
    coef_wake<-data$fits[data$fits$Gene == Gene,"Wake"]
    phase<-data$fits[data$fits$Gene == Gene,"PhiSin"]
    amp<-data$fits[data$fits$Gene == Gene,"AmpSin"]
    
    # out_response<-SWDMrFit(model,data$fits[data$fits$Gene == Gene,],method="Solve")
    # totalreponse<-out_response$circ_sol + out_response$trans_sol + out_response$SW_response
    # data_response=cbind.data.frame(time=out_response$time[-1],
    #                                response=c(out_response$circ_sol,amp * sin(2*pi/24*out_response$time[-1]+phase),
    #                                           out_response$trans_sol + out_response$SW_response,
    #                                           model@SWdist$Sleep*coef_sleep + model@SWdist$Wake*coef_wake,
    #                                           totalreponse),
    #                                value=c(rep("Circadian response",length(out_response$time[-1])),
    #                                        rep("Circadian force",length(out_response$time[-1])),
    #                                        rep("Sleep-Wake response",length(out_response$time[-1])),
    #                                        rep("Sleep-Wake force",length(out_response$time[-1])),
    #                                        rep("Total response",length(out_response$time[-1]))))
    #
    # gg_resp<-ggplot(aes(time,y=response,color=value),data=data_response)+geom_line(size=1)
    # gg_resp<-gg_resp+scale_x_continuous(breaks=seq(-24,96,by=12),limits = c(-24,96))+ theme_bw() + ggtitle(paste(Gene,"- Responses"))
    # gg_resp <- gg_resp + scale_color_manual(values=c("Circadian response" = cols[["SCNforce"]],
    #                                                  "Sleep-Wake response" = cols[["SWforce"]],
    #                                                  "Total response" = cols[["Modelfit"]],
    #                                                  "Circadian force" =alpha(cols[["SCNforce"]],.5),
    #                                                  "Sleep-Wake force" = alpha(cols[["SWforce"]],.5)))
    
    out_response<-SWDMrFit(model,data$fits[data$fits$Gene == Gene,],method="Solve")
    totalreponse<-out_response$circ_sol + out_response$trans_sol + out_response$SW_response
    data_response=cbind.data.frame(time=out_response$time[-1],
                                   response=c(out_response$circ_sol+model@intercept,
                                              out_response$trans_sol + out_response$SW_response+model@intercept,
                                              totalreponse+model@intercept),
                                   value=c(rep("Circadian response",length(out_response$time[-1])),
                                           rep("Sleep-Wake response",length(out_response$time[-1])),
                                           rep("Total response",length(out_response$time[-1]))))
    
    
    
    gg_resp<-ggplot(aes(time,y=response,color=value,linetype=value),data=data_response)+geom_line(size=1)
    
    # Add LD blocks
    gg_resp<-AddLDblocks(gg_resp,seq(12,204,by=24))
    
    rangev<-ggplot_build(gg_resp)$layout$panel_params[[1]]$y.range
    #gg_resp <- gg_resp + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = rangev[[1]]+(rangev[[2]]-rangev[[1]])*0.1, fill="red",alpha = 0.3)
    gg_resp <- gg_resp + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = Inf, fill="red",alpha = 0.3)
    
    
    gg_resp<-gg_resp+scale_x_continuous(trans = squish_trans(108, 203, 30),
                                        breaks = c(seq(0, 96, by = 12),204,216,228),limits=xlim)+ theme_bw() + ggtitle(paste(Gene,":",Tissue," - Responses",sep=""))
    gg_resp <- gg_resp + scale_color_manual(values=c("Circadian response" = cols[["SCNforce"]],
                                                     "Sleep-Wake response" = cols[["SWforce"]],
                                                     "Total response" = alpha(cols[["Modelfit"]],1)))
    gg_resp<-gg_resp + scale_linetype_manual(values=c("Circadian response" = "solid",
                                                      "Sleep-Wake response" = "solid",
                                                      "Total response" = "solid"))
    
    gg_resp<-gg_resp + annotate("rect",xmin=105,xmax=204, ymin = -Inf, ymax = Inf, fill="white",alpha = 1)
    
    gg_resp<-gg_resp + theme_classic() + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA),
                                               legend.position = "none")
    gg_resp<-gg_resp+ylab("Responses")
    return(gg_resp)
    
  }
}
