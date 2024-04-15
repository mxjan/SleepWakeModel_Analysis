############# FUNCTION FOR PLOT FITTING IN HUMAN and MOUSE ###############
# FOR TEST
ProbeGene<-"Homer1"
Dataset<-"Cortex"
DaysShows =c(1,2,3,4,5,10)
npretty=3
ylim=NULL
MeanPointSize=1
AddCISegment=T
ModelFitLine=1
FitMethod="Solve"
AddASkeldonSol=T
BSL_RepLine=1
AddFits=NULL
ColsAddFits=NULL
LtyAddFits=NULL
AddLowerUpperEquilibrium=F
AddSleepWake=T
ScaleSW=T
ReplicatesSamples=F
PreviousHourShow=6


##################################################
#---------- LINEAR PLOT FUNCTION ----------------#

LinearPlot<-function(ProbeGene,
                     Dataset,# Cortex,Liver,Desync,Res,Ctr
                     DaysShows=seq(1,10),
                     npretty=3,
                     ylim=NULL,
                     MeanPointSize=1,
                     AddCISegment=T, # Expression points
                     ModelFitLine=1,
                     FitMethod="Solve",
                     AddASkeldonSol=F,
                     BSL_RepLine=0, # model fits
                     AddFits=NULL,
                     ColsAddFits=NULL,
                     LtyAddFits=NULL,
                     AddLowerUpperEquilibrium=F,
                     AddSleepWake=F,
                     ScaleSW=T,
                     colorFit="black",
                     AddIntercept=T,
                     PlotForces=F){
  
  require(ggplot2)
  require(SWDMr)
  require(scales)
  require(dplyr)
  
  # Load data given dataset
  data<-LoadData(Dataset)
  
  # Color Codes / See RFunctions/color.R
  source(COLORFUN)
  cols<-ColorCode()
  
  # Get Model
  model<-GetModel(Dataset,data,ProbeGene)
  
  # Get Model fitted line(s)
  Model_outputs<-GetModelOutput(model = model,Dataset = Dataset,data = data,
                                ProbeGene =  ProbeGene,FitMethod = FitMethod,
                                AddASkeldonSol = AddASkeldonSol,
                                DaysShows = DaysShows,rescale=F)
  
  
  # Get Confidence Interval and mean time point of expression
  if (Dataset %in% c("Liver","Cortex")){
    CI_expr<-GetCI(data,ProbeGene)
  }
  if (Dataset %in% c("Desync","Res","Ctr")){
    # Add mean and SE based on emmeans
    CI_expr<-GetEMMEANS(ProbeGene,data)
    CI_expr$lw<-CI_expr$Mean-1.96*CI_expr$SE
    CI_expr$up<-CI_expr$Mean+1.96*CI_expr$SE
    CI_expr<-CI_expr[CI_expr$Exp == Dataset,]
  }
  
  dfd<-as.data.frame(cbind(CI_expr$Time,CI_expr[,"Mean"],floor(CI_expr$Time / 24)+1))
  colnames(dfd)<-c("Time","Gene","Day")
  dfd<-dfd[dfd$Day %in% DaysShows,]
  
  ####### PLOT ########
  p_Gene <- ggplot(dfd,aes(x=Time,y=Gene))
  
  # Add sleep opportunity or LD block
  p_Gene<-AddLDblockLinear(p_Gene,Dataset,DaysShows)
  
  # Add ZT0 line
  p_Gene<-p_Gene+geom_vline(xintercept=(DaysShows-1)*24,color="black",size=0.5,linetype="21")
  
  # Theme
  p_Gene<-p_Gene+ theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  
  # Scale x-axis
  Rdays<-seq(range(DaysShows)[1],range(DaysShows)[2])
  idx_Rdays<-! Rdays %in% ((DaysShows))
  labelsToPlot<-sort(unique(c((DaysShows-1)*24+12,(DaysShows)*24,(DaysShows-1)*24)))
  if (any(idx_Rdays)){
    p_Gene <- p_Gene + scale_x_continuous(trans = squish_trans(Rdays,idx_Rdays),#,
                                          breaks=labelsToPlot,
                                          labels=labelsToPlot,
                                          limits=c(min(labelsToPlot),max(labelsToPlot))) # Add y-axis labels on all days
  }else{
    p_Gene <- p_Gene + scale_x_continuous(breaks=labelsToPlot,
                                          labels=labelsToPlot,
                                          limits=c(min(labelsToPlot),max(labelsToPlot)))
  }
  
  # scale y-axis
  pretty_fun<-scales::breaks_pretty(n=npretty)
  if (is.null(ylim)){
    Expr_labels<-pretty_fun(CI_expr$Mean)
  }else{
    Expr_labels<-pretty_fun(ylim)
  }
  p_Gene<-p_Gene + scale_y_continuous(breaks=Expr_labels)#,limits=c(min(Expr_labels),max(Expr_labels)))
  
  # Add some line for expression
  #p_Gene<-p_Gene + geom_hline(yintercept = Expr_labels, color=rgb(0,0,0,.1),linetype="solid",size=0.25)
  
  
  if (PlotForces == T){
    SW<-model@SWdist
    fits<-data$fits[data$fits$Gene == Gene,]
    dfd<-as.data.frame(cbind(Time=data$SWdf$Time,SW=data$SWdf$Sleep*fits$Sleep+data$SWdf$Wake*fits$Wake,
                             Circ=fits$AmpSin*sin(2*pi/24*data$SWdf$Time+fits$PhiSin)))
    p_Gene<-p_Gene + annotate("line",x=dfd$Time,y=dfd$Circ,color=cols[["SCNforce"]])
    p_Gene<-p_Gene + annotate("line",x=dfd$Time,y=dfd$SW,color=cols[["SWforce"]])
    
  }
  
  if (AddSleepWake == T){
    if (is.null(ylim)){
      rangeV<-c(min(CI_expr$lw),max(CI_expr$up))
    }else{
      rangeV<-ylim
    }
    
    SW<-model@SWdist
    SW[which.min(abs(SW$Time)),"Time"]<-0 # Be sure it is 0, and not 0.0000001 -> ceiling will make it 1 instead of 0
    SW<-SW %>% group_by(ceiling(Time)) %>% mutate_at("Wake",mean)
    SW<-SW %>% group_by(ceiling(Time)) %>% mutate_at("Sleep",mean)
    SW$Wake<-SW$Wake
    #SW$Wake<-SW$Wake/.1 * (0.6 - 0.01) + 0.01
    dd<-cbind.data.frame(Time=SW$Time,y1=SW$Wake,Day_axis=floor(SW$Time/24)+1)
    dd<-dd[dd$Day_axis %in% DaysShows,]
    
    if (ScaleSW == T){
      midp<-rangeV[1]+(rangeV[2]-rangeV[1])/3
      dd$y1<-rangeV[1]+(midp-rangeV[1])*dd$y1
    }
    # 
    # 
    for (i in DaysShows){
      p_Gene<-p_Gene +
        annotate("ribbon",x=dd$Time[dd$Day_axis==i],ymax=dd$y1[dd$Day_axis==i],ymin=-Inf,
                 fill=alpha(cols[["SWforce"]],.2))+
        annotate("line",x=dd$Time[dd$Day_axis==i],y=dd$y1[dd$Day_axis==i],color=alpha(cols[["SWforce"]],.7))
    }
    
  }
  
  # BSL replicated line
  if (BSL_RepLine != 0){
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_bsl$Time,y=Model_outputs$out_response_bsl$y1,
                              color=alpha("black",.65),size=BSL_RepLine,linetype="42")
  }
  
  # Fitted lines
  if (ModelFitLine != 0 & AddASkeldonSol ==F){
    p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response$Time,y=Model_outputs$out_response$y1,
                              color=colorFit,size=ModelFitLine,alpha=1)
  }
  
  # Given other fits ?
  if (! is.null(AddFits)){
    for (i in 1:length(AddFits)){
      p_Gene<-p_Gene + annotate("path",x=AddFits[[i]]$time,y=AddFits[[i]]$y1,
                                color=ColsAddFits[i],size=.75,linetype=LtyAddFits[[i]])
      
    }
  }
  
  # Add line for Anne Skeldon solution
  if (AddASkeldonSol ==T){
    if (Dataset %in% c("Cortex","Liver","Desync")){
      intercept<-model@intercept
      ddint<-as.data.frame(cbind("Time"=12+((DaysShows-1)*24),"Intercept"=rep(model@intercept,length(DaysShows))))
    }else{
      intercept<-data$fits[data$fits$ProbeName == ProbeGene,"intercept"]
      ddint<-as.data.frame(cbind("Time"=12+((DaysShows-1)*24),"Intercept"=rep(intercept,length(DaysShows))))
    }
    if (AddIntercept == T){
      p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_circ_response$Time,y=Model_outputs$out_response_circ_response$y1
                                ,color=alpha(cols[["SCNforce"]],1),
                                size=ModelFitLine,linetype="solid")
      p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_sw_response$Time,y=Model_outputs$out_response_sw_response$y1,
                                color=alpha(cols[["SWforce"]],1),
                                size=ModelFitLine,linetype="solid")
      if (ModelFitLine != 0){
        p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response$Time,y=Model_outputs$out_response$y1,
                                  color="black",size=ModelFitLine,linetype="solid",alpha=1)
      }
    }else{
      p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_circ_response$Time,y=Model_outputs$out_response_circ_response$y1-intercept
                                ,color=alpha(cols[["SCNforce"]],1),
                                size=ModelFitLine,linetype="solid")
      p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response_sw_response$Time,y=Model_outputs$out_response_sw_response$y1-intercept,
                                color=alpha(cols[["SWforce"]],1),
                                size=ModelFitLine,linetype="solid")
      if (ModelFitLine != 0){
        p_Gene<-p_Gene + annotate("path",x=Model_outputs$out_response$Time,y=Model_outputs$out_response$y1-intercept,
                                  color="black",size=ModelFitLine,linetype="solid",alpha=1)
      }
    }

    
    
  }
  
  
  # Add white rect hide 
  if (length(idx_Rdays)>0){
    p_Gene<-p_Gene + annotate("rect",xmin=(Rdays[idx_Rdays]-1)*24,xmax=(Rdays[idx_Rdays])*24,
                              ymin=-Inf,ymax=Inf,fill="white")
  }

  
  
  # Add points and segments
  if (AddCISegment){
    TissueCols<-list("Cortex"="Cortex","Liver"="Liver","Desync"="Blood","Ctr"="Blood","Res"="Blood")
    Tissue<-TissueCols[[Dataset]]
    p_Gene <- p_Gene + annotate(x=CI_expr$Time,
                                xend=CI_expr$Time,"segment",
                                y=CI_expr$lw,yend=CI_expr$up,
                                colour=cols[[Tissue]],size=1.2,
                                lineend="round",linejoin="round")
  }
  if (MeanPointSize > 0){
    p_Gene <- p_Gene + geom_point(shape=21,fill=cols[[Tissue]],col="black",size=MeanPointSize)
  }
  
  # hide some line
  #p_Gene<-p_Gene + annotate("rect",xmin=24,xmax=Inf,ymin=-Inf,ymax=Inf,fill="white")
  
  # # Add limits
  if (is.null(ylim)){
    p_Gene<-p_Gene+coord_cartesian(ylim=c(min(CI_expr$lw),max(CI_expr$up)),clip = "off",expand = F)
  }else{
    p_Gene<-p_Gene+coord_cartesian(ylim=ylim,clip = "off",expand = F)
  }
  
  
  # X-lab
  p_Gene <- p_Gene + xlab("Time [h]")
  
  return(p_Gene)

  #  p_Gene<-p_Gene+scale_x_continuous(breaks=seq(min((DaysShows-1)*24),max((DaysShows-1)*24),by=12))
}





#---------- LINEAR PLOT FUNCTION ----------------#
##################################################


##################################################
#------------------ FUNCTIONS -------------------#


GetRangeExpressionHS<-function(ProbeGene){
  
  
  # Load data given dataset
  data<-LoadData("Res")
  
  # Color Codes / See RFunctions/color.R
  source(COLORFUN)
  cols<-ColorCode()
  
  # Get Model
  model0<-GetModel("Desync",data,ProbeGene)
  diffIntercept<- data$fits$intercept[data$fits$ProbeName == ProbeGene]-model0@intercept
  
  CI_expr<-GetEMMEANS(ProbeGene,data)
  CI_expr$lw<-CI_expr$Mean-1.96*CI_expr$SE
  CI_expr$up<-CI_expr$Mean+1.96*CI_expr$SE
  
  CI_expr$lwCI[CI_expr$Exp %in% c("Ctr","Res")]<-CI_expr$lwCI[CI_expr$Exp %in% c("Ctr","Res")]-diffIntercept
  CI_expr$upCI[CI_expr$Exp %in% c("Ctr","Res")]<-CI_expr$upCI[CI_expr$Exp %in% c("Ctr","Res")]-diffIntercept
  
  Range<-max(CI_expr$upCI)-min(CI_expr$lwCI)
  
  return(list(ylim1=c(min(CI_expr$lwCI),max(CI_expr$upCI)),
         ylim2=c(min(CI_expr$lwCI),max(CI_expr$upCI))+diffIntercept))
  
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


AddLDblockLinear<-function(p_Gene,Dataset,DaysShows){
  rectgrey<-"grey80"
  if (Dataset %in% c("Liver","Cortex")){
    for (i in DaysShows){
      p_Gene<-p_Gene+annotate("rect",xmin=c((DaysShows-1) * 24 + 12),xmax=c((DaysShows) * 24),ymin=-Inf,ymax=Inf,fill=rectgrey)
      if (3 %in% DaysShows){
        p_Gene<-p_Gene + annotate("rect",xmin=48,xmax=54,ymin=-Inf,ymax=Inf,fill=rgb(1,0,0,.05)) # SD block
      }
    }
  }
  if (Dataset == "Desync"){
    p_Gene<-AddRectDesyncLinear(p_Gene,DaysShows,rectgrey)
  }else if (Dataset == "Res"){
    p_Gene<-AddRectRestLinear(p_Gene,DaysShows,rectgrey)
  }else if (Dataset == "Ctr"){
    p_Gene<-AddRectExtLinear(p_Gene,DaysShows,rectgrey)
  }
  return(p_Gene)
}


GetModelOutput<-function(model,Dataset,data,ProbeGene,FitMethod,AddASkeldonSol,rescale_range,DaysShows,PreviousHourShow,rescale=T){
  
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
    message("Intercept:",intercept)
  }
  out_response<-SWDMrFit(model,params = param_model,method=FitMethod)
  
  # Get unmodified version
  out_response_unrescaled<-out_response
  out_response_unrescaled$intercept<-intercept
  
  
  out_response_circ_response<-NULL
  out_response_sw_response<-NULL
  if (AddASkeldonSol ==T){
    
    
    out_response_circ_response<-cbind.data.frame(Time=out_response$time[-1],y1=out_response$circ_sol+intercept)
    if (rescale == T){
      out_response_circ_response<-ProcessDataRaster(out_response_circ_response,value="y1",DoReplicate = T,Dorescale = T,
                                                    rescalefrom = rescale_range,replastH = PreviousHourShow)
      out_response_circ_response<-out_response_circ_response[out_response_circ_response$Day_axis %in% DaysShows,]
    }
    
    
    out_response_sw_response<-cbind.data.frame(Time=out_response$time[-1],y1=out_response$SW_response+out_response$trans_sol+intercept)
    if (rescale == T){
      out_response_sw_response<-ProcessDataRaster(out_response_sw_response,value="y1",DoReplicate = T,Dorescale = T,
                                                rescalefrom = rescale_range,replastH = PreviousHourShow)
      out_response_sw_response<-out_response_sw_response[out_response_sw_response$Day_axis %in% DaysShows,]
    }
    
  }
  
  out_response<-cbind.data.frame(Time=out_response$time,y1=out_response$y1)
  
  out_response_bsl<-cbind.data.frame(Time=cumsum(rep(.1,2640)),
                                     y1=rep(out_response$y1[out_response$Time > 0.05 & out_response$Time <=24.0],11))
  if (rescale == T){
    out_response<-ProcessDataRaster(out_response,value="y1",DoReplicate = T,Dorescale = T,
                                  rescalefrom = rescale_range,replastH = PreviousHourShow)
    out_response<-out_response[out_response$Day_axis %in% DaysShows,]
  }
  
  if (rescale == T){
    out_response_bsl<-ProcessDataRaster(out_response_bsl,value="y1",DoReplicate = T,Dorescale = T,
                                      rescalefrom = rescale_range,replastH = PreviousHourShow)
    out_response_bsl<-out_response_bsl[out_response_bsl$Day_axis %in% DaysShows & ! out_response_bsl$Day_axis %in% c(1,2),]
  }
  
  
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


AddRectDesyncLinear<-function(ggDF2,DaysShows,rectgrey){
  for (i in DaysShows){
    if (i == 1){
      ggDF2<-ggDF2+annotate("rect",xmin=0,xmax=8,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }
    else if (i == 2){
      ggDF2<-ggDF2+annotate("rect",xmin=24,xmax=32,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 3){
      ggDF2<-ggDF2+annotate("rect",xmin=48,xmax=(48+9.2),ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 4){
      ggDF2<-ggDF2+annotate("rect",xmin=72+4,xmax=72+9.2+4,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 5){
      ggDF2<-ggDF2+annotate("rect",xmin=96+8,xmax=96+9.2+8,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 6){
      ggDF2<-ggDF2+annotate("rect",xmin=120+12,xmax=120+9.2+12,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 7){
      ggDF2<-ggDF2+annotate("rect",xmin=144+16,xmax=144+16+9.2,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 8){
      ggDF2<-ggDF2+annotate("rect",xmin=168+20,xmax=168+20+9.2,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i == 9){
      ggDF2<-ggDF2+annotate("rect",xmin=192+24,xmax=192+24+9.2,ymin=-Inf,ymax=Inf,fill=rectgrey)
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


AddRectExtLinear<-function(ggDF2,DaysShows,rectgrey){
  for (i in DaysShows){
    print(i)
    if (i <= 2){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24),xmax=(i*24-24)+8,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i <= 9){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24)-1,xmax=(i*24-24)+9,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i >= 11){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24),xmax=(i*24-24)+12,ymin=-Inf,ymax=Inf,fill=rectgrey)
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

AddRectRestLinear<-function(ggDF2,DaysShows,rectgrey){
  
  for (i in DaysShows){
    print(i)
    if (i <= 2){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24),xmax=(i*24-24)+8,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i <= 9){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24)+1,xmax=(i*24-24)+7,ymin=-Inf,ymax=Inf,fill=rectgrey)
    }else if (i >= 11){
      ggDF2<-ggDF2+annotate("rect",xmin=(i*24-24),xmax=(i*24-24)+12,ymin=-Inf,ymax=Inf,fill=rectgrey)
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

# x<-seq(1,120)+.1
squish_trans<-function(Rdays,idx_Rdays){
  require(scales)
  trans<-function(x){
    
    if (any(is.na(x))) return(x)
    
    sf<-rep(1,length(Rdays));sf[idx_Rdays]<-.1
    transd<-cumsum(24*rep(1,length(Rdays))*sf)+(abs(Rdays[1]-1)*24)
    apfun<-approxfun(x=(Rdays)*24,y=transd,method="linear")
    
    # Max give value
    idxmax<-x>max((Rdays)*24)
    idxmin<-x<min((Rdays)*24)
    idxranged<-x<=max((Rdays)*24) & x>=min((Rdays)*24)
    x[idxmax]<-x[idxmax]-((max((Rdays)*24))-(max(transd)))
    x[idxmin]<-x[idxmin]
    x[idxranged]<-apfun(x[idxranged])
    
    return(x)
  }
  
  inv<-function(x){
    
    if (any(is.na(x))) return(x)
    
    sf<-rep(1,length(Rdays));sf[idx_Rdays]<-.1
    transd<-cumsum(24*rep(1,length(Rdays))*sf)+(abs(Rdays[1]-1)*24)
    apfun<-approxfun(y=(Rdays)*24,x=transd,method="linear")
    
    # Max give value
    idxmax<-x>max(transd)
    idxmin<-x<min(transd)
    idxranged<-x<=max(transd) & x>=min(transd)
    x[idxmax]<-max((Rdays)*24) + (x[idxmax]-max(transd))
    x[idxmin]<-x[idxmin]
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
  ord<-order(unique(Time))
  return(list(time=unique(Time)[ord],y1=predict(model,newdata=cbind.data.frame(ZT=unique(ZT)))[ord]))
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
  ord<-order(unique(Time_dataset))
  return(list(time=Time_dataset[ord],y1=predict(model,newdata=cbind.data.frame(ZT=ZTime_dataset))[ord]))
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



################### FUNCTIONS ####################
##################################################


