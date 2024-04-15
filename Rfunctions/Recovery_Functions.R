# Functions to view and analyze recovery from SD

ExpressionAndEffectSize<-function(Gene,Dataset,centerSW=T){
  gg1<-PlotMouseFit(Gene,Dataset,
                    ASContribution = F,
                    bsldashed = T,
                    xlim = c(24,240),
                    AddSleepWake=T,
                    PointSize = .5,
                    lineSize = .75,
                    MeanPointSize = .75,
                    lineSizeBSL = .75)
  
  gg2<-PlotFoldChange(Gene,Dataset,xlim=c(24,96),centerSW=T)
  
  gg1<-RemoveXaxis(gg1)+ 
    theme(plot.margin = unit(c(0,1,0,0), "mm"))+
    scale_x_continuous(position = "top")+
    coord_cartesian(xlim=c(24,96))+
    ylab("Expression\n[log2 CPM]")+
    theme(text = element_text(size = 8))+
    ggtitle(bquote(italic(.(Gene))))
  
  gg2<-gg2 + 
    theme(plot.margin = unit(c(0,1,0,0), "mm"))+
    coord_cartesian(xlim=c(24,96))+
    ylab("Effect Size")+
    theme(text = element_text(size = 8))
  
  # cowplot::plot_grid(gg1, 
  #                    gg2, ncol = 1, align = "v",axis="b")
  gg1/gg2
}


EffectSize_RecSleepReplaced<-function(Gene,Dataset,doylim=c(-3,3)){
  library(patchwork)
  cols<-ColorCode()
  
  dd<-ViewModelFoldChange(ProbeGene = Gene,Dataset = Dataset,xlimPlot=c(0,48))
  NicePlot<-function(dd,doxlim=NULL,title,CenterSW=T){
    # as df
    ddf<-cbind.data.frame(Time=dd$Time,FC=dd$FC,RecoverySleep=dd$BSL_sleep_start)
    gg<-ggplot(aes(Time,FC),data=ddf)
    idx<-dd$BSL_sleep_start == T
    gg<-gg+geom_line(color=colorCode["Modelfit"],size=0.75)+theme_classic()+scale_x_continuous(breaks=seq(0,72,by=12),labels=seq(0,72,by=12)+54)
    gg<-gg+annotate("line",dd$Time[idx],dd$Tau1[idx],color="red",size=0.5,linetype="22")
    if (dd$Zeta<1){
      gg<-gg+annotate("line",dd$Time[idx],dd$Tau2[idx],color="red",size=0.5,linetype="22")
    }
    gg<-gg+theme(legend.position = "none")
    gg<-gg+scale_y_continuous(limits=doylim)
    # gg<-gg+geom_hline(yintercept = c(-1,1),linetype="22",col="red")
    gg<-gg + annotate("rect",xmin = -6,ymin = -Inf,xmax=0,ymax=Inf,fill=rgb(1,0,0,.3))
    gg<-gg+coord_cartesian(xlim=doxlim)
    gg<-gg+geom_vline(xintercept = dd$res,color="steelblue",linetype="solid")
    gg<-gg+ylab("Effect Size") +theme(text=element_text(size = 8),plot.title = element_text(size=8),plot.margin = unit(c(0,1,0,0), "mm"))
    gg<-gg + annotate("text",x=36,y=2,label=title,hjust=1)
    
    # Add SW
    SWdiff_Sleepbsl<-dd$SWdist[dd$SWdist$Time>24 & dd$SWdist$Time<=48,]
    SWdiff_Sleeprec<-dd$SWdist[dd$SWdist$Time>48 & dd$SWdist$Time<=120,]
    SWdiff_Sleepdiff<-SWdiff_Sleeprec
    SWdiff_Sleepdiff$Sleep<-SWdiff_Sleeprec$Sleep-rep(SWdiff_Sleepbsl$Sleep,3)
    Tsleepdiff<-seq(54.1,120,by=.1)
    
    SWdiff_Sleepdiff<-SWdiff_Sleepdiff %>% group_by(ceiling(Time)) %>% mutate_at("Sleep",mean)
    SWdiff_Sleepdiff<-SWdiff_Sleepdiff$Sleep#/.1
    SWdiff_Sleepdiff<-SWdiff_Sleepdiff[seq(48.1,120,by=.1)>=54.1]
    
    yrange<-layer_scales(gg)$y$range$range
    if (CenterSW == T){
      yrange[2]<-max(abs(yrange))
      yrange[1]<- -1*yrange[2]
    }
    
    
    SWdiff_Sleepdiff_yr<-((SWdiff_Sleepdiff - (-1)) / (2)) *(yrange[2]*2 - yrange[1]*2) + yrange[1]*2
    
    gg<-gg+annotate("ribbon",x=Tsleepdiff-54,ymax=SWdiff_Sleepdiff_yr,ymin=SWdiff_Sleepdiff_yr[length(SWdiff_Sleepdiff_yr)],
                    fill=alpha(cols[["SWforce"]],.2))+
      annotate("line",x=Tsleepdiff-54,y=SWdiff_Sleepdiff_yr,color=alpha(cols[["SWforce"]],.7))
    
    gg<-gg+xlab("Time [h]")
    
    gg
  }
  
  NicePlot6<-function(dd,title){
    ddf1<-cbind.data.frame(RecSleep=dd[[1]],Time1SD=dd[[2]]+54)
    ddf2<-cbind.data.frame(TimeLoessfit=seq(0,42,by=.1),LoessFit=dd[[3]])
    gg<-ggplot(aes(RecSleep,Time1SD),data=ddf1)+geom_point(size=.75)
    gg<-gg+scale_x_continuous(breaks=seq(0,42,by=6))
    #gg<-gg+annotate("path",ddf2$TimeLoessfit,ddf2$LoessFit,col="red")
    gg<-gg+ylab("Time")+xlab("Rec. Sleep duration [h]")
    gg<-gg+theme_classic()+ggtitle(title)+theme(legend.position = "none",text=element_text(size = 8),plot.title = element_text(size=8))
    gg
  }
  
  doxlim<-c(0,36)
  # gff<-((RemoveXlab(NicePlot(dd[[2]],doxlim,title="18h Rec. sleep"))  |  RemoveXlab(RemoveYaxis(NicePlot(dd[[3]],doxlim,title="12h  Rec. sleep")))) /
  #            ((NicePlot(dd[[4]],doxlim,title="6h Rec. sleep")) |  RemoveYaxis(NicePlot(dd[[5]],doxlim,title="0h Rec. sleep"))))
  
  # Title
  p_Title<-ggplot() + 
    annotate("text",
             x = 1,
             y = 1,
             size = 0.01,
             label = "")+
    theme_void()+ggtitle("Recovery sleep duration:")+theme(plot.margin = margin(0,0,0,0),plot.title = element_text(hjust = 0.5,size=12))
  
  gff<-(p_Title/(RemoveXaxis(NicePlot(dd[[5]],doxlim,title="0h"))  |  RemoveXaxis(RemoveYaxis(NicePlot(dd[[4]],doxlim,title="6h")))) /
          ((NicePlot(dd[[3]],doxlim,title="12h")) |  RemoveYaxis(NicePlot(dd[[2]],doxlim,title="18h"))))+plot_layout(heights = c(1,50,50))
  
  gff2<-NicePlot6(dd[[6]],"Time to reach 1\nEffect size")
  
  return(list(gff,gff2))
  
}

GetFoldChange<-function(ProbeGene,Dataset="Cortex"){
  # Load data given dataset
  data<-LoadData(Dataset)
  
  model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)
  
  # Simulate expression as Z-score
  out_response<-GetModelFitted_Zscore(ProbeGene,Dataset,data,model)
  
  # Fold-change after SD
  DiffSWresponse<-GetDiffSWresponse(out_response,Dataset)
  
  return(DiffSWresponse)
  
}

GetFCendSD<-function(ProbeGene,Dataset,data){
  
  model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)
  
  if (Dataset %in% c("Liver","Cortex")){
    param_model<-data$fits[data$fits$Gene == ProbeGene,][1,]
    intercept<-model@intercept
  }
  
  out<-SWDMrFit(model,param_model)
  
  FC<-out$y1[out$time == 54]-out$y1[out$time == 30]

  sdexpr<-sqrt(param_model$RSS/(nrow(model@Gexp)-length(SWDMr:::GetFreeFixedParams(model)$FreeParams$subparameter)))
  return(FC/sdexpr)
  
  # model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)
  # 
  # FC<-mean(data$rna_expr[,ProbeGene][data$rna_expr$Time == 54])-mean(data$rna_expr[,ProbeGene][data$rna_expr$Time == 30])
  # if (Dataset %in% c("Liver","Cortex")){
  #   param_model<-data$fits[data$fits$Gene == ProbeGene,][1,]
  #   intercept<-model@intercept
  # }
  # sdexpr<-sqrt(param_model$RSS/(nrow(model@Gexp)-length(SWDMr:::GetFreeFixedParams(model)$FreeParams$subparameter)))
  # return(FC/sdexpr)
  # FCdata<-GetFoldChange(ProbeGene,Dataset)
  #return(FCdata$FC[1])
}


GetmaxFCafterSD<-function(ProbeGene,Dataset,data){
  model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)
  
  if (Dataset %in% c("Liver","Cortex")){
    param_model<-data$fits[data$fits$Gene == ProbeGene,][1,]
    intercept<-model@intercept
  }
  
  out<-SWDMrFit(model,param_model)
  
  FC<-(out$y1[out$time > 54 & out$time<= 72]-intercept)/sdexpr-
    (out$y1[out$time > 30 & out$time<= 48]-intercept)/sdexpr
  FC<-max(abs(FC))
  
  sdexpr<-sqrt(param_model$RSS/(nrow(model@Gexp)-length(SWDMr:::GetFreeFixedParams(model)$FreeParams$subparameter)))
  return(FC/sdexpr)
  # FCdata<-GetFoldChange(ProbeGene,Dataset)
  #return(FCdata$FC[1])
}

PlotFoldChange<-function(ProbeGene,Dataset="Cortex",xlim=c(54,96),centerSW=T){
  require(dplyr)
  cols<-ColorCode()
  
  FCdata<-GetFoldChange(ProbeGene,Dataset)
  
  dd<-cbind.data.frame(FC=FCdata$FC,Time=FCdata$time)
  gg<-ggplot(aes(x=Time,y=FC),data=dd)+theme_classic()
  gg<-gg + scale_x_continuous(trans = squish_trans(108, 203, 30),
                              breaks = c(seq(0, 96, by = 12),204,216,228))
  gg<-gg + annotate("rect",xmin = 48,ymin = -Inf,xmax=54,ymax=Inf,fill=rgb(1,0,0,.3),color=rgb(1,0,0,.3))
  gg<-gg + coord_cartesian(xlim=xlim)
  gg<-AddLDblocks(gg,starts=seq(60,222,by=24))
  gg<-gg+ylab("Fold-Change [Z-score]")+xlab("Time [h]")
  gg<-gg + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA)) 
  
  # Add sleep-wake diff
  data<-LoadData(Dataset)
  SWdiff_Sleepbsl<-data$SWdf[data$SWdf$Time>24 & data$SWdf$Time<=48,]
  SWdiff_Sleeprec<-data$SWdf[data$SWdf$Time>48 & data$SWdf$Time<=120,]
  SWdiff_Sleepdiff<-SWdiff_Sleeprec
  SWdiff_Sleepdiff$Sleep<-SWdiff_Sleeprec$Sleep-rep(SWdiff_Sleepbsl$Sleep,3)
  Tsleepdiff<-seq(54.1,120,by=.1)
  
  SWdiff_Sleepdiff<-SWdiff_Sleepdiff %>% group_by(ceiling(Time)) %>% mutate_at("Sleep",mean)
  SWdiff_Sleepdiff<-SWdiff_Sleepdiff$Sleep#/.1
  SWdiff_Sleepdiff<-SWdiff_Sleepdiff[seq(48.1,120,by=.1)>=54.1]
  
  gg2<-gg +geom_line(size=.75,color=cols[["Modelfit"]])
  yrange<-layer_scales(gg2)$y$range$range
  
  if (centerSW == T){
    yrange[2]<-max(abs(yrange))
    yrange[1]<- -1*yrange[2]
    SWdiff_Sleepdiff_yr<-((SWdiff_Sleepdiff - (-1)) / (2)) *(yrange[2]*2 - yrange[1]*2) + yrange[1]*2
  }else{
    idx_min<-which.min(yrange)
    idx_max<-which.max(yrange)
    yrange[idx_min]<-min(abs(yrange))
    yrange[idx_max]<- -1*yrange[idx_min]
    SWdiff_Sleepdiff_yr<-((SWdiff_Sleepdiff - (-1)) / (2)) *(yrange[2]*2 - yrange[1]*2) + yrange[1]*2
  }
  
  
  
  gg<-gg+annotate("ribbon",x=Tsleepdiff,ymax=SWdiff_Sleepdiff_yr,ymin=SWdiff_Sleepdiff_yr[length(SWdiff_Sleepdiff_yr)],
           fill=alpha(cols[["SWforce"]],.2))+
  annotate("line",x=Tsleepdiff,y=SWdiff_Sleepdiff_yr,color=alpha(cols[["SWforce"]],.7))
  
  gg<-gg +geom_line(size=.75,color=cols[["Modelfit"]])
  
  return(gg)
  
}


ViewModelFoldChange<-function(ProbeGene,Dataset="Cortex",xlimPlot=NULL){
  
  par(mfrow=c(2,5))
  
  # Load data given dataset
  data<-LoadData(Dataset)
  # Color Codes 
  cols<-ColorCode()
  
  # Get Model
  model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)

  # Simulate expression as Z-score
  out_response<-GetModelFitted_Zscore(ProbeGene,Dataset,data,model)
  DiffSWresponse<-GetDiffSWresponse(out_response,Dataset) # Fold-change after SD
  
  # some ranges for plots
  ranges<<-c(min(DiffSWresponse$FC),max(DiffSWresponse$FC))
  xlimPlot<<-xlimPlot
  
  # If FC is below a theshold, take this value:
  PrevThresh<<-0
  
  # Get Decay lines, save results
  dd<-list()
  dd[[1]]<-GetDecayLine(DiffSWresponse,data,HourOfRECsleep = 42,DoPolot=T,model=model)
  rm(DiffSWresponse)
  ################
  
  PlotDecay<-function(HourOfRECsleep=0,DoPolot=F){
    
    model<-GetModel(Dataset = Dataset,data = data,ProbeGene = ProbeGene)
    
    # Avoid taking SD
    if (HourOfRECsleep>18){
      BSLSleepToTake<- HourOfRECsleep - 24
    }else{
      BSLSleepToTake<-HourOfRECsleep
    }
    
    idx_bsl<-model@SWdist$Time>(BSLSleepToTake+6) & model@SWdist$Time<=((BSLSleepToTake+6)+24)
    idx_rec1<-model@SWdist$Time>(48+(HourOfRECsleep+6)) & model@SWdist$Time<=(48+(HourOfRECsleep+6)+24)
    idx_rec2<-model@SWdist$Time>(48+(HourOfRECsleep+6)+24) & model@SWdist$Time<=(48+(HourOfRECsleep+6)+48)
    
    model@SWdist[idx_rec1,]<-model@SWdist[idx_bsl,]
    model@SWdist[idx_rec2,]<-model@SWdist[idx_bsl,]
    model@SWdist$Time<-round(seq(model@SWdist$Time[1],tail(model@SWdist$Time,1),by=.1),1)
    model<-SetParametersModel(model)
    # Simulate expression as Z-score
    out_response<-GetModelFitted_Zscore(ProbeGene,Dataset,data,model)
    
    # Fold-change after SD
    DiffSWresponse<-GetDiffSWresponse(out_response,Dataset)
    
    # Get Decay lines
    res<-GetDecayLine(DiffSWresponse = DiffSWresponse,data,HourOfRECsleep = HourOfRECsleep,DoPolot=DoPolot,model=model)
    rm(DiffSWresponse)
    return(res)
  }
  
  
  
  PlotDecay(HourOfRECsleep=30,DoPolot=T)
  PlotDecay(HourOfRECsleep=25,DoPolot=T)
  PlotDecay(HourOfRECsleep=20,DoPolot=T)
  dd[[2]]<-PlotDecay(HourOfRECsleep=18,DoPolot=T)
  dd[[3]]<-PlotDecay(HourOfRECsleep=12,DoPolot=T)
  dd[[4]]<-PlotDecay(HourOfRECsleep=6,DoPolot=T)
  #PlotDecay(HourOfRECsleep=4,DoPolot=T)
  PlotDecay(HourOfRECsleep=2,DoPolot=T)
  dd[[5]]<-PlotDecay(HourOfRECsleep=0,DoPolot=T)
  
  Trec_res<-as.numeric(sapply(seq(0,42,by=1),function(x){PlotDecay(x)}))
  Rec_TP<-seq(0,42,by=1)
  #plot(Rec_TP,Trec_res,main=ProbeGene)
  mod1<-loess(Trec_res~Rec_TP,span=.5)
  yfit1=predict(mod1,newdata=seq(0,42,by=.1))
  #points(seq(0,42,by=.1),yfit1,type="l",lwd=2,col="red")
  dd[[6]]<-list(Rec_TP,Trec_res,yfit1)
  
  
  return(dd)
}


GetDecayLine<-function(DiffSWresponse,data,HourOfRECsleep=42,DoPolot=F,model=model){
  
  # Threshold
  FCThresh<-1
  
  # Get Zeta value (Overdamped or Underdamped)
  Zeta<-GetZeta(DiffSWresponse$param_model$Gene,data$fits)
  
  if (Zeta < 1){
    # Time constant
    tau <- 2/exp(DiffSWresponse$param_model$loggamma)
    
    # Get amplitude
    idx<-DiffSWresponse$TimeSinceEndSD>=(HourOfRECsleep)
    Time<-DiffSWresponse$TimeSinceEndSD[idx]-DiffSWresponse$TimeSinceEndSD[idx][1]
    FC<-DiffSWresponse$FC[idx]
    
    # Get start amp phase
    omega<-DiffSWresponse$param_model$omega

    # Use linear model to help start fitting; accurate if initial speed = 0, but not always the case
    CosP<-cos(omega*Time)*exp(-(Time)/tau)
    SinP<-sin(omega*Time)*exp(-(Time)/tau)
    model_exp<-lm(FC~0+CosP+SinP)
    cm<-coef(model_exp)
    Amp<-sqrt(cm[1]^2+cm[2]^2)
    
    model_exp_nls<-nls(FC~exp(-Time/tau)*A*sin(w*Time+phi),start=list(A=Amp,phi=atan2(cm[1],cm[2]),w=omega),algorithm = "port",control = list(warnOnly=T))
    Amp<-coef(model_exp_nls)["A"]
    #plot(Time,FC);lines(Time,fitted(model_exp_nls),col="red")
    
    Time<-DiffSWresponse$TimeSinceEndSD-DiffSWresponse$TimeSinceEndSD[idx][1]
    #app<-approxfun(cm["A"]*exp(-(Time)/tau),DiffSWresponse$TimeSinceEndSD)
    app<-approxfun(Amp*exp(-(Time)/tau),DiffSWresponse$TimeSinceEndSD)
    
    
    res<-c(app(FCThresh),app(-FCThresh))
    if (max(FC)<FCThresh & min(FC)> -1*FCThresh){
      res<-PrevThresh
    }else{
      PrevThresh<<-res[!is.na(res)]
      res<-res[!is.na(res)]
    }
    
    if (DoPolot == T){
      plot(DiffSWresponse$TimeSinceEndSD[DiffSWresponse$TimeSinceEndSD>HourOfRECsleep],
           DiffSWresponse$FC[DiffSWresponse$TimeSinceEndSD>HourOfRECsleep],
           ylim=ranges,main=HourOfRECsleep,xlim=xlimPlot,pch=19,
           type="l",lwd=3,col="black")
      lines(DiffSWresponse$TimeSinceEndSD[DiffSWresponse$TimeSinceEndSD<HourOfRECsleep],
           DiffSWresponse$FC[DiffSWresponse$TimeSinceEndSD<HourOfRECsleep],
           ylim=ranges,main=HourOfRECsleep,xlim=xlimPlot,pch=19,
           type="l",lwd=3,col="green")#col=c(rgb(0,1,0,.5),"black")[factor(DiffSWresponse$TimeSinceEndSD>HourOfRECsleep,levels=c("FALSE","TRUE"))]
      lines(DiffSWresponse$TimeSinceEndSD,
            Amp*exp(-(Time)/tau),
            #cm["A"]*exp(-(Time)/tau),
            col="red")
      lines(DiffSWresponse$TimeSinceEndSD,
            #-1*cm["A"]*exp(-(Time)/tau),
            -1*Amp*exp(-(Time)/tau),
            col="red")
      # lines(DiffSWresponse$TimeSinceEndSD,
      #       cm["A"]*exp(-(Time)/tau)*sin(cm["omega"]*(Time)+cm["phi"]),col="red")
      abline(h=FCThresh,v=c(app(FCThresh),app(-FCThresh)))
      
      dd<-list(Time=DiffSWresponse$TimeSinceEndSD,
               FC=DiffSWresponse$FC,
               BSL_sleep_start=DiffSWresponse$TimeSinceEndSD>HourOfRECsleep,
               Zeta=Zeta,
               Tau1=Amp*exp(-(Time)/tau),
               Tau2=-1*Amp*exp(-(Time)/tau),
               res=res,
               SWdist=model@SWdist)
      return(dd)
      
    }
    
    return(res)
    
    
    
  }else{
    gamma<-exp(DiffSWresponse$param_model$loggamma)
    omega<-DiffSWresponse$param_model$omega
    
    idx<-DiffSWresponse$TimeSinceEndSD>HourOfRECsleep
    Time<-DiffSWresponse$TimeSinceEndSD[idx]-DiffSWresponse$TimeSinceEndSD[idx][1]
    FC<-DiffSWresponse$FC[idx]
    
    C1part<-exp( (-gamma/2*Time) + sqrt( (gamma/2)^2-omega^2)*Time )
    C2part<-exp( (-gamma/2*Time) + -sqrt( (gamma/2)^2-omega^2)*Time )
    
    model_exp<-lm(FC~0+C1part+C2part)
    cm<-coef(model_exp)
    Time<-DiffSWresponse$TimeSinceEndSD-DiffSWresponse$TimeSinceEndSD[idx][1]
    C1part<-exp( (-gamma/2*Time) + sqrt( (gamma/2)^2-omega^2)*Time )
    C2part<-exp( (-gamma/2*Time) + -sqrt( (gamma/2)^2-omega^2)*Time )
    app<-approxfun(cm[1] * C1part,# +
                   #cm[2] * C2part,
                   DiffSWresponse$TimeSinceEndSD)
    
    res<-c(app(FCThresh),app(-FCThresh))
    if (max(FC)<FCThresh & min(FC)> -1*FCThresh){
      res<-PrevThresh
    }else{
      PrevThresh<<-res[!is.na(res)]
      res<-res[!is.na(res)]
    }
    
    if (DoPolot == T){
      # plot(DiffSWresponse$TimeSinceEndSD[DiffSWresponse$TimeSinceEndSD>HourOfRECsleep],
      #      DiffSWresponse$FC[DiffSWresponse$TimeSinceEndSD>HourOfRECsleep],
      #      ylim=ranges,main=HourOfRECsleep,xlim=xlimPlot,pch=19,
      #      type="l",lwd=3,col="black")
      # lines(DiffSWresponse$TimeSinceEndSD[DiffSWresponse$TimeSinceEndSD<HourOfRECsleep],
      #       DiffSWresponse$FC[DiffSWresponse$TimeSinceEndSD<HourOfRECsleep],
      #       ylim=ranges,main=HourOfRECsleep,xlim=xlimPlot,pch=19,
      #       type="l",lwd=3,col="green")#col=c(rgb(0,1,0,.5),"black")[factor(DiffSWresponse$TimeSinceEndSD>HourOfRECsleep,levels=c("FALSE","TRUE"))]
      # lines(DiffSWresponse$TimeSinceEndSD,
      #       cm[2] * C2part,
      #       col="blue")
      # lines(DiffSWresponse$TimeSinceEndSD,
      #       cm[1] * C1part,
      #       col="red")
      # abline(h=FCThresh,v=c(app(FCThresh),app(-FCThresh)))

      dd<-list(Time=DiffSWresponse$TimeSinceEndSD,
               FC=DiffSWresponse$FC,
               BSL_sleep_start=DiffSWresponse$TimeSinceEndSD>HourOfRECsleep,
               Zeta=Zeta,
               Tau1=cm[1] * C1part,
               Tau2=cm[2] * C2part,
               res=res,
               SWdist=model@SWdist)
      return(dd)
      
    }

    return(res)
    
    
  }
  
}


# Expression after SD (T54, expression goes up to T222, 7 days)
GetDiffSWresponse<-function(fitted,Dataset){
  if (Dataset %in% c("Liver","Cortex")){
    DiffSWresponse<-list()
    fitted_ASD<-fitted$y1[fitted$time > 54.0 & fitted$time <= 222]
    # Replicate baseline oscillation for 7 days
    fitted_BSL<-rep(fitted$y1[fitted$time > 6.0 & fitted$time <= 30.0],7)
    # Differential of sleep-wake response
    DiffSWresponse$FC<-fitted_ASD - fitted_BSL
    DiffSWresponse$time<-seq(54.1,222,by=.1)
    DiffSWresponse$TimeSinceEndSD<-DiffSWresponse$time-54
    DiffSWresponse$param_model<-fitted$param_model
    return(DiffSWresponse)
  }
}


######## Simulate expression
GetModelFitted_Zscore<-function(ProbeGene,Dataset,data,model){
  if (Dataset %in% c("Liver","Cortex")){
    param_model<-data$fits[data$fits$Gene == ProbeGene,][1,]
    intercept<-model@intercept
  }
  if (Dataset %in% c("Desync","Res","Ctr")){
    param_model<-data$fits[data$fits$ProbeName == ProbeGene,][1,]
    if (Dataset == "Desync"){
      intercept<-model@intercept
    }else{
      intercept<-data$fits[data$fits$ProbeName == ProbeGene,"intercept"][1,]
    }
  }
  out_response<-SWDMrFit(model,params = param_model,method="Solve")
  
  # Get intercept of the model and standard error of residual to compute effect size
  meanexpr<-model@intercept
  sdexpr<-sqrt(param_model$RSS/(nrow(model@Gexp)-length(SWDMr:::GetFreeFixedParams(model)$FreeParams$subparameter)))
  out_response$y1<-(out_response$y1-meanexpr)/sdexpr # scale fitted value
  out_response$param_model<-param_model
  
  return(out_response)
}


