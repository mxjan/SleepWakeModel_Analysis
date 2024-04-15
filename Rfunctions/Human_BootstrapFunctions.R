################## BOOTSTRAP FUNCTIONS #####################
# Use parameter of the first fit as starting position, then try a few other possibilities
Bootstrap_FasterFitDataFull<-function(ProbeID, MeanSWdf, HT_Desync, MeanSWdfExt, HT_Ext, MeanSWdfRest, HT_Res,initparams){

  # Get intercept desynchronization protocols
  Expr<-HT_Desync[HT_Desync$Time < 72,ProbeID]
  Time<-HT_Desync[HT_Desync$Time < 72,"Time"]
  MeanGeneExprInBaseline<-coef(lm(Expr~sin(2*pi/24*Time) + cos(2*pi/24*Time) ))[["(Intercept)"]]
  
  # Build models
  model_desync<-Human_TimeCourse_GetFullModel(MeanSWdf,HT_Desync,ProbeID,intercept=MeanGeneExprInBaseline)
  model_ext<-Human_TimeCourse_GetFullModel(MeanSWdfExt,HT_Ext,ProbeID)
  model_res<-Human_TimeCourse_GetFullModel(MeanSWdfRest,HT_Res,ProbeID)
  
  # Prepare objective function
  objFdesync<-SWDMrGetEvalFun(model_desync)
  objFext<-SWDMrGetEvalFun(model_ext)
  objFres<-SWDMrGetEvalFun(model_res)
  
  # Prepare models limits
  paramsL<-Human_TimeCourse_FullModelBoundaries(model_desync,model_ext,model_res)
  
  params_used<-initparams[names(paramsL$WiSd$params)]
  
  # Combine RSS for all fits
  objFAll<-function(params){
    paramsdesync<-params[c(1,2,3,4,5,6)]
    paramsextres<-params[c(1,2,3,4,5,6,7)]
    
    RSS_desync<-objFdesync(paramsdesync)
    RSS_ext<-objFext(paramsextres)
    RSS_res<-objFres(paramsextres)
    
    RSSAll<-RSS_desync+RSS_ext+RSS_res
    
    return(RSSAll)
    
  }
  
  # Start with known parameters
  if (params_used$Wake > 0 | params_used$Sleep < 0){
    fit_boot<-nlminb(params_used,objective = objFAll,hessian = F,lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
    fits<-list(c(fit_boot$par,value=fit_boot$objective))
    for (i in 1:nrow(paramsL$WiSd$params_multi)){
      fit_boot<-nlminb(paramsL$WiSd$params_multi[i,],objective = objFAll,hessian = F,lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
      fit_boot<-c(fit_boot$par,value=fit_boot$objective)
      fits<-c(fits,list(fit_boot))
    }
  }else{
    fit_boot<-nlminb(params_used,objective = objFAll,hessian = F,lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper)
    fits<-list(c(fit_boot$par,value=fit_boot$objective))
    for (i in 1:nrow(paramsL$WdSi$params_multi)){
      fit_boot<-nlminb(paramsL$WdSi$params_multi[i,],objective = objFAll,hessian = F,lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper)
      fit_boot<-c(fit_boot$par,value=fit_boot$objective)
      fits<-c(fits,list(fit_boot))
    }
  }
  
  fits_df<-as.data.frame(do.call("rbind",fits))
  optimxres<-fits_df[which.min(fits_df$value),,drop=F][1,]
  
  # Split between desync and rest-ext sleep (intercept change)
  resdesync<-optimxres[,c(1,2,3,4,5,6)]
  resextres<-optimxres[,c(1,2,3,4,5,6,7)]
  
  # Get Statistics
  stats<-GetStatsFull(resdesync,resextres,model_desync,model_ext,model_res,ProbeID)
  
  return(list(paramsDesync=resdesync,paramsExtRes=resextres,stats=stats,
              model_desync=model_desync,model_ext=model_ext, model_res = model_res,
              resid_desync=stats$resid_desync,fitted_desync=stats$fitted_desync,
              resid_rest=stats$resid_rest,fitted_rest=stats$fitted_rest,
              resid_ext=stats$resid_ext,fitted_ext=stats$fitted_ext))
  
}

DoBootstrap<-function(B,fit,params,Gene,ProbeID,boot="samples",parametric=F,Replace=F,seed=42,speed="fast"){
  set.seed(seed+B)
  HT_Desync_boot<-HT_Desync[,c(ProbeID,"Time")]
  HT_Ext_boot<-HT_Ext[,c(ProbeID,"Time")]
  HT_Res_boot<-HT_Res[,c(ProbeID,"Time")]
  
  if (boot == "resid"){
    if (parametric == T){
      rand_resid_desync<-fit$fitted_desync +rnorm(n = length(fit$resid_desync),mean=0,sd=sd(fit$resid_desync))
      rand_resid_rest<-fit$fitted_rest +rnorm(n = length(fit$resid_rest),mean=0,sd=sd(fit$resid_rest))
      rand_resid_ext<-fit$fitted_ext +rnorm(n = length(fit$resid_ext),mean=0,sd=sd(fit$resid_ext))
    }else{
      rand_resid_desync<-fit$fitted_desync +sample(fit$resid_desync,replace=Replace)
      rand_resid_rest<-fit$fitted_rest +sample(fit$resid_rest,replace=Replace)
      rand_resid_ext<-fit$fitted_ext +sample(fit$resid_ext,replace=Replace)
    }
    
    NAdesync<-is.na(HT_Desync_boot[,ProbeID])
    if (any(NAdesync)){
      HT_Desync_boot[,ProbeID][!NAdesync]<-rand_resid_desync
    }else{
      HT_Desync_boot[,ProbeID]<-rand_resid_desync
    }
    
    NAext<-is.na(HT_Ext_boot[,ProbeID])
    if (any(NAext)){
      HT_Ext_boot[,ProbeID][!NAext]<-rand_resid_ext
    }else{
      HT_Ext_boot[,ProbeID]<-rand_resid_ext
    }
    
    NArest<-is.na(HT_Res_boot[,ProbeID])
    if (any(NArest)){
      HT_Res_boot[,ProbeID][!NArest]<-rand_resid_rest
    }else{
      HT_Res_boot[,ProbeID]<-rand_resid_rest
    }
  }else if (boot == "samples"){
    
    NAdesync<-is.na(HT_Desync_boot[,ProbeID])
    if (any(NAdesync)){
      HT_Desync_boot[,ProbeID][!NAdesync]<-(HT_Desync_boot[!NAdesync,] %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }else{
      HT_Desync_boot[,ProbeID]<-(HT_Desync_boot %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }
    
    NAext<-is.na(HT_Ext_boot[,ProbeID])
    if (any(NAext)){
      HT_Ext_boot[,ProbeID][!NAext]<-(HT_Ext_boot[!NAext,] %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }else{
      HT_Ext_boot[,ProbeID]<-(HT_Ext_boot %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }
    
    NArest<-is.na(HT_Res_boot[,ProbeID])
    if (any(NArest)){
      HT_Res_boot[,ProbeID][!NArest]<-(HT_Res_boot[!NArest,] %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }else{
      HT_Res_boot[,ProbeID]<-(HT_Res_boot %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(ProbeID),replace=T)))$Resampling
    }
  }else if (boot == "TSsamples"){
    
    GSM_sam_desync<-unlist(sapply(sample(unique(GSE48113_metaD$`subject:ch1`),replace=T),function(ID){rownames(GSE48113_metaD[GSE48113_metaD$`subject:ch1` == ID,])}))
    HT_Desync_boot<-HT_Desync_boot[GSM_sam_desync,]
    
    GSE39445_metaD_res<-GSE39445_metaD[GSE39445_metaD$`sleepprotocol:ch1` == "Sleep Restriction",]
    GSM_sam_res<-unlist(sapply(sample(unique(GSE39445_metaD_res$`subject:ch1`),replace=T),function(ID){rownames(GSE39445_metaD_res[GSE39445_metaD_res$`subject:ch1` == ID,])}))
    HT_Res_boot<-HT_Res_boot[GSM_sam_res,]
    
    GSE39445_metaD_ext<-GSE39445_metaD[GSE39445_metaD$`sleepprotocol:ch1` == "Sleep Extension",]
    GSM_sam_ext<-unlist(sapply(sample(unique(GSE39445_metaD_ext$`subject:ch1`),replace=T),function(ID){rownames(GSE39445_metaD_ext[GSE39445_metaD_ext$`subject:ch1` == ID,])}))
    HT_Ext_boot<-HT_Ext_boot[GSM_sam_ext,]
  }

  if (speed=="fast"){
    res_boot<-Bootstrap_FasterFitDataFull(ProbeID , MeanSWdf, HT_Desync_boot, MeanSWdfExt, HT_Ext_boot, MeanSWdfRest, HT_Res_boot,initparams = fit$paramsExtRes)
  }else{
    res_boot<-FitDataFull(ProbeID , MeanSWdf, HT_Desync_boot, MeanSWdfExt, HT_Ext_boot, MeanSWdfRest, HT_Res_boot)
  }
  
  optimxres<-res_boot$paramsExtRes
  optimxres$Gene<-Gene
  optimxres$ProbeName<-ProbeID
  res_boot$stats$tau<-Gettau(ProbeID,optimxres)
  res_boot$stats$SWrc<-GetSWrcHuman(ProbeID,MeanSWdf,HT_Desync = HT_Desync,fits = optimxres)
  res_boot$stats$zeta<-GetZeta(ProbeID,optimxres)
  res_boot$stats$phi_lag<-GetPhaseLag(ProbeID,optimxres,PosAmp = T)
  BootstrapFitting_res<-list(ProbeName=ProbeID,Gene=Gene,FittingType=paste("Bootstrap_",B,sep=""),
                        Wake=optimxres$Wake,Sleep=optimxres$Sleep,loggamma=optimxres$loggamma,
                        omega=optimxres$omega,AmpSin=optimxres$AmpSin,PhiSin=optimxres$PhiSin,intercept=optimxres$intercept,
                        tau=res_boot$stats$tau,SWrc=res_boot$stats$SWrc,zeta=res_boot$stats$zeta,phi_lag=res_boot$stats$phi_lag,
                        RSS=res_boot$stats$RSS_model,BF=res_boot$stats$BF,
                        BIC=res_boot$stats$BIC_model,BIC_null=res_boot$stats$BIC_flat,
                        KendallTau=res_boot$stats$KendallTau)
  
  return(BootstrapFitting_res)
}

PlotBootstrap<-function(res,CI=c(0.025,0.975)){
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  
  PFunc<-function(x){
    tmp_dd<-dd[dd$variable == x,]
    quant<-quantile(tmp_dd$value,prob=CI)
    x.dens <- density(tmp_dd$value)
    df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
    gg<-ggplot(aes(value),data=tmp_dd)+geom_histogram(aes(y = after_stat(density)))
    gg<-gg + geom_density() + ggtitle(paste(x," CI [",signif(quant[1],2),";",signif(quant[2],2),"]",sep="")) + theme_bw()
    gg <-  gg + geom_area(data = subset(df.dens, x >= quant[1] & x <= quant[2]), 
                          aes(x=x,y=y), fill = alpha("red",.25))
  }
  res$zeta<-log10(res$zeta)
  dd<-melt(res[,c("Wake","Sleep","loggamma","omega","AmpSin","PhiSin","RSS","tau","SWrc","zeta","phi_lag")])
  # pl<-lapply(c("Wake","Sleep","loggamma","omega","AmpSin","PhiSin"),PFunc)
  # print(wrap_plots(pl))
  pl<-lapply(c("tau","SWrc","zeta","phi_lag"),PFunc)
  print(wrap_plots(pl))
}
