################## BOOTSTRAP FUNCTIONS #####################
# Use parameter of the first fit as starting position, then try a few other possibilities
Bootstrap_FasterFitFullModel<-function(swdmr,Gene,initparams=NULL){
  
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  objfun<-SWDMrGetEvalFun(model)
  paramsL<-C57BL6_TimeCourse_FullModelBoundaries()
  
  params_used<-initparams[names(paramsL$WiSd$params)]
  
  # Start with known parameters
  if (params_used$Wake > 0 | params_used$Sleep < 0){
    fit<-nlminb(params_used,objective = objfun,hessian = F,lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
    fits<-list(c(fit$par,value=fit$objective))
  }else{
    fit<-nlminb(params_used,objective = objfun,hessian = F,lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper)
    fits<-list(c(fit$par,value=fit$objective))
  }
  
  # Other possible start
  fit<-nlminb(paramsL$WiSd$params,objective = objfun,hessian = F,lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
  fit<-c(fit$par,value=fit$objective)
  fits<-c(fits,list(fit))
  
  
  fit<-nlminb(paramsL$WdSi$params,objective = objfun,hessian = F,lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper)
  fit<-c(fit$par,value=fit$objective)
  fits<-c(fits,list(fit))
  
  fits_df<-as.data.frame(do.call("rbind",fits))
  
  optimxres<-fits_df[which.min(fits_df$value),,drop=F][1,]
  
  out<-SWDMrFit(model,params = optimxres)
  stats<-SWDMrStats(model,out,detailed = T)
  optimxres$Gene<-Gene

  return(list(optimxres = optimxres,stats = stats))
}

DoBootstrap<-function(B,fit,params,swdmr,Gene,boot="samples",parametric=F,Replace=F,seed=42,speed="fast"){
  set.seed(seed+B)
  swdmr_boot<-swdmr
  
  # Bootstrap using residuals
  # Or bootstrap on samples, within time-points with replacement
  if (boot == "resid"){
    if (parametric == T){
      swdmr_boot@Gexp[,Gene]<-fit$stats$fitted+rnorm(n = length(fit$stats$residuals),mean=0,sd=sd(fit$stats$residuals))
    }else{
      swdmr_boot@Gexp[,Gene]<-fit$stats$fitted+sample(fit$stats$residuals,replace=Replace)
    }
  }else if (boot == "samples"){
    swdmr_boot@Gexp[,Gene]<-(swdmr_boot@Gexp[c(Gene,"Time")] %>% group_by(Time) %>% mutate(Resampling=sample(!!sym(Gene),replace=T)))$Resampling
  }
  
  

  if (speed=="fast"){
    res_boot<-Bootstrap_FasterFitFullModel(swdmr_boot,Gene = Gene,initparams = fit$optimxres)
  }else{
    res_boot<-FitFullModel(swdmr_boot,Gene)
  }
  
  optimxres<-res_boot$optimxres
  optimxres$Gene<-Gene
  res_boot$stats$stats$tau<-Gettau(Gene,optimxres)
  res_boot$stats$stats$SWrc<-GetSWrcMouse(Gene,swdmr_boot,fits = optimxres)
  res_boot$stats$stats$zeta<-GetZeta(Gene,optimxres)
  res_boot$stats$stats$phi_lag<-GetPhaseLag(Gene,optimxres,PosAmp = T)
  
  BootstrapFitting_res<-list(Gene=Gene,FittingType=paste("Bootstrap_",B,sep=""),
                        Wake=optimxres$Wake,Sleep=optimxres$Sleep,loggamma=optimxres$loggamma,
                        omega=optimxres$omega,AmpSin=optimxres$AmpSin,PhiSin=optimxres$PhiSin,
                        tau=res_boot$stats$stats$tau,SWrc=res_boot$stats$stats$SWrc,zeta=res_boot$stats$stats$zeta,phi_lag=res_boot$stats$stats$phi_lag,
                        RSS=res_boot$stats$stats$RSS,BF=res_boot$stats$stats$BayesFactor,
                        BIC=res_boot$stats$stats$BIC,BIC_null=res_boot$stats$stats$BIC_flat,
                        KendallTau=res_boot$stats$stats$KendalTau)
  
  
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
