###################################################
# - Functions to build model for human transcriptomics data
# - Data fit merge the datasets of extended, restricted sleep 
#   and desynchronization protocols
# 
#



###################################################
########### FULL MODEL #########

Human_TimeCourse_GetFullModel<-function(SWdf,Exprdf,ProbeID,intercept=NULL){
 
  swdmr<-SWDMr(SWdist=SWdf, Gexp=Exprdf[,c(ProbeID,"Time")])
  
  ############# Build model ###############
  model<-initDDHOmodel(swdmr, VarExp = ProbeID, UseDampingRatio = F)
  
  # Set an intercept if given
  if (! is.null(intercept)){
    model<-FixIntercept(model,intercept)
  }

  # We replicate baseline for some days
  model<-ReplicateDrivingForce(model,c(12.0,36.0),20)
  
  # Add sleep-wake force
  model<-AddForce(model,"Wake")
  model<-AddForce(model,"Sleep")
  
  # A sin-wave force is applied with a period of 24h
  model<-AddSinF(model,FixPer = 24)
  
  # Start is set at intercept with speed of 0
  model<-SetYinitMode(model,mode = "Intercept_0")
  
  # Penalize the fitting for unstable value for 10 replicated days
  model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10,weight = 10)
  
  # Compute the fit using RSS
  model<-SetFittingValue(model,value = "RSS")
  
  model<-SetParametersModel(model)
  
  return(model)
  
   
}


# Full Model limits
Human_TimeCourse_FullModelBoundaries<-function(model_desync,model_ext,model_res){
  
  MeanExpr<-mean(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  SDExpr<-sd(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WiSd_upper<-c(Wake=10,Sleep=0,loggamma=log(100),omega=2*pi/9,AmpSin=10,PhiSin=2*pi,intercept = MeanExpr + 2*SDExpr)
  WiSd_lower<-c(Wake=0,Sleep=-10,loggamma=log(0.01),omega=0,AmpSin=-10,PhiSin=0,intercept = MeanExpr - 2*SDExpr)
  
  WiSd_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr)
  # Some additional initial parameters: SW driven
  WiSd_initparams1<-c(Wake=0.01,Sleep=-0.01,loggamma=log(0.1),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr)
  WiSd_initparams1<-rbind(WiSd_initparams1,
                          c(Wake=0.05,Sleep=-0.05,loggamma=log(2),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr))
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WdSi_upper<-c(Wake=0,Sleep=10,loggamma=log(100),omega=2*pi/9,AmpSin=10,PhiSin=2*pi,intercept = MeanExpr + 2*SDExpr)
  WdSi_lower<-c(Wake=-10,Sleep=0,loggamma=log(0.01),omega=0,AmpSin=-10,PhiSin=0,intercept = MeanExpr - 2*SDExpr)
  
  WdSi_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr)
  # Some additional initial parameters: SW driven
  WdSi_initparams1<-c(Wake=-0.01,Sleep=0.01,loggamma=log(0.1),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr)
  WdSi_initparams1<-rbind(WdSi_initparams1,
                         c(Wake=-0.05,Sleep=0.05,loggamma=log(2),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr))

  return(list(WiSd = list(upper = WiSd_upper, lower = WiSd_lower, params = WiSd_initparams0, params_multi = WiSd_initparams1 ),
              WdSi = list( upper = WdSi_upper, lower = WdSi_lower, params = WdSi_initparams0, params_multi = WdSi_initparams1  ) ))
}


# Some stats of fitting 
GetStatsFull<-function(resdesync,resextres,model_desync,model_ext,model_res,ProbeID){
  
  k <-7 + 1 # 7 free parameters + var
  
  # Desync
  out_desync<-SWDMrFit(model_desync,params = resdesync)
  stats_desync<-SWDMrStats(model_desync,out_desync,detailed = T)
  
  # Ext Rest sleep
  out_rest<-SWDMrFit(model_res,params = resextres)
  stats_rest<-SWDMrStats(model_res,out_rest,detailed = T)
  out_ext<-SWDMrFit(model_ext,params = resextres)
  stats_ext<-SWDMrStats(model_ext,out_ext,detailed = T)
  
  # Stats model
  ntot<-stats_desync$stats$n+stats_ext$stats$n+stats_rest$stats$n
  RSStot<-stats_desync$stats$RSS+stats_ext$stats$RSS+stats_rest$stats$RSS
  var<-RSStot/ntot
  NLL<- (ntot/2)*(log(2*pi)+log(var)+1)
  BIC <- -2*(-NLL)+(k)*log(ntot)
  AIC <- 2*(k) - 2*(-NLL)
  
  #BIC<-ntot*log(RSStot/ntot)+k*log(ntot)
  #AIC<- 2*k+ntot*log(RSStot/ntot)
  
  # flat
  xflat<-c(model_desync@Gexp[,ProbeID],model_ext@Gexp[,ProbeID],model_res@Gexp[,ProbeID])
  xexp<-c(rep("Desync",length(model_desync@Gexp[,ProbeID])),
          rep("ExtRest",length(model_ext@Gexp[,ProbeID])),
          rep("ExtRest",length(model_res@Gexp[,ProbeID])))
  RSSflat<-sum(resid(lm(xflat~xexp))^2)
  varflat<-RSSflat/ntot
  NLLflat<- (ntot/2)*(log(2*pi)+log(varflat)+1)
  BICflat <- -2*(-NLLflat)+(2+1)*log(ntot)
  
  #BICflat<-ntot*log(RSSflat/ntot)+(2+1)*log(ntot)
  
  BF<-exp((BICflat - BIC)/2)
  
  # Kendall tau
  GeneExp<-na.omit(c(stats_desync$fitted,stats_rest$fitted,stats_ext$fitted))
  predval<-na.omit(c(model_desync@Gexp[,model_desync@VarExp ], model_res@Gexp[,model_res@VarExp ], model_ext@Gexp[,model_ext@VarExp ] ))
  KendallT<-cor(GeneExp,predval,method="kendall")
  
  # likelihood model
  #var<-RSStot/ntot
  #NLL<- -1* sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
  
  return(list(RSS_model= RSStot, AIC_model= AIC, BIC_model=BIC,
              RSS_flat =RSSflat, BIC_flat = BICflat, BF =  BF,KendallTau=KendallT,NLL = NLL,
              resid_desync=stats_desync$residuals,fitted_desync=stats_desync$fitted,
              resid_rest=stats_rest$residuals,fitted_rest=stats_rest$fitted,
              resid_ext=stats_ext$residuals,fitted_ext=stats_ext$fitted))
}

# FITTING 
FitDataFull<-function(ProbeID, MeanSWdf, HT_Desync, MeanSWdfExt, HT_Ext, MeanSWdfRest, HT_Res ){
  
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
 
  fits<-optimx(paramsL$WiSd$params,fn = objFAll,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
  for (i in 1:nrow(paramsL$WiSd$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WiSd$params_multi[i,],fn = objFAll,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper) )
  }
  
  fits<-rbind(fits, optimx(paramsL$WdSi$params,fn = objFAll,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  for (i in 1:nrow(paramsL$WdSi$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WdSi$params_multi[i,],fn = objFAll,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  }  
  
  
  optimxres<-fits[order(fits$value),,drop=F][1,]
  
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

###################################################
########### CIRCADIAN MODEL #########

Human_TimeCourse_GetCircadianModel<-function(SWdf,Exprdf,ProbeID,intercept=NULL){
  
  swdmr<-SWDMr(SWdist=SWdf, Gexp=Exprdf[,c(ProbeID,"Time")])
  
  ############# Build model ###############
  model<-initDDHOmodel(swdmr, VarExp = ProbeID, UseDampingRatio = F)
  
  # Set an intercept if given
  if (! is.null(intercept)){
    model<-FixIntercept(model,intercept)
  }
  
  
  # We replicate baseline for some days
  model<-ReplicateDrivingForce(model,c(12.0,36.0),20)
  
  # A sin-wave force is applied with a period of 24h
  model<-AddSinF(model,FixPer = 24)
  
  # Start is set at intercept with speed of 0
  model<-SetYinitMode(model,mode = "Intercept_0")
  
  # Penalize the fitting for unstable value for 10 replicated days
  model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10,weight = 10)
  
  # Compute the fit using RSS
  model<-SetFittingValue(model,value = "RSS")
  
  model<-SetParametersModel(model)
  
  return(model)
  
  
}

# circadian Model limits
Human_TimeCourse_CircadianModelBoundaries<-function(model_desync,model_ext,model_res){
  
  MeanExpr<-mean(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  SDExpr<-sd(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  
  F_upper<-c(loggamma=log(100),omega=2*pi/9,AmpSin=10,PhiSin=2*pi,intercept = MeanExpr + 2*SDExpr)
  F_lower<-c(loggamma=log(0.01),omega=0,AmpSin=-10,PhiSin=0,intercept = MeanExpr - 2*SDExpr)
  F_initparams0<-c(loggamma=log(1e-1),omega=2*pi/24,AmpSin=0,PhiSin=pi,intercept = MeanExpr)
  
  F_initparams1<-c(loggamma=log(1e-1),omega=2*pi/24,AmpSin=0.01,PhiSin=pi,intercept = MeanExpr)
  F_initparams1<-rbind(F_initparams1,
                        c(loggamma=log(1e-1),omega=2*pi/24,AmpSin=-0.01,PhiSin=pi,intercept = MeanExpr))  
  F_initparams1<-rbind(F_initparams1,
                       c(loggamma=log(2),omega=2*pi/24,AmpSin=0.05,PhiSin=pi,intercept = MeanExpr))  
  F_initparams1<-rbind(F_initparams1,
                       c(loggamma=log(2),omega=2*pi/24,AmpSin=-0.05,PhiSin=pi,intercept = MeanExpr))  
  

  return(list(upper = F_upper, lower = F_lower, params = F_initparams0, params_multi = F_initparams1 ))
  
}

# Some stats of fitting 
GetStatsCircadian<-function(resdesync,resextres,model_desync,model_ext,model_res,ProbeID){
  
  k <- 5 + 1
  
  # Desync
  out_desync<-SWDMrFit(model_desync,params = resdesync)
  stats_desync<-SWDMrStats(model_desync,out_desync,detailed = T)
  
  # Ext Rest sleep
  out_rest<-SWDMrFit(model_res,params = resextres)
  stats_rest<-SWDMrStats(model_res,out_rest,detailed = T)
  out_ext<-SWDMrFit(model_ext,params = resextres)
  stats_ext<-SWDMrStats(model_ext,out_ext,detailed = T)
  
  # Stats model
  ntot<-stats_desync$stats$n+stats_ext$stats$n+stats_rest$stats$n
  RSStot<-stats_desync$stats$RSS+stats_ext$stats$RSS+stats_rest$stats$RSS
  var<-RSStot/ntot
  NLL<- (ntot/2)*(log(2*pi)+log(var)+1)
  BIC <- -2*(-NLL)+(k)*log(ntot)
  AIC <- 2*(k) - 2*(-NLL)
  # BIC<-ntot*log(RSStot/ntot)+k*log(ntot)
  # AIC<- 2*k+ntot*log(RSStot/ntot)
  
  # flat
  xflat<-c(model_desync@Gexp[,ProbeID],model_ext@Gexp[,ProbeID],model_res@Gexp[,ProbeID])
  xexp<-c(rep("Desync",length(model_desync@Gexp[,ProbeID])),
          rep("ExtRest",length(model_ext@Gexp[,ProbeID])),
          rep("ExtRest",length(model_res@Gexp[,ProbeID])))
  RSSflat<-sum(resid(lm(xflat~xexp))^2)
  #BICflat<-ntot*log(RSSflat/ntot)+(2+1)*log(ntot)
  varflat<-RSSflat/ntot
  NLLflat<- (ntot/2)*(log(2*pi)+log(varflat)+1)
  BICflat <- -2*(-NLLflat)+(2+1)*log(ntot)
  
  BF<-exp((BICflat - BIC)/2)
  
  
  # Kendall tau
  GeneExp<-na.omit(c(stats_desync$fitted,stats_rest$fitted,stats_ext$fitted))
  predval<-na.omit(c(model_desync@Gexp[,model_desync@VarExp ], model_res@Gexp[,model_res@VarExp ], model_ext@Gexp[,model_ext@VarExp ] ))
  KendallT<-cor(GeneExp,predval,method="kendall")
  
  # likelihood model
  var<-RSStot/ntot
  NLL<- -1* sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
  
  return(list(RSS_model= RSStot, AIC_model= AIC, BIC_model=BIC,
              RSS_flat =RSSflat, BIC_flat = BICflat, BF =  BF, KendallTau=KendallT,NLL = NLL))
}


# FITTING 
FitDataCircadian<-function(ProbeID, MeanSWdf, HT_Desync, MeanSWdfExt, HT_Ext, MeanSWdfRest, HT_Res ){
  
  # Get intercept desynchronization protocols
  Expr<-HT_Desync[HT_Desync$Time < 72,ProbeID]
  Time<-HT_Desync[HT_Desync$Time < 72,"Time"]
  MeanGeneExprInBaseline<-coef(lm(Expr~sin(2*pi/24*Time) + cos(2*pi/24*Time) ))[["(Intercept)"]]
  
  # Build models
  model_desync<-Human_TimeCourse_GetCircadianModel(MeanSWdf,HT_Desync,ProbeID,intercept = MeanGeneExprInBaseline)
  model_ext<-Human_TimeCourse_GetCircadianModel(MeanSWdfExt,HT_Ext,ProbeID)
  model_res<-Human_TimeCourse_GetCircadianModel(MeanSWdfRest,HT_Res,ProbeID)
  
  # Prepare objective function
  objFdesync<-SWDMrGetEvalFun(model_desync)
  objFext<-SWDMrGetEvalFun(model_ext)
  objFres<-SWDMrGetEvalFun(model_res)
  
  # Prepare models limits
  paramsL<-Human_TimeCourse_CircadianModelBoundaries(model_desync,model_ext,model_res)
  
  # Combine RSS for all fits
  objFAll<-function(params){
    paramsdesync<-params[c(1,2,3,4)]
    paramsextres<-params[c(1,2,3,4,5)]
    
    RSS_desync<-objFdesync(paramsdesync)
    RSS_ext<-objFext(paramsextres)
    RSS_res<-objFres(paramsextres)
    
    RSSAll<-RSS_desync+RSS_ext+RSS_res
    
    return(RSSAll)
    
  }
  
  # Under damping fits
  fits<-optimx(paramsL$params,fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$lower,upper=paramsL$upper)
  
  for (i in 1:nrow(paramsL$params_multi)){
    fits<-rbind(fits, optimx(paramsL$params_multi[i,],fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$lower,upper=paramsL$upper) )
  }
  
  optimxres<-fits[order(fits$value),,drop=F][1,]
  

  # Split between desync and rest-ext sleep (intercept change)
  resdesync<-optimxres[,c(1,2,3,4)]
  resextres<-optimxres[,c(1,2,3,4,5)]
  
  # Get Statistics
  stats<-GetStatsCircadian(resdesync,resextres,model_desync,model_ext,model_res,ProbeID)
  
  return(list(paramsDesync=resdesync,paramsExtRes=resextres,stats=stats,
              model_desync=model_desync,model_ext=model_ext, model_res = model_res))
}


###################################################
########### Drop SinF model #########

Human_TimeCourse_DropSinFModel<-function(SWdf,Exprdf,ProbeID,intercept=NULL){
  
  swdmr<-SWDMr(SWdist=SWdf, Gexp=Exprdf[,c(ProbeID,"Time")])
  
  ############# Build model ###############
  model<-initDDHOmodel(swdmr, VarExp = ProbeID, UseDampingRatio = F)
  
  # Set an intercept if given
  if (! is.null(intercept)){
    model<-FixIntercept(model,intercept)
  }
  
  # We replicate baseline for some days
  model<-ReplicateDrivingForce(model,c(12.0,36.0),20)
  
  # Add sleep-wake force
  model<-AddForce(model,"Wake")
  model<-AddForce(model,"Sleep")
  
  # Start is set at intercept with speed of 0
  model<-SetYinitMode(model,mode = "Intercept_0")
  
  # Penalize the fitting for unstable value for 10 replicated days
  model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10,weight = 10)
  
  # Compute the fit using RSS
  model<-SetFittingValue(model,value = "RSS")
  
  model<-SetParametersModel(model)
  
  return(model)
  
}


# Full Model limits
Human_TimeCourse_DropSinFModelBoundaries<-function(model_desync,model_ext,model_res){
  
  MeanExpr<-mean(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  SDExpr<-sd(c(model_ext@Gexp[,model_ext@VarExp],model_res@Gexp[,model_res@VarExp]),na.rm=T)
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WiSd_upper<-c(Wake=10,Sleep=0,loggamma=log(100),omega=2*pi/9,intercept = MeanExpr + 2*SDExpr)
  WiSd_lower<-c(Wake=0,Sleep=-10,loggamma=log(0.01),omega=0,intercept = MeanExpr - 2*SDExpr)
  
  WiSd_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24,intercept = MeanExpr)
  # Some additional initial parameters: SW driven
  WiSd_initparams1<-c(Wake=0.01,Sleep=-0.01,loggamma=log(0.1),omega=2*pi/24,intercept = MeanExpr)
  WiSd_initparams1<-rbind(WiSd_initparams1,
                          c(Wake=0.05,Sleep=-0.05,loggamma=log(2),omega=2*pi/24,intercept = MeanExpr))
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WdSi_upper<-c(Wake=0,Sleep=10,loggamma=log(100),omega=2*pi/9,intercept = MeanExpr + 2*SDExpr)
  WdSi_lower<-c(Wake=-10,Sleep=0,loggamma=log(0.01),omega=0,intercept = MeanExpr - 2*SDExpr)
  
  WdSi_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24,intercept = MeanExpr)
  # Some additional initial parameters: SW driven
  WdSi_initparams1<-c(Wake=-0.01,Sleep=0.01,loggamma=log(0.1),omega=2*pi/24,intercept = MeanExpr)
  WdSi_initparams1<-rbind(WdSi_initparams1,
                          c(Wake=-0.05,Sleep=0.05,loggamma=log(2),omega=2*pi/24,intercept = MeanExpr))
  
  return(list(WiSd = list(upper = WiSd_upper, lower = WiSd_lower, params = WiSd_initparams0, params_multi = WiSd_initparams1 ),
              WdSi = list( upper = WdSi_upper, lower = WdSi_lower, params = WdSi_initparams0, params_multi = WdSi_initparams1  ) ))
}


# Some stats of fitting 
GetStatsDropSinF<-function(resdesync,resextres,model_desync,model_ext,model_res,ProbeID){
  
  k <- 5 + 1
  
  # Desync
  out_desync<-SWDMrFit(model_desync,params = resdesync)
  stats_desync<-SWDMrStats(model_desync,out_desync,detailed = T)
  
  # Ext Rest sleep
  out_rest<-SWDMrFit(model_res,params = resextres)
  stats_rest<-SWDMrStats(model_res,out_rest,detailed = T)
  out_ext<-SWDMrFit(model_ext,params = resextres)
  stats_ext<-SWDMrStats(model_ext,out_ext,detailed = T)
  
  # Stats model
  ntot<-stats_desync$stats$n+stats_ext$stats$n+stats_rest$stats$n
  RSStot<-stats_desync$stats$RSS+stats_ext$stats$RSS+stats_rest$stats$RSS
  var<-RSStot/ntot
  NLL<- (ntot/2)*(log(2*pi)+log(var)+1)
  BIC <- -2*(-NLL)+(k)*log(ntot)
  AIC <- 2*(k) - 2*(-NLL)
  # BIC<-ntot*log(RSStot/ntot)+k*log(ntot)
  # AIC<- 2*k+ntot*log(RSStot/ntot)
  
  # flat
  xflat<-c(model_desync@Gexp[,ProbeID],model_ext@Gexp[,ProbeID],model_res@Gexp[,ProbeID])
  xexp<-c(rep("Desync",length(model_desync@Gexp[,ProbeID])),
          rep("ExtRest",length(model_ext@Gexp[,ProbeID])),
          rep("ExtRest",length(model_res@Gexp[,ProbeID])))
  RSSflat<-sum(resid(lm(xflat~xexp))^2)
  #BICflat<-ntot*log(RSSflat/ntot)+(2+1)*log(ntot)
  varflat<-RSSflat/ntot
  NLLflat<- (ntot/2)*(log(2*pi)+log(varflat)+1)
  BICflat <- -2*(-NLLflat)+(2+1)*log(ntot)
  
  BF<-exp((BICflat - BIC)/2)
  
  # Kendall tau
  GeneExp<-na.omit(c(stats_desync$fitted,stats_rest$fitted,stats_ext$fitted))
  predval<-na.omit(c(model_desync@Gexp[,model_desync@VarExp ], model_res@Gexp[,model_res@VarExp ], model_ext@Gexp[,model_ext@VarExp ] ))
  KendallT<-cor(GeneExp,predval,method="kendall")
  
  # likelihood model
  var<-RSStot/ntot
  NLL<- -1* sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
  
  return(list(RSS_model= RSStot, AIC_model= AIC, BIC_model=BIC,
              RSS_flat =RSSflat, BIC_flat = BICflat, BF =  BF,KendallTau=KendallT,NLL = NLL))
}

# FITTING 
FitDataDropSinF<-function(ProbeID, MeanSWdf, HT_Desync, MeanSWdfExt, HT_Ext, MeanSWdfRest, HT_Res ){
  
  # Get intercept desynchronization protocols
  Expr<-HT_Desync[HT_Desync$Time < 72,ProbeID]
  Time<-HT_Desync[HT_Desync$Time < 72,"Time"]
  MeanGeneExprInBaseline<-coef(lm(Expr~sin(2*pi/24*Time) + cos(2*pi/24*Time) ))[["(Intercept)"]]
  
  # Build models
  model_desync<-Human_TimeCourse_DropSinFModel(MeanSWdf,HT_Desync,ProbeID,intercept=MeanGeneExprInBaseline)
  model_ext<-Human_TimeCourse_DropSinFModel(MeanSWdfExt,HT_Ext,ProbeID)
  model_res<-Human_TimeCourse_DropSinFModel(MeanSWdfRest,HT_Res,ProbeID)
  
  # Prepare objective function
  objFdesync<-SWDMrGetEvalFun(model_desync)
  objFext<-SWDMrGetEvalFun(model_ext)
  objFres<-SWDMrGetEvalFun(model_res)
  
  # Prepare models limits
  paramsL<-Human_TimeCourse_DropSinFModelBoundaries(model_desync,model_ext,model_res)
  
  # Combine RSS for all fits
  objFAll<-function(params){
    paramsdesync<-params[c(1,2,3,4)]
    paramsextres<-params[c(1,2,3,4,5)]
    
    RSS_desync<-objFdesync(paramsdesync)
    RSS_ext<-objFext(paramsextres)
    RSS_res<-objFres(paramsextres)
    
    RSSAll<-RSS_desync+RSS_ext+RSS_res
    
    return(RSSAll)
    
  }
  
  fits<-optimx(paramsL$WiSd$params,fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
  for (i in 1:nrow(paramsL$WiSd$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WiSd$params_multi[i,],fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper) )
  }
  
  fits<-rbind(fits, optimx(paramsL$WdSi$params,fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  for (i in 1:nrow(paramsL$WdSi$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WdSi$params_multi[i,],fn = objFAll,method=c("nlminb"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  }  
  
  
  optimxres<-fits[order(fits$value),,drop=F][1,]
  
  # Split between desync and rest-ext sleep (intercept change)
  resdesync<-optimxres[,c(1,2,3,4)]
  resextres<-optimxres[,c(1,2,3,4,5)]
  
  # Get Statistics
  stats<-GetStatsDropSinF(resdesync,resextres,model_desync,model_ext,model_res,ProbeID)
  
  return(list(paramsDesync=resdesync,paramsExtRes=resextres,stats=stats,
              model_desync=model_desync,model_ext=model_ext, model_res = model_res))
}