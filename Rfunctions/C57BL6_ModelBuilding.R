###################################################

########### FULL MODEL #########
# This model is driven by Sleep-Wake and a second force that is circadian (i.e. undisturbed by sleep deprivation)

C57BL6_TimeCourse_GetFullModel<-function(swdmr,Gene){
  
  ############# Build model ###############
  model<-initDDHOmodel(swdmr, VarExp = Gene, UseDampingRatio = F)
  
  # Fix the intercept using sin-wave fitting in baseline
  Time<-swdmr@Gexp[,"Time"]
  Expr<-swdmr@Gexp[,Gene]
  MeanGeneExprInBaseline<-coef(lm(Expr[Time<=48]~sin(2*pi/24*Time[Time<=48]) + cos(2*pi/24*Time[Time<=48]) ))[["(Intercept)"]]
  model<-FixIntercept(model,MeanGeneExprInBaseline)
  
  # Add sleep-wake force
  model<-AddForce(model,"Wake")
  model<-AddForce(model,"Sleep")
  
  # A sin-wave force is applied with a period of 24h
  model<-AddSinF(model,FixPer = 24)
  
  # Start is set at intercept with speed of 0
  model<-SetYinitMode(model,mode = "Intercept_0",values = c(0,48))
  
  # We replicate the first 24h baseline for 20 day
  model<-ReplicateDrivingForce(model,c(0.1,24.0),25)
  
  # Penalize the fitting for unstable value for 10 replicated days
  model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10,weight = 100)
  
  # Compute the fit using RSS
  model<-SetFittingValue(model,value = "RSS")
  
  # Set parameters
  model<-SetParametersModel(model)
  
  return(model)
}


##### Full Model limits ######
# Upper, Lower and initial values of the full model
# Limits on omega set to avoid very high natural frequency of the oscillator (> RNA-seq sampling frequency)
# Limits on Wake, Sleep and AmpSin set to avoid very large change, very fast. RK4 timestep is 6min (0.1h)
C57BL6_TimeCourse_FullModelBoundaries<-function(){
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WiSd_upper<-c(Wake=10,Sleep=0,loggamma=log(100),omega=2*pi/12,AmpSin=10,PhiSin=2*pi)
  WiSd_lower<-c(Wake=0,Sleep=-10,loggamma=log(0.01),omega=2*pi/72,AmpSin=-10,PhiSin=0)
  WiSd_initparams0<-c(Wake=0,Sleep=0,loggamma=log(1e-1),omega=2*pi/24,AmpSin=0,PhiSin=pi)
  
  # Some additional initial parameters: SW driven
  WiSd_initparams1<-c(Wake=0.01,Sleep=-0.01,loggamma=1e-1,omega=2*pi/24,AmpSin=0,PhiSin=pi)
  WiSd_initparams1<-rbind(WiSd_initparams1,
                          c(Wake=0.05,Sleep=-0.05,loggamma=log(2),omega=2*pi/24,AmpSin=0,PhiSin=pi))
  
  ######  Wake decrease mRNA speed, Sleep increase mRNA speed  
  WdSi_upper<-c(Wake=0,Sleep=10,loggamma=log(100),omega=2*pi/12,AmpSin=10,PhiSin=2*pi)
  WdSi_lower<-c(Wake=-10,Sleep=0,loggamma=log(0.01),omega=2*pi/72,AmpSin=-10,PhiSin=0)
  WdSi_initparams0<-c(Wake=0,Sleep=0,loggamma=log(1e-1),omega=2*pi/24,AmpSin=0,PhiSin=pi)
  
  # Some additional initial parameters: SW driven
  WdSi_initparams1<-c(Wake=-0.01,Sleep=0.01,loggamma=log(1e-1),omega=2*pi/24,AmpSin=0,PhiSin=pi)
  WdSi_initparams1<-rbind(WdSi_initparams1,
                          c(Wake=-0.05,Sleep=0.05,loggamma=log(2),omega=2*pi/24,AmpSin=0,PhiSin=pi))
  
  
  return(list(WiSd = list(upper = WiSd_upper, lower = WiSd_lower, params = WiSd_initparams0, params_multi = WiSd_initparams1 ),
              WdSi = list(upper = WdSi_upper, lower = WdSi_lower, params = WdSi_initparams0, params_multi = WdSi_initparams1 ) ))
  
}



###################################################
###################################################

########### CIRCADIAN MODEL ###########
# This model is driven only by a force that is circadian (i.e. undisturbed by sleep deprivation)
C57BL6_TimeCourse_GetCircadianModel<-function(swdmr,Gene){
  
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  
  # Sleep and Wake have a force of 0
  model<-AddForce(model,"Sleep",0)
  model<-AddForce(model,"Wake",0)
  
  model<-SetParametersModel(model)
  
  return(model)
}


##### Circadian Model limits ######
# Upper, Lower and initial values of the circadian model
C57BL6_TimeCourse_CircadianModelBoundaries<-function(){
  
  ######  Underdamping condition 
  upper<-c(loggamma=log(100),omega=2*pi/12,AmpSin=10,PhiSin=2*pi)
  lower<-c(loggamma=log(0.01),omega=2*pi/72,AmpSin=-10,PhiSin=0)
  initparams0<-c(loggamma=log(0.1),omega=2*pi/24,AmpSin=0,PhiSin=pi)
  
  # Some additional initial parameters
  initparams1<-c(loggamma=log(0.1),omega=2*pi/24,AmpSin=0.01,PhiSin=pi)
  initparams1<-rbind(initparams1,
                     c(loggamma=log(0.1),omega=2*pi/24,AmpSin=-0.01,PhiSin=pi))
  initparams1<-rbind(initparams1,
                     c(loggamma=log(2),omega=2*pi/24,AmpSin=0.01,PhiSin=pi))
  initparams1<-rbind(initparams1,
                     c(loggamma=log(2),omega=2*pi/24,AmpSin=-0.01,PhiSin=pi))
  
  return(list(upper = upper, lower = lower, params = initparams0, params_multi = initparams1 ))
  
}


########### NOSinWave MODEL #########
# This model is driven by Sleep-Wake

C57BL6_TimeCourse_GetNoSinWaveModel<-function(swdmr,Gene){
  
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  
  # SinF have an amplitude of 0
  model<-AddSinF(model,FixPer = 24,FixPhi = 0,FixAmp = 0)
  
  model<-SetParametersModel(model)
  
  return(model)
}


##### Full Model limits ######
# Upper, Lower and initial values of the full model
# Limits on omega set to avoid very high natural frequency of the oscillator (> RNA-seq sampling frequency)
# Limits on Wake, Sleep and AmpSin set to avoid very large change, very fast. RK4 timestep is 6min (0.1h)
C57BL6_TimeCourse_NoSinWaveModelBoundaries<-function(){
  
  ######  Wake increase mRNA speed, Sleep decrease mRNA speed 
  WiSd_upper<-c(Wake=10,Sleep=0,loggamma=log(100),omega=2*pi/12)
  WiSd_lower<-c(Wake=0,Sleep=-10,loggamma=log(0.01),omega=2*pi/72)
  WiSd_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24)
  
  # Some additional initial parameters: SW driven
  WiSd_initparams1<-c(Wake=0.01,Sleep=-0.01,loggamma=log(0.1),omega=2*pi/24)
  WiSd_initparams1<-rbind(WiSd_initparams1,
                          c(Wake=0.05,Sleep=-0.05,loggamma=log(2),omega=2*pi/24))
  WiSd_initparams1<-rbind(WiSd_initparams1,
                          c(Wake=.1,Sleep=-.1,loggamma=log(5),omega=2*pi/12))
  
  ######  Wake decrease mRNA speed, Sleep increase mRNA speed  
  WdSi_upper<-c(Wake=0,Sleep=10,loggamma=log(100),omega=2*pi/12)
  WdSi_lower<-c(Wake=-10,Sleep=0,loggamma=log(0.01),omega=2*pi/72)
  WdSi_initparams0<-c(Wake=0,Sleep=0,loggamma=log(0.1),omega=2*pi/24)
  
  # Some additional initial parameters: SW driven
  WdSi_initparams1<-c(Wake=-0.01,Sleep=0.01,loggamma=log(0.1),omega=2*pi/24)
  WdSi_initparams1<-rbind(WdSi_initparams1,
                          c(Wake=-0.05,Sleep=0.05,loggamma=log(2),omega=2*pi/24))
  WdSi_initparams1<-rbind(WdSi_initparams1,
                          c(Wake=-.1,Sleep=.1,loggamma=log(5),omega=2*pi/12))
  
  return(list(WiSd = list(upper = WiSd_upper, lower = WiSd_lower, params = WiSd_initparams0, params_multi = WiSd_initparams1 ),
              WdSi = list(upper = WdSi_upper, lower = WdSi_lower, params = WdSi_initparams0, params_multi = WdSi_initparams1 ) ))
  
}



FitFullModel<-function(swdmr,Gene){
  
  model<-C57BL6_TimeCourse_GetFullModel(swdmr,Gene)
  objfun<-SWDMrGetEvalFun(model)
  paramsL<-C57BL6_TimeCourse_FullModelBoundaries()
  
  # Wake increase mRNA speed, Sleep decrease mRNA speed 
  fits<-optimx(paramsL$WiSd$params,fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
  for (i in 1:nrow(paramsL$WiSd$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WiSd$params_multi[i,],fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper) )
  }
  
  # Wake decrease mRNA speed, Sleep increase mRNA speed
  fits<-rbind(fits, optimx(paramsL$WdSi$params,fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  for (i in 1:nrow(paramsL$WdSi$params_multi)){
    fits<-rbind(fits, optimx(paramsL$WdSi$params_multi[i,],fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  }  
  
  optimxres<-fits[order(fits$value),,drop=F][1,]
  out<-SWDMrFit(model,params = optimxres)
  stats<-SWDMrStats(model,out,detailed = T)
  
  return(list(optimxres = optimxres,stats = stats))
}



FitDropSW<-function(swdmr,Gene){
  
  modelC<-C57BL6_TimeCourse_GetCircadianModel(swdmr,Gene)
  objfun<-SWDMrGetEvalFun(modelC)
  
  paramsL<-C57BL6_TimeCourse_CircadianModelBoundaries()
  
  fitsC<-optimx(paramsL$params,fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$lower,upper=paramsL$upper)
  for (i in 1:nrow(paramsL$params_multi)){
    fitsC<-rbind(fitsC, optimx(paramsL$params_multi[i,],fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$lower,upper=paramsL$upper) )
  }
  
  optimxresC<-fitsC[order(fitsC$value),,drop=F][1,]
  outC<-SWDMrFit(modelC,params = optimxresC)
  statsC<-SWDMrStats(modelC,outC,detailed = T)
  
  return(list(optimxres = optimxresC,stats = statsC))
}


FitDropSinF<-function(swdmr,Gene){
  
  modelSWo<-C57BL6_TimeCourse_GetNoSinWaveModel(swdmr,Gene)
  objfun<-SWDMrGetEvalFun(modelSWo)
  paramsL<-C57BL6_TimeCourse_NoSinWaveModelBoundaries()
  
  # Wake increase mRNA speed, Sleep decrease mRNA speed 
  fitsSWo<-optimx(paramsL$WiSd$params,fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper)
  for (i in 1:nrow(paramsL$WiSd$params_multi)){
    fitsSWo<-rbind(fitsSWo, optimx(paramsL$WiSd$params_multi[i,],fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WiSd$lower,upper=paramsL$WiSd$upper) )
  }
  
  # Wake decrease mRNA speed, Sleep increase mRNA speed
  fitsSWo<-rbind(fitsSWo, optimx(paramsL$WdSi$params,fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  for (i in 1:nrow(paramsL$WdSi$params_multi)){
    fitsSWo<-rbind(fitsSWo, optimx(paramsL$WdSi$params_multi[i,],fn = objfun,method=c("nlminb","L-BFGS-B"),control=list(maxit=1000),lower=paramsL$WdSi$lower,upper=paramsL$WdSi$upper) )
  }  
  
  optimxresSWo<-fitsSWo[order(fitsSWo$value),,drop=F][1,]
  outCSWo<-SWDMrFit(modelSWo,params = optimxresSWo)
  statsSWo<-SWDMrStats(modelSWo,outCSWo,detailed = T)
  
  return(list(optimxres = optimxresSWo,stats = statsSWo))
}
