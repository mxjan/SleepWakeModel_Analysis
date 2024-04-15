# Calculate Sleep-Wake response contribution in mouse dataset
GetSWrcMouse<-function(Genes,swdmr,fits){
  SWrcs<-c()
  for (Gene in Genes){
    model<-C57BL6_TimeCourse_GetFullModel(swdmr =swdmr ,Gene=Gene)
    out<-SWDMrFit(model,fits[fits$Gene %in% Gene,][1,],method="Solve")
    # Baseline index
    idx<-out$time < 24 & out$time >0
    # Circadian amplitude (peak to trough)
    CC<-max(out$circ_sol[idx])-min(out$circ_sol[idx])
    # Sleep-wake amplitude (peak to trough)
    SS<-max(out$SW_response[idx]+out$trans_sol[idx])-min(out$SW_response[idx]+out$trans_sol[idx])
    # SWrc
    SWrc<-SS/(SS+CC)
    SWrcs<-c(SWrcs,SWrc)
    names(SWrcs)[length(SWrcs)]<-Gene
  }
  return(SWrcs)
}


# Calculate Sleep-Wake response contribution in human dataset
GetSWrcHuman<-function(ProbeIDs,SWdf_human,HT_Desync,fits,AfterCR=F,humanmodels_fun="Human_ModelBuilding.R"){
  if (exists("HUMANMODELS_FUN")){
    source(HUMANMODELS_FUN)
  }else{
    source(humanmodels_fun)
  }
  
  SWrcs<-c()
  for (ProbeID in ProbeIDs){
    # Get intercept desynchronization protocols
    Expr<-HT_Desync[HT_Desync$Time < 72,ProbeID]
    Time<-HT_Desync[HT_Desync$Time < 72,"Time"]
    MeanGeneExprInBaseline<-coef(lm(Expr~sin(2*pi/24*Time) + cos(2*pi/24*Time) ))[["(Intercept)"]]
    
    if (AfterCR == F){
      model<-Human_TimeCourse_GetFullModel(SWdf = SWdf_human,Exprdf = HT_Desync,ProbeID = ProbeID,intercept = MeanGeneExprInBaseline )
      out<-SWDMrFit(model,fits[fits$ProbeName %in% ProbeID,][1,],method="Solve")
      
      if (! any(is.na(out$trans_sol))){
        idx<-out$time < 24 & out$time >0
        
        CC<-max(out$circ_sol[idx])-min(out$circ_sol[idx])
        SS<-max(out$SW_response[idx]+out$trans_sol[idx])-min(out$SW_response[idx]+out$trans_sol[idx])
        
        SWrc<-SS/(SS+CC)
      }else{
        SWrc<-NA
      }
    }else{
      
      model<-Human_TimeCourse_GetFullModel(SWdf = SWdf_human,Exprdf = HT_Desync,ProbeID = ProbeID)
      out<-SWDMrFit(model,fits[fits$ProbeName %in% ProbeID,][1,],method="Solve")
      if (! any(is.na(out$trans_sol))){
        idx<-out$time > 192 & out$time < 258
        
        CC<-max(out$circ_sol[idx])-min(out$circ_sol[idx])
        SS<-max(out$SW_response[idx]+out$trans_sol[idx])-min(out$SW_response[idx]+out$trans_sol[idx])
        
        SWrc<-SS/(SS+CC)
      }else{
        SWrc<-NA
      }
    }
    

    
    SWrcs<-c(SWrcs,SWrc)
    names(SWrcs)[length(SWrcs)]<-ProbeID
  }

  return(SWrcs)
}


# Get time constant to get back to equilibrium position or baseline oscillation
Gettau<-function(Genes,fits){
  taus<-c()
  # Human dataset
  for (Gene in Genes){
    if ("ProbeName" %in% colnames(fits)){
    
      idx<-which(fits$ProbeName %in% Gene)[1]
    }else{
      idx<-which(fits$Gene %in% Gene)[1]
    }
    omega<-fits[idx,"omega"]
    gamma<-exp(fits[idx,"loggamma"])
    
    if (gamma/(2*omega)  > 1){
      tau<- -1/ (-gamma/2 + sqrt( (gamma/2)^2-omega^2))
    }else{
      tau <- 2/gamma
    }
    taus<-c(taus,tau)
    names(taus)[length(taus)]<-Gene
  }
  return(taus)
}

GetZeta<-function(Genes,fits){
  zetas<-c()
  for (Gene in Genes){
    if ("ProbeName" %in% colnames(fits)){
      idx<-which(fits$ProbeName %in% Gene)[1]
    }else{
      idx<-which(fits$Gene %in% Gene)[1]
    }
    omega<-fits[idx,"omega"]
    gamma<-exp(fits[idx,"loggamma"])
    
    zetas<-c(zetas,gamma/(2*omega))
    names(zetas)[length(zetas)]<-Gene
  }
  return(zetas)
}

GetUpperLowerEquilibrium<-function(Genes,fits,Exprdata,asSD=F){
  
  uppers<-c()
  lowers<-c()
  # Human dataset
  for (Gene in Genes){
    if ("ProbeName" %in% colnames(fits)){
      idx<-which(fits$ProbeName %in% Gene)[1]
    }else{
      idx<-which(fits$Gene %in% Gene)[1]
      
      # Equilibrium 
      Time<-Exprdata[,"Time"]
      Expr<-Exprdata[,Gene]
      MeanGeneExprInBaseline<-coef(lm(Expr[Time<=48]~sin(2*pi/24*Time[Time<=48]) + cos(2*pi/24*Time[Time<=48]) ))[["(Intercept)"]]
    }
    
    SleepWake<-fits[idx,c("Wake","Sleep")]
    omega<-fits[idx,"omega"]
    
    upper<-(as.numeric(SleepWake[which(SleepWake > 0)])*0.1)/omega^2
    lower<-(as.numeric(SleepWake[which(SleepWake < 0)])*0.1)/omega^2
    
    if (length(upper) == 0){upper<-0}
    if (length(lower) == 0){lower<-0}
    
    upper<-upper + MeanGeneExprInBaseline
    lower<-lower + MeanGeneExprInBaseline
    
    if (asSD == T){
      upper<-(upper-mean(Expr))/sd(Expr)
      lower<-(lower-mean(Expr))/sd(Expr)
    }
    
    
    uppers<-c(uppers,upper)
    lowers<-c(lowers,lower)
    
    names(uppers)[length(uppers)]<-Gene
    names(lowers)[length(lowers)]<-Gene
    
  }
  return(list(lower=lowers,upper=uppers))
}


# Simple phase lag if driving force is sinusoidal, answer in the form of A*cos(wt - theta)
# For ω < ω0 the sign of the position and the force are the same
# For ω > ω0, the sign of the solution flips
GetPhaseLag<-function(Genes,fits,PosAmp=F){
  thetas<-c()
  drivingfreq<-2*pi/24
  for (Gene in Genes){
    if ("ProbeName" %in% colnames(fits)){
      idx<-which(fits$ProbeName %in% Gene)[1]
    }else{
      idx<-which(fits$Gene %in% Gene)[1]
    }
    omega<-fits[idx,"omega"]
    gamma<-exp(fits[idx,"loggamma"])
    theta<-atan( (drivingfreq*gamma)/((omega)^2-drivingfreq^2) )
    
    if (PosAmp == T){
      if (drivingfreq > omega){
        theta<-theta+pi
      }
    }
    
    thetas<-c(thetas,theta)
    names(thetas)[length(thetas)]<-Gene
  }
  return(thetas)
}

GetWakeSleepEquilibrium<-function(Genes,fits,Exprdata){
  
  WakeEqs<-c()
  SleepEqs<-c()
  # Human dataset
  for (Gene in Genes){
    if ("ProbeName" %in% colnames(fits)){
      idx<-which(fits$ProbeName %in% Gene)[1]
    }else{
      idx<-which(fits$Gene %in% Gene)[1]
      
      # Equilibrium 
      Time<-Exprdata[,"Time"]
      Expr<-Exprdata[,Gene]
      MeanGeneExprInBaseline<-coef(lm(Expr[Time<=48]~sin(2*pi/24*Time[Time<=48]) + cos(2*pi/24*Time[Time<=48]) ))[["(Intercept)"]]
    }
    
    SleepWake<-fits[idx,c("Wake","Sleep")]
    omega<-fits[idx,"omega"]
    
    WakeEq<-(as.numeric(SleepWake["Wake"])*0.1)/omega^2
    SleepEq<-(as.numeric(SleepWake["Sleep"])*0.1)/omega^2
    
    if (length(WakeEq) == 0){WakeEq<-0}
    if (length(SleepEq) == 0){SleepEq<-0}
    
    WakeEq<-WakeEq + MeanGeneExprInBaseline
    SleepEq<-SleepEq + MeanGeneExprInBaseline
    
    
    WakeEqs<-c(WakeEqs,WakeEq)
    SleepEqs<-c(SleepEqs,SleepEq)
    
    names(WakeEqs)[length(WakeEqs)]<-Gene
    names(SleepEqs)[length(SleepEqs)]<-Gene
    
  }
  return(list(WakeEq=WakeEqs,SleepEq=SleepEqs))
}


# Get the p-value of a cosine model
PvalCosineModel<-function(x,Time){
  idx<-is.na(x)
  x<-x[! idx]
  Time<-Time[!idx]
  cosmodel<-lm(x~cos(2*pi/24*Time)+sin(2*pi/24*Time))
  fstat<-summary(cosmodel)$fstatistic
  pval<-1-pf(fstat[1],fstat[2],fstat[3])
  return(pval)
}

# Get BIC statistics of a flat model
BICFlatModel<-function(x,Experiment=NULL){
  if (is.null(Experiment)){
    return(BIC(lm(x~1)))
  }else{
    return(BIC(lm(x~Experiment)))
  }
}

# Get BIC statistics of an ANOVA model
BICANOVAModel<-function(x,Time,Experiment=NULL){
  Time<-as.factor(Time)
  if (is.null(Experiment)){
    return(BIC(lm(x~Time)))
  }else{
    return(BIC(lm(x~Experiment+ Time %in% Experiment)))
  }
}

# Get kendall tau statistics of an ANOVA model
KTANOVAModel<-function(x,Time,Experiment=NULL){
  Time<-as.factor(Time)
  if (is.null(Experiment)){
    return(cor(fitted(lm(x~Time)),x,method="kendall"))
  }else{
    return(cor(fitted(lm(x~Time + Experiment)),x[! is.na(x)],method="kendall"))
  }
}


# Get BIC statistics of a masking model
BICMaskginModel<-function(x,Time,Sleep,Wake,Experiments=NULL){
  
  if (is.null(Experiments)){
    fit<-lm(x~0+cos(2*pi/24*Time)+sin(2*pi/24*Time)+Sleep+Wake)
  }else{
    fit<-lm(x~0+Experiments+cos(2*pi/24*Time)+sin(2*pi/24*Time)+Sleep+Wake)
  }
  return(BIC(fit))
}

# Get KT statistics of a masking model
KTMaskginModel<-function(x,Time,Sleep,Wake,Experiments=NULL){
  
  if (is.null(Experiments)){
    fit<-lm(x~0+cos(2*pi/24*Time)+sin(2*pi/24*Time)+Sleep+Wake)
  }else{
    fit<-lm(x~0+Experiments+cos(2*pi/24*Time)+sin(2*pi/24*Time)+Sleep+Wake)
  }
  return(cor(fitted(fit),x[! is.na(x)],method="kendall"))
}


