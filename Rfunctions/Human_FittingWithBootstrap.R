# RUN SBATCH AS ARRAY
args = commandArgs(trailingOnly=TRUE) 

# Args 1 is Rdata with Gene expression and Sleep-Wake data.frame
ExprMetaF<-args[1]

# Args 2 is models functions and required function, comma separated
FunSourceF<-args[2]

# Args 3 is file output
OF<-args[3]

# The split to do
SplitToDO<-as.numeric(args[4])

# Number of split
NSplit<-as.numeric(args[5])

# 200 bootstrap are done, edit if you want to change
nboot<-200

# Example: 
# Rscript Human_FittingWithBootstrap.R HS_BloodTranscriptomeAndSleepWake.Rdata Human_ModelBuilding.R,Functions.R,Human_BootstrapFunctions.R Human_BloodGenesFitsWithBootstrap_1_2125.txt 1 2125

print(paste("Rdata frame given:",ExprMetaF))
print(paste("Rfunctions given:",FunSourceF))
print(paste("Output:",OF))
print(paste("Number of chunk to do:",NSplit))
print(paste("chunk analyzed:",SplitToDO))
print(paste("Number of bootstrap:",nboot))

################## RUN #################

library(SWDMr)
library(optimx)
library(foreach)
library(dplyr)

# Load data
load(ExprMetaF)

# Load model building functions for mouse fitting
for (funs in strsplit(FunSourceF,",")[[1]]){
  print(paste("Get Functions from:",funs))
  source(funs)
}


# Get All Probes
Probes<-colnames(HT_Desync)[which(! colnames(HT_Desync) %in% c("Time", "Subject"))]
# Split Genes
ProbesSplit<-split(Probes,  cut(seq_along(Probes), NSplit, labels = FALSE))
# Use only selected Gene
Probes<-ProbesSplit[[SplitToDO]]
print(paste("N Probes in chunk:",length(Probes)))

# Write output file header
outF <- file(OF)
writeLines(paste("ProbeName","Gene","FittingType", # Gene and DataFitting or BootstrapX
                 "Wake","Sleep","loggamma","omega","AmpSin","PhiSin","intercept", # Parameters
                 "tau","SWrc","zeta","phi_lag", # derived paramters
                 "RSS","BF","BIC","BIC_null","KendallTau", # statistics
                 sep="\t"), outF)
close(outF)

# Add lines to output file
outF <- file(OF, "a")

# Function used to write data
wparams <- function(outF, df) {
  
  for (j in 1:nrow(df)){
    writeLines(paste(df$ProbeName[j],df$Gene[j],df$FittingType[j],
                     df$Wake[j],df$Sleep[j],df$loggamma[j],df$omega[j],df$AmpSin[j],df$PhiSin[j],df$intercept[j],
                     df$tau[j],df$SWrc[j],df$zeta[j],df$phi_lag[j],
                     df$RSS[j],df$BF[j],df$BIC[j],df$BIC_null[j],df$KendallTau[j],
                     sep="\t"), outF)
  }
  return(outF)
}


# Fit for all Genes in chunck
fits <- foreach(i = 1:length(Probes),.init=outF,.packages=c("optimx","SWDMr"),.combine="wparams") %do% {

  ProbeID<-Probes[i]
  Gene<-ExprMeta[ExprMeta$ProbeName %in% ProbeID,"GeneName"][1]
  
  print(paste("DO PROBE:",ProbeID,"[",Gene,"]"))

  # Fit data
  tf<-try(fit<-FitDataFull(ProbeID , MeanSWdf, HT_Desync, MeanSWdfExt, HT_Ext, MeanSWdfRest, HT_Res))
  if (is(tf,"try-error")){
    DataFitting_res<-list(ProbeName=ProbeID,Gene=Gene,FittingType="DataFitting",
                          Wake=NA,Sleep=NA,loggamma=NA,
                          omega=NA,AmpSin=NA,PhiSin=NA,intercept=NA,
                          tau=NA,SWrc=NA,zeta=NA,phi_lag=NA,
                          RSS=NA,BF=NA,
                          BIC=NA,BIC_null=NA,
                          KendallTau=NA)
    results<-do.call("rbind.data.frame",c(list(DataFitting_res)))
  }else{
    # Add derived paramters
    optimxres<-fit$paramsExtRes
    optimxres$Gene<-Gene
    optimxres$ProbeName<-ProbeID
    fit$stats$tau<-Gettau(ProbeID,optimxres)
    fit$stats$SWrc<-GetSWrcHuman(ProbeID,MeanSWdf,HT_Desync = HT_Desync,fits = optimxres)
    fit$stats$zeta<-GetZeta(ProbeID,optimxres)
    fit$stats$phi_lag<-GetPhaseLag(ProbeID,optimxres,PosAmp = T)
    DataFitting_res<-list(ProbeName=ProbeID,Gene=Gene,FittingType="DataFitting",
                Wake=optimxres$Wake,Sleep=optimxres$Sleep,loggamma=optimxres$loggamma,
                omega=optimxres$omega,AmpSin=optimxres$AmpSin,PhiSin=optimxres$PhiSin,intercept=optimxres$intercept,
                tau=fit$stats$tau,SWrc=fit$stats$SWrc,zeta=fit$stats$zeta,phi_lag=fit$stats$phi_lag,
                RSS=fit$stats$RSS_model,BF=fit$stats$BF,
                BIC=fit$stats$BIC_model,BIC_null=fit$stats$BIC_flat,
                KendallTau=fit$stats$KendallTau)

    if (fit$stats$BF >= exp(1)){
      # Run non_parameteric bootstrap of residuals with replacement
      res_nonparametricR<-lapply(1:nboot,function(B){DoBootstrap(B,fit = fit,params = optimxres,Gene=Gene,
                                                                 ProbeID=ProbeID,boot="TSsamples",parametric = F,Replace = T,speed = "fast")})
      
      results<-do.call("rbind.data.frame",c(list(DataFitting_res),res_nonparametricR)) 
    }else{
      results<-do.call("rbind.data.frame",c(list(DataFitting_res)))
    }
  }

  return(results)
}

close(outF)


