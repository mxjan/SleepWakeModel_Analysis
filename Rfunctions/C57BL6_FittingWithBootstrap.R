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
# Rscript C57BL6_FittingWithBootstrap.R B6Cortex_FilePrepSWDMr.Rdata C57BL6_ModelBuilding.R,Functions.R,C57BL6_BootstrapFunctions.R C57BL6_CortexGenesFitsWithBootstrap_1_2125.txt 1 2125

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

swdmr <- SWDMr(SWdist=SWdf, Gexp=rna_expr)

# Get All Genes
Genes<-colnames(swdmr@Gexp)[! colnames(swdmr@Gexp) == "Time"]
# Split Genes
GenesSplit<-split(Genes,  cut(seq_along(Genes), NSplit, labels = FALSE))
# Use only selected Gene
Genes<-GenesSplit[[SplitToDO]]
print(paste("N genes in chunk:",length(Genes)))

# Write output file header
outF <- file(OF)
writeLines(paste("Gene","FittingType", # Gene and DataFitting or BootstrapX
                 "Wake","Sleep","loggamma","omega","AmpSin","PhiSin", # Parameters
                 "tau","SWrc","zeta","phi_lag", # derived paramters
                 "RSS","BF","BIC","BIC_null","KendallTau", # statistics
                 sep="\t"), outF)
close(outF)

# Add lines to output file
outF <- file(OF, "a")

# Function used to write data
wparams <- function(outF, df) {
  
  for (j in 1:nrow(df)){
    writeLines(paste(df$Gene[j],df$FittingType[j],
                     df$Wake[j],df$Sleep[j],df$loggamma[j],df$omega[j],df$AmpSin[j],df$PhiSin[j],
                     df$tau[j],df$SWrc[j],df$zeta[j],df$phi_lag[j],
                     df$RSS[j],df$BF[j],df$BIC[j],df$BIC_null[j],df$KendallTau[j],
                     sep="\t"), outF)
  }
  return(outF)
}


# Fit for all Genes in chunck
fits <- foreach(i = 1:length(Genes),.init=outF,.packages=c("optimx","SWDMr"),.combine="wparams") %do% {

  Gene<-Genes[i]
  
  print(paste("DO GENE:",Gene))

  # Fit data
  tf<-try(fit<-FitFullModel(swdmr,Gene))
  if (is(tf,"try-error")){
    DataFitting_res<-list(Gene=Gene,FittingType="DataFitting",
                          Wake=NA,Sleep=NA,loggamma=NA,
                          omega=NA,AmpSin=NA,PhiSin=NA,
                          tau=NA,SWrc=NA,zeta=NA,phi_lag=NA,
                          RSS=NA,BF=NA,
                          BIC=NA,BIC_null=NA,
                          KendallTau=NA)
    results<-do.call("rbind.data.frame",c(list(DataFitting_res)))
  }else{
    # Add derived paramters
    optimxres<-fit$optimxres
    optimxres$Gene<-Gene
    fit$stats$stats$tau<-Gettau(Gene,optimxres)
    fit$stats$stats$SWrc<-GetSWrcMouse(Gene,swdmr,fits = optimxres)
    fit$stats$stats$zeta<-GetZeta(Gene,optimxres)
    fit$stats$stats$phi_lag<-GetPhaseLag(Gene,optimxres,PosAmp = T)
    DataFitting_res<-list(Gene=Gene,FittingType="DataFitting",
                Wake=optimxres$Wake,Sleep=optimxres$Sleep,loggamma=optimxres$loggamma,
                omega=optimxres$omega,AmpSin=optimxres$AmpSin,PhiSin=optimxres$PhiSin,
                tau=fit$stats$stats$tau,SWrc=fit$stats$stats$SWrc,zeta=fit$stats$stats$zeta,phi_lag=fit$stats$stats$phi_lag,
                RSS=fit$stats$stats$RSS,BF=fit$stats$stats$BayesFactor,
                BIC=fit$stats$stats$BIC,BIC_null=fit$stats$stats$BIC_flat,
                KendallTau=fit$stats$stats$KendalTau)

    # Run non_parameteric bootstrap of residuals with replacement
    res_nonparametricR<-lapply(1:nboot,function(B){DoBootstrap(B,fit = fit,params = optimxres,
                                                               swdmr,Gene,boot="samples",parametric = F,Replace = T,speed = "fast")})

    results<-do.call("rbind.data.frame",c(list(DataFitting_res),res_nonparametricR))
  }

  return(results)
}

close(outF)


