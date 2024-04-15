args = commandArgs(trailingOnly=TRUE) 

# Args 1 is Rdata
ExprMetaF<-args[1]
# Args 2 is models
FunSourceF<-args[2]
# Args 3 is file output
OF<-args[3]
# Args 4 is ncors number
NCORES<-as.numeric(args[4])

# Example: 
# Rscript C57BL6_Fitting_Full.R B6Cortex_FilePrepSWDMr.Rdata C57BL6_ModelBuilding.R C57BL6_CortexGenesFits.txt 10

################## RUN #################

library(SWDMr)
library(optimx)
library(doSNOW)

########### Functions

# Write returned parameters from foreach
wparams <- function(outF, d) {
  writeLines(paste(d$Gene,d$params$Wake,d$params$Sleep,d$params$loggamma,d$params$omega,
                   d$params$RSS,d$params$BF,d$params$BIC,d$params$KendallTau,
                   sep="\t"), outF)
  return(outF)
}

# Fit datas and return fitting values to write
FitData<-function(swdmr,Gene){
  
  # Full model
  full<-FitDropSinF(swdmr,Gene)
  
  retv<-c()
  retv<-c(retv,full$optimxres[1,c("Wake","Sleep","loggamma","omega")])
  retv<-c(retv,c(full$stats$stats$RSS, full$stats$stats$BayesFactor, full$stats$stats$BIC, full$stats$stats$KendalTau ))
  names(retv)[seq(5,8)]<-c("RSS","BF","BIC","KendallTau")
  return(retv)
}

###########

# Load data
load(ExprMetaF)

# Load model building functions for mouse fitting
source(FunSourceF)

swdmr <- SWDMr(SWdist=SWdf, Gexp=rna_expr)

# Get All Genes
Genes<-colnames(swdmr@Gexp)[! colnames(swdmr@Gexp) == "Time"]

# # Reduce genes number for testing
#Genes<-Genes[1:30]

# Write output file header
outF <- file(OF)
writeLines(paste("Gene","Wake","Sleep","loggamma","omega","RSS","BF","BIC","KendallTau",
                 sep="\t"), outF)
close(outF)

# Add lines to output file
outF <- file(OF, "a")

# Create parallel Sockets
cl <- makeCluster(NCORES)
clusterExport(cl,c("swdmr"))
registerDoSNOW(cl)

# Fit for all Genes
fits <- foreach(i = 1:length(Genes),.init=outF,.packages=c("optimx","SWDMr"),.combine="wparams") %dopar% {
  
  Gene<-Genes[i]
  
  tf<-try(fit<-FitData(swdmr,Gene))
  
  if (is(tf,"try-error")){fit<-list(Wake=NA,Sleep=NA,loggamma=NA,omega=NA,RSS=NA,BF=NA,BIC=NA,
                                    KendallTau=NA)}
  
  return(list(Gene=Gene,params=fit))
  
}

stopCluster(cl)

close(outF)
