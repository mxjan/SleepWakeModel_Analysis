# ALL PATH FOR FILES

# PROJECT DIRECTORY
#PROJECTDIR<-"C:/Users/Maxime Jan/Projects/SleepWake_Model/"
PROJECTDIR<-"E:/Projects/SleepWake_Model/"
ANALYSISDIR<-paste(PROJECTDIR,"Analysis_SleepWakeModel/",sep="")
DATADIR<-paste(PROJECTDIR,"Data/",sep="")
FIGUREDIR<-paste(PROJECTDIR,"Analysis_SleepWakeModel/Figures/",sep="")


######### ANALYSIS FILES ##########
COLORFUN<-paste(ANALYSISDIR,"RFunctions/color.R",sep="") # Contain all colors for liver, cortex, sleep-wake, circadian, Ztime etc...
MOUSEMODELS_FUN<-paste(ANALYSISDIR,"RFunctions/C57BL6_ModelBuilding.R",sep="")
HUMANMODELS_FUN<-paste(ANALYSISDIR,"RFunctions/Human_ModelBuilding.R",sep="")
RASTERPLOT_FUN<-paste(ANALYSISDIR,"RFunctions/RasterPlot.R",sep="")
LINEARPLOT_FUN<-paste(ANALYSISDIR,"RFunctions/PlotFits.R",sep="")
TEMPLATE_FUN<-paste(ANALYSISDIR,"RFunctions/Templates.R",sep="")
GENERAL_FUN<-paste(ANALYSISDIR,"RFunctions/Functions.R",sep="")
PCA_FUN<-paste(ANALYSISDIR,"RFunctions/PCA_Functions.R",sep="")
TOPGO_FUN<-paste(ANALYSISDIR,"RFunctions/EnrichmentAnalysis.R",sep="")
TIMESIGNATR_FUN<-paste(ANALYSISDIR,"RFunctions/TimeStampFns.R",sep="")
TISSUESYNC_FUN<-paste(ANALYSISDIR,"RFunctions/TissueSync_Functions.R",sep="")
RECOVERY_FUN<-paste(ANALYSISDIR,"RFunctions/Recovery_Functions.R",sep="")

######### DATA FILES #########
# CORTEX
CORTEXDATADIR<-paste(DATADIR,"C57BL6J_Cortex_TimeCourse_CHN2019/",sep="")
CORTEXMETADATA<-paste(CORTEXDATADIR,"CHN2019_metadata.txt",sep="")
CORTEXREADCOUNTSDIR<-paste(CORTEXDATADIR,"STAR_mm10_ReadCounts/",sep="")
CORTEX_METADATA_RDATA<-paste(CORTEXDATADIR,"CHN2019_TimeCourse_MetaData.RData",sep="")
CORTEX_NORMEXPR_RDATA<-paste(CORTEXDATADIR,"CHN2019_TimeCourse_NormalizedData_Genes.RData",sep="")
CORTEX_SWDMRDATA_RDATA<-paste(CORTEXDATADIR,"B6Cortex_FilePrepSWDMr.Rdata",sep="")

CORTEXFITS<-paste(CORTEXDATADIR,"C57BL6_CortexGenesFits.txt",sep="") # Full model
CORTEXFITSBOOT<-paste(CORTEXDATADIR,"C57BL6_CortexGenesFitsWithBootstrap.txt.gz",sep="") # Full model
CORTEXFITS_DropSin<-paste(CORTEXDATADIR,"C57BL6_CortexGenesFits_DropSin_061221.txt",sep="") # Driven by Sleep-wake only
CORTEXFITS_DropSW<-paste(CORTEXDATADIR,"C57BL6_CortexGenesFits_DropSW_061221.txt",sep="") # Driven by Circadian force only

CORTEXFITTED_RDATA<-paste(CORTEXDATADIR,"CortexFittedValues.Rdata",sep="") # Fitted data of the full model 

# LIVER
LIVERDATADIR<-paste(DATADIR,"C57BL6J_Liver_TimeCourse_MJCHN2020/",sep="")
LIVERMETADATA<-paste(LIVERDATADIR,"MetaData.txt",sep="")
LIVERREADCOUNTSDIR<-paste(LIVERDATADIR,"STAR_mm10_ReadCounts/",sep="")
LIVER_METADATA_RDATA<-paste(LIVERDATADIR,"MJCHN2020Liver_TimeCourse_MetaData.RData",sep="")
LIVER_NORMEXPR_RDATA<-paste(LIVERDATADIR,"MJCHN2020Liver_TimeCourse_NormalizedData_Genes.RData",sep="")
LIVER_SWDMRDATA_RDATA<-paste(LIVERDATADIR,"B6Liver_FilePrepSWDMr.Rdata",sep="")

LIVERFITS<-paste(LIVERDATADIR,"C57BL6_LiverGenesFits.txt",sep="") # Full model
LIVERFITSBOOT<-paste(LIVERDATADIR,"C57BL6_LiverGenesFitsWithBootstrap.txt.gz",sep="") # Full model
LIVERFITS_DropSin<-paste(LIVERDATADIR,"C57BL6_LiverGenesFits_DropSin_061221.txt",sep="") # Driven by Sleep-wake only
LIVERFITS_DropSW<-paste(LIVERDATADIR,"C57BL6_LiverGenesFits_DropSW_061221.txt",sep="") # Driven by Circadian force only

LIVERFITTED_RDATA<-paste(LIVERDATADIR,"LiverFittedValues.Rdata",sep="") # Fitted data of the full model saved

# Mouse Sleep-Wake
MOUSESLEEPWAKE_DIR<-paste(DATADIR,"BXD_SleepWake_SDMJ2018/",sep="")
GO_ANN<-paste(DATADIR,"GOannot.txt",sep="")

# Human Blood Constant Routine
HUMANCRDATADIR<-paste(DATADIR,"Human_ConstantRoutine/",sep="")
HUMANCR_MICROARRAY<-paste(HUMANCRDATADIR,"GSE39445/",sep="")
HUMANCR_MICROARRAY_CELL<-paste(HUMANCR_MICROARRAY,"GSE39445/",sep="")
HUMANCR_MICROARRAY_SOFT<-paste(HUMANCR_MICROARRAY,"GSE39445_SOFT/",sep="")

# Human Blood Forced Desynchrony
HUMANFDDATADIR<-paste(DATADIR,"Human_ForcedDesynchrony/",sep="")
HUMANFD_MICROARRAY<-paste(HUMANFDDATADIR,"GSE48113/",sep="")
HUMANFD_MICROARRAY_CELL<-paste(HUMANFD_MICROARRAY,"GSE48113/",sep="")
HUMANFD_MICROARRAY_SOFT<-paste(HUMANFD_MICROARRAY,"GSE48113_SOFT/",sep="")

# Human transcriptome and sleep directory and data once processed
HUMANPROCESSEDDIR<-paste(DATADIR,"Human_ProcessData/",sep="")
HUMANTRANSCRIPTOME_RDATA<-paste(HUMANPROCESSEDDIR,"HT_NormalizedData.Rdata",sep="") # First normalization
HUMANTRANSCRIPTOMECORRECT_RDATA<-paste(HUMANPROCESSEDDIR,"HT_NormCorrData.Rdata",sep="") # Correct for subject effect
HUMANTRANSCRIPTOMEEMEANS_RDATA<-paste(HUMANPROCESSEDDIR,"Data4Emmeans.Rdata",sep="") # used for visualization
HUMANTRANSCRIPTOMESWDMR_RDATA<-paste(HUMANPROCESSEDDIR,"HT_SWDMr.Rdata",sep="") # will be used for fitting

HUMANCR_SLEEPDIR<-paste(HUMANCRDATADIR,"CRC259/",sep="")
HUMANCR_SLEEP_RDATA<-paste(HUMANPROCESSEDDIR,"CRC259_SWdf.Rdata",sep="")
HUMANFD_SLEEPDIR<-paste(HUMANFDDATADIR,"CRC262/",sep="")
HUMANFD_SLEEP_RDATA<-paste(HUMANPROCESSEDDIR,"CRC262_SWdf.Rdata",sep="")
HUMANFD_CIRCADIANINFO_RDATA<-paste(HUMANPROCESSEDDIR,"CircadianInfo.Rdata",sep="")

HUMANDATASET_RDATA<-paste(HUMANPROCESSEDDIR,"HS_BloodTranscriptomeAndSleepWake.Rdata",sep="")

HUMANFITS_DropSin<-paste(HUMANPROCESSEDDIR,"Human_BloodProbesFits_DropSin_061221.txt",sep="")
HUMANFITS_DropSW<-paste(HUMANPROCESSEDDIR,"Human_BloodProbesFits_DropSW_061221.txt",sep="")
HUMANFITS<-paste(HUMANPROCESSEDDIR,"Human_BloodGenesFits.txt",sep="")
HUMANFITSBOOT<-paste(HUMANPROCESSEDDIR,"Human_BloodGenesFitsWithBootstrapSamples.txt.gz",sep="") # Full model

HUMANFITTED_RDATA<-paste(HUMANPROCESSEDDIR,"HumanFittedValues.Rdata",sep="")

ARCHER2014FIG6<-paste(HUMANPROCESSEDDIR,"Results_Fig6_Archer2014.csv",sep="")

