# Rename files for RNA-seq SDrec from liver

FASTQ_DIR<-"/media/md0/TimeCourseC57BL6J/Filtered_fastq_Liver/"
OUT_DIR<-"/home/maxime/Projects/RNA_SDrec_Liver/STAR_ReadFiltered/"

Meta_D<-read.table("MetaData.txt",header=T,stringsAsFactors = F)
Meta_D$FileID <- as.character(Meta_D$FileID)

cat("module add Alignment/STAR/2.7.0e")
cat("\n")

for (ID in Meta_D$FileID){
  
  # unmodified ID
  oID<-ID
  
  if (nchar(ID) == 1){ ID <- paste("0",ID,sep="")}
  fl <- list.files(path = FASTQ_DIR,pattern = paste("^",ID,"_L",sep=""))
  # print(ID)
  # print(fl)
  
  cmd <- c("")
  cmd <- paste(cmd , "STAR --runThreadN 20 --genomeDir  /index/STAR/2.7.0e_mm10")
  cmd <- paste(cmd , "--readFilesIn" , paste(paste(FASTQ_DIR,fl,sep=""),collapse = ",") )
  cmd <- paste(cmd , "--readFilesCommand zcat" )
  cmd <- paste(cmd , "--outFileNamePrefix" , paste(OUT_DIR,Meta_D$SampleID[Meta_D$FileID == oID],sep=""))
  cmd <- paste(cmd , "--outSAMtype BAM SortedByCoordinate")
  cmd <- paste(cmd , "--quantMode GeneCounts")
  cat(cmd)
  cat("\n")
}
