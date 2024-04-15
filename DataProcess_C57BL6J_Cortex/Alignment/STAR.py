#!/usr/bin/python3

import os
import sys

#RawDir = "/media/md0/TimeCourseC57BL6J/RAW_fastq/"
FiltDir = "/home/maxime/Projects/Gene_Expression_Modeling_by_Sleep/Filtered_fastq/"

# Fragment size is calculated from the peak of the bioanalyser trace. It includes adapters (= 122 bp).
FragSize = {"T54S2":406,"T48S2":407,"T24S1":389,"T48S3":389,"T6N3":379,"ZT0E":411,"T48S4":407,"ZT0A":448,"T6S1":385,"ZT0C":371,
	"T30S4":387,"T3N3":419,"T12N2":380,"T12S1":411,"T12S2":373,"ZT0B":399,"T24S4":413,"T18N4":403,"T12N4":370,"T12N1":387,"T24S3":396,
	"T3S5":404,"T18S2":406,"T18N2":395,"T12S4":385,"T3N4":399,"T3N2":422,"T18S4":390,"T6S3":452,"T18N3":418,"T3S3":388,"T6N4":402,
	"T30S1":362,"T6S4":380,"T6N2":406,"T54S3":378,"T30S2":407,"T54S4":422,"T3S2":392,"T18S6":439,"T54S4":422,"T6S3":452,"T12N4":370,
	"T12S2":373,"J76S2":357,"T42S2":359,"J70N1":365,"T36S4":364,"T18N7":369,"T12N6":387,"J76S3":377,"ZT0H":373,"J70S2":360,"T42S5":384,
	"J70N2":384,"J76N2":375,"T36S3":382,"ZT0E2":354,"T36S2":364,"J76N1":364,"T12N5":382,"J76S1":384,"J70S3":371,"T18N8":371,
	"T42S3":384,"J70S1":384,"ZT0I":354}

#print(FragSize)

fdict = {}
for f in os.listdir(FiltDir):
	if 'fastq.gz' in f:
		fs = f.split("_")
		fname = fs[1]
		if fname not in fdict:
			fdict[fname] = []
		fdict[fname].append(FiltDir+f)

# Controls
for fname in FragSize:
	if fname not in fdict:
		print(fname+" in FragSize but not in files")
		sys.exit(1)
		
for fname in fdict:
	if fname not in FragSize:
		print(fname+' in files but not in FragSize')
		sys.exit(1)

# Print command:
print("module add Alignment/STAR/2.7.0e")
for fname in fdict:
	cmd = ''
	cmd += 'STAR --runThreadN 20 --genomeDir  /index/STAR/2.7.0e_mm10 '
	cmd += '--readFilesIn '+','.join(sorted(fdict[fname]))+' '
	cmd += '--readFilesCommand zcat '
	cmd += '--outFileNamePrefix /home/maxime/Projects/Gene_Expression_Modeling_by_Sleep/STAR_ReadFiltered/'+fname+' '
	cmd += '--outSAMtype BAM SortedByCoordinate '
	cmd += '--quantMode GeneCounts '
	print(cmd)
