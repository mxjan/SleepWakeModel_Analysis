#!/usr/bin/python3

# Run using  ~/Projects/Gene_Expression_Modeling_by_Sleep/filter_fastq.py > filter.sh
# nohup parallel -j 10 < filter.sh &

import os

RawDir = "/media/md0/TimeCourseC57BL6J/RAW_fastq_Liver/"
FiltDir = "/media/md0/TimeCourseC57BL6J/Filtered_fastq_Liver/"

flist = []
for f in os.listdir(RawDir):
	if 'fastq.gz' in f:
		flist.append(f)

# Print commands
print("module load Filtering/fastq_illumina_filter-0.1")
for f in flist:
	cmd = ''
	cmd += 'zcat '+RawDir+f+' '
	cmd += '| fastq_illumina_filter '
	cmd += '--keep N -v '
	cmd += '| gzip > '+FiltDir+f.replace(".fastq.gz",'.filtered.fastq.gz')+' '
	print(cmd)
