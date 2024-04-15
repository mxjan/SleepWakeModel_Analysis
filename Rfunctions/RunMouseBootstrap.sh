#!/bin/bash

#SBATCH --account pfranken_harmonic_model_ci
#SBATCH --time=15:00:00 
#SBATCH --mem-per-cpu=500M
#SBATCH --array=1-100

NCHUNK=100

module load gcc r

Rscript C57BL6_FittingWithBootstrap.R B6Cortex_FilePrepSWDMr.Rdata C57BL6_ModelBuilding.R,Functions.R,C57BL6_BootstrapFunctions.R /work/FAC/FBM/CIG/pfranken/harmonic_model_ci/C57BL6_CortexGenesFitsWithBootstrapSamples_${SLURM_ARRAY_TASK_ID}_$NCHUNK.txt ${SLURM_ARRAY_TASK_ID} $NCHUNK

Rscript C57BL6_FittingWithBootstrap.R B6Liver_FilePrepSWDMr.Rdata C57BL6_ModelBuilding.R,Functions.R,C57BL6_BootstrapFunctions.R /work/FAC/FBM/CIG/pfranken/harmonic_model_ci/C57BL6_LiverGenesFitsWithBootstrapSamples_${SLURM_ARRAY_TASK_ID}_$NCHUNK.txt ${SLURM_ARRAY_TASK_ID} $NCHUNK


