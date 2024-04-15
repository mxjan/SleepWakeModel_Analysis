#!/bin/bash

#SBATCH --account pfranken_harmonic_model_ci
#SBATCH --time=6:00:00 
#SBATCH --mem-per-cpu=750M
#SBATCH --array=1-1000%100

NCHUNK=1000

module load gcc r

Rscript Human_FittingWithBootstrap.R HS_BloodTranscriptomeAndSleepWake.Rdata Human_ModelBuilding.R,Functions.R,Human_BootstrapFunctions.R /work/FAC/FBM/CIG/pfranken/harmonic_model_ci/Human_BloodGenesFitsWithBootstrapSamples_${SLURM_ARRAY_TASK_ID}_$NCHUNK.txt ${SLURM_ARRAY_TASK_ID} $NCHUNK


# 35sec by fit * 200 boot * 41000 Probes = 80'000h CPU with normal fit
# 11sec by fit * 200 boot * 41000 Probes = 25'000h CPU with normal fit
