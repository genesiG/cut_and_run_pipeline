#!/bin/bash
#BSUB -P GG_CUTANDTAG_NUCLEATION
#BSUB -J run_R_steps
#BSUB -o Analysis_Data/peaks/seacr/log/run_R.log
#BSUB -e Analysis_Data/peaks/seacr/log/run_R.error
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 2
#BSUB -M 32000
#BSUB -R "rusage[mem=16000] span[hosts=1]"

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "--- Running step_5e ---"
Rscript --vanilla Scripts/step_5e_define_nucleation_spreading.R

echo "--- Running step_7 ---"
Rscript --vanilla Scripts/step_7_inhibitor_resistance.R

echo "--- DONE ---"
