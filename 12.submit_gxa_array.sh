#!/bin/bash
#SBATCH --partition=parallel
#SBATCH --time=108:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --array=1-10
#SBATCH --job-name=longer_brms_gxa
#SBATCH --output=logs/%x-%A_%a.log

module purge
module load R
module load HTSlib
module load GCC/13.2.0
echo "options(bitmapType='cairo')" > ~/.Rprofile

export OMP_NUM_THREADS=1

FILES=(
  model_A1_PE0.R
  model_A1_PE1.R
  model_A2_PE0.R
  model_A2_PE1.R
  model_A2_PE2.R
  model_A3_PE0.R
  model_A3_PE1.R
  model_A3_PE2.R
  model_A3_PE3.R
  character-state_model_20260505.R
)

SCRIPT=${FILES[$SLURM_ARRAY_TASK_ID - 1]}

echo "Array task : $SLURM_ARRAY_TASK_ID"
echo "Running: $SCRIPT"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

Rscript "$SCRIPT"

echo "End: $(date)"
