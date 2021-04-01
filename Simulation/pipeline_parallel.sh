#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH --array=1-13

#SBATCH --output=job.%A_%a.out
#SBATCH --error=job.%A_%a.err

module load gcc/10.2.0-zuvaafu 
module load r/4.0.4-py3-4khjixy
Rscript --vanilla pipeline_wrapper.R --args $SLURM_ARRAY_TASK_ID
