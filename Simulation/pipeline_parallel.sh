#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --partition=biocrunch # partition (from guide in pronto website)
#SBATCH --time=02-00:00:00 # 48h - probably will take less but just to be safe
#SBATCH --array=1-13 # array of jobs (13 parameter combinations)

#SBATCH --output=job.%A_%a.out # names have master job id and array job id
#SBATCH --error=job.%A_%a.err

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd /work/LAS/phylo-lab/petrucci/EvolProject2021/Simulation # go to the correct directory

module load gcc/10.2.0-zuvaafu 
module load r/4.0.4-py3-4khjixy # load required modules

Rscript --vanilla pipeline_wrapper.R --args $SLURM_ARRAY_TASK_ID # run job number x
