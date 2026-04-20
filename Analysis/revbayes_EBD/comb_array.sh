#!/bin/bash

#SBATCH --nodes=1 # 1 node per job
#SBATCH --time=02-00:00:00 # 48h - probably will take less but just to be safe
#SBATCH --array=1-13 # array of jobs (13 parameter combinations)

#SBATCH --output=output/jobs/array/job_%A_%a.out
#SBATCH --output=output/jobs/array/job_%A_%a.out

#SBATCH --job-name="array_EBD"

#SBATCH --mail-user=petrucci@iastate.edu   # my e-mail
#SBATCH --mail-type=BEGIN # get notifications for all job cases
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

mkdir output/jobs/job_${SLURM_ARRAY_TASK_ID}
sbatch --time=2-00:00:00 --array=1-100 --output="output/jobs/job_${SLURM_ARRAY_TASK_ID}/ebd_${SLURM_ARRAY_TASK_ID}_%A_%a.out" --error="output/jobs/job_${SLURM_ARRAY_TASK_ID}/ebd_${SLURM_ARRAY_TASK_ID}_%A_%a.err" --mail-user=petrucci@iastate.edu --mail-type=BEGIN --mail-type=END --mail-type=FAIL --job-name="${SLURM_ARRAY_TASK_ID}_array" --wrap="sh ebd.sh ${SLURM_ARRAY_TASK_ID}"
