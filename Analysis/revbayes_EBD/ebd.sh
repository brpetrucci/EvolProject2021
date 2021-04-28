# load modules
module load gcc/7.3.0-xegsmw4
module load revbayes/1.1.0-b6wyspw

# go to the correct directory
cd /work/LAS/phylo-lab/petrucci/EvolProject2021/Analysis/revbayes_EBD/

# create a file to hole comb and rep definitions
touch aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev

# echo the definitions on it
printf "comb <- " > aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev
printf $1 >> aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev
printf "\nrep <- " >> aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev
printf $SLURM_ARRAY_TASK_ID >> aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev

# source it and the actual script
rb aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev scripts/ebd.Rev

# remove file
rm aux/aux_${1}_${SLURM_ARRAY_TASK_ID}.Rev
