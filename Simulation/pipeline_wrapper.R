# source pipeline (parameter definitions, simulation functions etc)
source(file = "sim_pipeline.R")

# recover arguments
args <- commandArgs(trailingOnly = TRUE)

###
# check arguments

# need to have at least one
if (length(args) < 1) {
  stop("At least one argument must be suplied.")
}

# cannot have more than 3
if (length(args) > 3) {
  stop("Too many arguments. Please supply a comb and optionally
       an nReps argument only.")
}

# make the first argument numeric
arg1 <- as.numeric(args[1])

# assumed to be the combination of parameters
if (arg1 > nrow(key) || arg1 < 1 ||
    round(floor(arg1) - arg1) != 0) {
  stop("First argument is not a valid combination of parameters.")
}

# set comb to arg1
comb <- arg1

# if there is only one, set nReps to default 50 and simDir 
# to default /work/LAS/phylo-lab/petrucci/EvolProject2021/Simulation/replicates
if (length(args) == 1) {
  nReps <- 50
  
  # create default simDir if it doesn't exist
  simDir <- "/work/LAS/phylo-lab/petrucci/EvolProject2021/Simulation/replicates"
  smart.dir.create(simDir)
} else if (length(args) == 2) {
  # assume the second argument is nReps
  nReps <- as.numeric(args[2])
  
  # create simDir if it doesn't exist
  simDir <- "/work/LAS/phylo-lab/petrucci/EvolProject2021/Simulation/replicates"
  smart.dir.create(simDir)
} else {
  nReps <- as.numeric(args[2])
  simDir <- args[3]
  
  # create simDir if it doesn't exist
  smart.dir.create(simDir)
}

# save key
write_tsv(key, file = paste0(simDir, "key.tsv"))

# run simulations
simulate(nReps, comb, key, simDir)