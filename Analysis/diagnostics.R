###############################
##  Checking convergence for ##
## MCMC analysis in RevBayes ##
## Bruno do Rosario Petrucci ##
##       Evolution 2021      ##
###############################

###
# packages we will need

# coda
library(coda)

# readr
library(readr)

# ggplot2
library(ggplot2)

# reshape2
library(reshape2)

###
# get command line arguments
args <- commandArgs(trailingOnly = TRUE)[-1]

# take the first one - base directory for the MCMC output
baseDir <- args[1]

# and the second - base name of the output files
baseName <- args[2]

# number of parameter combinations (null and TD)
nComb <- 26

# number of reps
nRep <- 50

# ESS data frame
ess <- data.frame()

# for each parameter combination
for (i in 1:nComb) {
  # null or trait?
  type <- ifelse(i %% 2 == 0, "traits", "null")

  # name of the directory for this combination
  combDir <- paste0(baseDir, type, "_comb_", ceiling(i / 2), "/")

  # for each rep 
  for (j in 1:nRep) {
    # rep directory
    repDir <- paste0(combDir, "rep_", j, "/")

    # read log file
    log <- suppressMessages(read_tsv(paste0(repDir, baseName, ".log")))

    # burnin value - 1/4 of sample
    burnin <- ceiling(nrow(log) / 4)
    
    # make it an mcmc object
    mcmcRep <- mcmc(data = log, start = burnin)

    # get the ESS vector (excluding first two)
    essRep <- effectiveSize(mcmcRep)[-1:-2]

    # append it to the data frame
    ess <- rbind(ess, essRep)
  }
}

# plot ess 
essPlot <- ggplot(data = melt(ess), aes(x = variable, y = value)) + 
    geom_boxplot(aes(fill = variable)) +
    theme_bw()

# save it
ggsave(filename = paste0(baseDir, "ess_plot.png"), plot = essPlot)

# save ess
write_tsv(ess, paste0(baseDir, "ess.tsv"))
