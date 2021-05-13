#################################
## Summarizing results of MCMC ##
##  Bruno do Rosario Petrucci  ##
##        Evolution 2021       ##
#################################

###
# packages we will need

# coda
library(coda)

# readr
library(readr)

# tidyverse
library(tidyverse)

# reshape2
library(reshape2)

# source hpd funtion
# did not want to have to download the whole package
source("/work/LAS/phylo-lab/petrucci/EvolProject2021/Analysis/hpd.R")

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

# data frame for summary numbers
summary <- data.frame()

# for each parameter combination
for (i in 1:nComb) {
  # print the comb we are at
  print(i)

  # null or trait?
  type <- ifelse(i %% 2 == 0, "traits", "null")

  # comb number
  nCombCur <- ceiling(i / 2)

  # name of the directory for this combination
  combDir <- paste0(baseDir, type, "_comb_", nCombCur, "/")

  # for each rep 
  for (j in 1:nRep) {
    # rep directory
    repDir <- paste0(combDir, "rep_", j, "/")

    # read log file
    log <- suppressMessages(read_tsv(paste0(repDir, baseName, ".log")))

    # burnin value - 1/2 of sample
    burnin <- ceiling(nrow(log) / 2)

    # apply burnin
    log <- log[(burnin + 1):nrow(log), -2:-5]

    # make log a data frame
    log <- as.data.frame(log)

    # vector for this rep
    summ <- c(nCombCur, type, j)

    # get HPD intervals
    hpds <- suppressWarnings(hpd(log))

    # for each column in the log
    for (k in 2:ncol(log)) {
      # get the data
      sample <- log[, k]
      
      # append the mean
      summ <- c(summ, mean(sample))

      # append the low, median, and high HPD
      summ <- c(summ, hpds[k - 1, 1], median(sample), hpds[k - 1, 2])
    }

    # append to data frame
    summary <- rbind(summary, summ)
  }
}

# names for parameter columns
parColNames <- unlist(lapply(2:ncol(log), function(x)
  paste0(colnames(log)[x], c("_mean", "_low95", "_median", "_high95"))))

# column names
colnames(summary) <- c("comb", "type", "rep", parColNames)

# save data frame
write_tsv(summary, paste0(baseDir, "summary.tsv"))
