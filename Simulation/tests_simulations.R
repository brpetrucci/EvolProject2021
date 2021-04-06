##################################
##  Tests on simulation results ##
##     Evolution 2021 Project   ## 
##   Bruno do Rosario Petrucci  ##
##################################

# need readr for the read_tsv function
library(readr)

# and ape to read nexus data
library(ape)

# and geiger to estimate bm variance
library(geiger)

# devtools and pb
library(devtools)
load_all()

# directory where all simulations are
simDir <- "/Users/petrucci/Documents/research/EvolProject2021/Simulation/replicates/"

# create objects to hold our diagnostics

# matrix to hold total number of species for null and trait sims
nullN <- traitN <- matrix(0, nrow = 13, ncol = 50)

# and same for extant number
nullExtantN <- traitExtantN <- matrix(0, nrow = 13, ncol = 50)

# rho for null and trait too
nullRho <- traitRho <- matrix(0, nrow = 13, ncol = 50)

# bm trait mean
nullBMMean <- traitBMMean <- matrix(0, nrow = 13, ncol = 50)

# and variance
nullBMVar <- traitBMVar <- matrix(0, nrow = 13, ncol = 50)

# same for discrete
nullSTMean <- traitSTMean <- nullSTVar <- traitSTVar <-
  matrix(0, nrow = 13, ncol = 50)


bmVals <- function(x, simList, traitFunc) {
  # list of trait values at the time of simulation's end
  traitList <- as.numeric(unlist(lapply(which(simList[[x]]$EXTANT), function(y) 
    traitFunc[[x]][[paste0("t", y)]](max(simList[[x]]$TS)))))
  
  if (length(traitList) < 5) return(NA)
  
  names(traitList) <- paste0("t", which(simList[[x]]$EXTANT))
  
  # phylogeny to get MRCA time
  phy <- drop.fossil(make.phylo(simList[[x]]))
  phy$root.time <- NULL
  
  # estimate BM
  estimate <- fitContinuous(phy, traitList)
  
  return(list(mean = estimate$opt$z0, var = estimate$opt$sigsq))
}

# and the discrete ratios
stVals <- function(x, simList, traitFunc) {
  traitList <- as.numeric(unlist(lapply(which(simList[[x]]$EXTANT), function(y) 
    traitFunc[[x]][[paste0("t", y)]](max(simList[[x]]$TS)))))
  
  names(traitList) <- paste0("t", which(simList[[x]]$EXTANT))
  
  return(list(mean = mean(traitList), var = var(traitList)))
}

# loop through combinations of parameters
for (i in 1:13) {
  print(i)
  
  # directory for the data
  combDir <- paste0(simDir, "comb_", i)
  
  # null directory
  nullDir <- paste0(combDir, "/null/")
  
  # and traits
  traitsDir <- paste0(combDir, "/traits/")
  
  # load sim lists
  load(paste0(nullDir, "sim_list.RData"))
  simListNull <- simList
  
  load(paste0(traitsDir, "sim_list.RData"))
  
  # populate nullN and nullExtantN
  nullN[i, ] <- unlist(lapply(1:length(simListNull), function(x)
    length(simListNull[[x]]$TS)))
  nullExtantN[i, ] <- unlist(lapply(1:length(simListNull), function(x)
    sum(simListNull[[x]]$EXTANT)))
  
  # and traitN and traitExtantN
  traitN[i, ] <- unlist(lapply(1:length(simList), function(x)
    length(simList[[x]]$TS)))
  traitExtantN[i, ] <- unlist(lapply(1:length(simList), function(x)
    sum(simList[[x]]$EXTANT)))
  
  # get function lists for null
  load(paste0(nullDir, "bm_traits.RData"))
  load(paste0(nullDir, "st_traits.RData"))
  
  nullBMTraitsFunc <- bmTraitsFunc
  nullSTTraitsFunc <- stTraitsFunc
  
  # and for traits
  load(paste0(traitsDir, "bm_traits.RData"))
  load(paste0(traitsDir, "st_traits.RData"))
  
  # get null bm trait values
  nullBMVals <- lapply(1:50, function(x)
    bmVals(x, simListNull, nullBMTraitsFunc))
  nullBMMean[i, ] <- unlist(lapply(1:50, function(x) 
    nullBMVals[[x]]$mean))
  nullBMVar[i, ] <- unlist(lapply(1:50, function(x) 
    nullBMVals[[x]]$var))
  
  # get trait bm trait values
  traitBMVals <- lapply(1:50, function(x)
    bmVals(x, simList, bmTraitsFunc))
  traitBMMean[i, ] <- unlist(lapply(1:50, function(x) 
    traitBMVals[[x]]$mean))
  traitBMVar[i, ] <- unlist(lapply(1:50, function(x) 
    traitBMVals[[x]]$var))
  
  # get null st trait values
  nullSTVals <- lapply(1:50, function(x)
    stVals(x, simListNull, nullSTTraitsFunc))
  nullSTMean[i, ] <- unlist(lapply(1:50, function(x) 
    nullSTVals[[x]]$mean))
  nullSTVar[i, ] <- unlist(lapply(1:50, function(x) 
    nullSTVals[[x]]$var))
  
  # get trait st trait values
  traitSTVals <- lapply(1:50, function(x)
    stVals(x, simList, stTraitsFunc))
  traitSTMean[i, ] <- unlist(lapply(1:50, function(x) 
    traitSTVals[[x]]$mean))
  traitSTVar[i, ] <- unlist(lapply(1:50, function(x) 
    traitSTVals[[x]]$var))
}

boxplot.matrix(nullBMMean, use.cols = FALSE, outline = FALSE)
abline(h=0)
boxplot.matrix(traitBMMean, use.cols = FALSE, outline = FALSE, add = TRUE, col = "green")

boxplot.matrix(nullBMVar, use.cols = FALSE, outline = FALSE)
abline(h = 0.021)
abline(h = 0.031, col = 'red')
abline(h = 0.011, col = 'blue')
boxplot.matrix(traitBMVar, use.cols = FALSE, outline = FALSE, col = "green", add = TRUE)

expected.trait <- function(q01, q10, t) {
  q <- q01 + q10
  q01 / q - q01 / q * exp(- q * t)
}

stMean1 <- expected.trait(0.0145, 0.0355, 60)
stMean2 <- expected.trait(0.0245, 0.0755, 60)
stMean3 <- expected.trait(0.01, 0.015, 50)

boxplot.matrix(nullSTMean, use.cols = FALSE, outline = FALSE)
abline(h = stMean1)
abline(h = stMean2, col = 'red')
abline(h = stMean3, col = 'blue')
boxplot.matrix(traitSTMean, use.cols = FALSE, outline = FALSE, add = TRUE, col = "green")

boxplot.matrix(nullSTVar, use.cols = FALSE, outline = FALSE)
abline(h = stMean1 * (1 - stMean1))
abline(h = stMean2 * (1 - stMean2), col = 'red')
abline(h = stMean3 * (1 - stMean3), col = 'blue')
boxplot.matrix(traitSTVar, use.cols = FALSE, outline = FALSE, col = "green", add = TRUE)
