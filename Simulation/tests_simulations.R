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

# directory where all simulations are
simDir <- "/Users/petrucci/Documents/research/EvolProject2021/Simulation/replicates/"

# create objects to hold our diagnostics

# matrix to hold total number of species for null and trait sims
nullN <- traitN <- matrix(0, nrow = 13, ncol = 50)

# and same for extant number
nullExtantN <- traitExtantN <- matrix(0, nrow = 13, ncol = 50)

# rho for null and trait too
nullRho <- traitRho <- matrix(0, nrow = 13, ncol = 50)

# bm trait variance at the present divided by tMax
bmVar <- matrix(0, nrow = 13, ncol = 50)

# and for discrete trait ratio of # of 1s/# of 0s
stRatio <- matrix(0, nrow = 13, ncol = 50)

# loop through combinations of parameters
for (i in 1:13) {
  # directory for the data
  combDir <- paste0(simDir, "comb_", i)
  
  # null directory
  nullDir <- paste0(combDir, "/null/")
  
  # and traits
  traitsDir <- paste0(combDir, "/traits/")
  
  # load sim lists
  load(paste0(nullDir, "sim_list.RData"))
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
  
  # get function lists
  load(paste0(traitsDir, "bm_traits.RData"))
  load(paste0(traitsDir, "st_traits.RData"))
  
  # get bm trait variances
  bmVar[i, ] <- unlist(lapply(1:50, function(x) {
    # list of trait values at the time of simulation's end
    traitList <- as.numeric(unlist(lapply(which(simList[[x]]$EXTANT), function(y) 
      bmTraitsFunc[[x]][[paste0("t", y)]](max(simList[[x]]$TS)))))
    
    if (length(traitList) < 5) return(NA)
    
    names(traitList) <- paste0("t", which(simList[[x]]$EXTANT))
    
    # phylogeny to get MRCA time
    phy <- drop.fossil(make.phylo(simList[[x]]))
    phy$root.time <- NULL
    
    # estimate BM variance
    return(fitContinuous(phy, traitList)$opt$sigsq)
  }))
  
  # and the discrete ratios
  stRatio[i, ] <- unlist(lapply(1:50, function(x) {
    traitList <- as.numeric(
      unlist(
        read.nexus.data(
          paste0(traitsDir, "trait_lists/st_mol_traits_", x, ".nex"))))
    sum(traitList) / (length(traitList) - sum(traitList))
  }))
}




