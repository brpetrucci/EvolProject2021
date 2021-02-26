# get the file name for the tree
dir <- "C:/Users/bruno/Documents/RESEARCH/PhD/EvolutionProject/Slater_2015/"
tree_filename <- paste0(dir, "canidae_dated_MCCT.tre")

tree <- ape::read.nexus(file = tree_filename)

extant_data <- read.csv("C:/Users/bruno/Documents/RESEARCH/PhD/EvolutionProject/Slater_2015/extant_ratios.csv")
extinct_data <- read.csv("C:/Users/bruno/Documents/RESEARCH/PhD/EvolutionProject/Slater_2015/fossil.ratios.csv")

taxon_list <- tree$tip.label

final_data <- read.csv(paste0(dir, "all.final.data.csv"))
final_data[, 1] <- as.character(final_data[, 1])
final_taxon <- final_data[, 1]

final_list <- final_taxon[final_taxon %in% taxon_list]

mass_data <- final_data[final_data[, 1] %in% final_list, c("X", "mass")]
diet_data <- final_data[final_data[, 1] %in% final_list, c("X", "diet")]
  
new_tree <- drop.tip(tree, which(!(taxon_list %in% final_list)))
write.nexus(new_tree, file = paste0(dir, "calibration_tree.tre"))

nexus_mass <- mass_data[, 2]
names(nexus_mass) <- mass_data[, 1]

nexus_diet <- diet_data[, 2]
names(nexus_diet) <- diet_data[, 1]

write.nexus.data(x = nexus_mass, file = paste0(dir, "mass_data.nex"))
write.nexus.data(x = nexus_diet, file = paste0(dir, "diet_data.nex"))

