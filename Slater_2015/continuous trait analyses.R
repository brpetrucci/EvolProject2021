###### run continuous trait analyses over posterior sample ######
rm(list=ls())
require(geiger);
require(phytools);
require(OUwie);
require(MASS)

load('canid_tree_data.rdata'); ## load trees and posterior param estimates

## a function to compute AIC and AICc, which are not reported for phytools diversity models ##

aic.comp <- function(lnl, k, n = NULL, small.sample = F) {
	
	AIC <- -2*lnl + 2*k
	if(small.sample ==  F){
		return(AIC)
	} else{
		if(is.null(n)==T) stop("you must specify the sample size to use AIC.c")
		weight <- (2*k * (k+1)) / (n-k-1)
		return(AIC + weight)	
	}
} 

## set up results matrices ##

bm.model.res.mass <- bm.model.res.rlga <- bm.var.res.mass <-bm.var.res.rlga <-  acdc.model.res.mass <- acdc.model.res.rlga <- div.model.res.mass <-div.model.res.rlga <- div.diet.model.res.mass <- div.diet.model.res.rlga <- trend.model.res.mass <- trend.model.res.rlga <-  ou.model.res.mass <- ou.model.res.rlga <- matrix(data =NA, nrow = 500, ncol = 5)

colnames(bm.model.res.mass) <- colnames(bm.var.res.mass) <-  colnames(acdc.model.res.mass) <- colnames(div.model.res.mass) <- colnames(div.diet.model.res.mass) <- colnames(trend.model.res.mass) <- colnames(ou.model.res.mass) <- colnames(bm.model.res.rlga) <- colnames(bm.var.res.rlga) <-  colnames(acdc.model.res.rlga) <- colnames(div.model.res.rlga) <- colnames(div.diet.model.res.rlga) <- colnames(trend.model.res.rlga) <- colnames(ou.model.res.rlga) <-c("sample", "LnL", "aicc", "Sig2", "param")


ou.optima.mass <- ou.optima.rlga <- matrix(data =NA, nrow = 500, ncol = 4)
colnames(ou.optima.mass) <- colnames(ou.optima.rlga) <- c("tree", "hypercarnivore", "mesocarnivore", "hypocarnivore")

bm.var.rates.mass <- bm.var.rates.rlga <- matrix(data =NA, nrow = 500, ncol = 4)
colnames(bm.var.rates.mass) <- colnames(bm.var.rates.rlga) <- c("tree", "hypercarnivore", "mesocarnivore", "hypocarnivore")


weights.mass <- weights.rlga <- matrix(data =NA, nrow = 500, ncol = 7)
colnames(weights.mass) <- colnames(weights.rlga) <-c("bm", "bm.var", "acdc", "diversity", "dietary_div", "trend", "OU")

## prepare body mass and rlga data / dietary data using "best tree" ##

d <- read.csv("all.final.data.csv", stringsAsFactors =F, row.names =1) ## 

d <- treedata(best.tree, d, sort =TRUE)$data

## create two sets of data; one for mass and one for RLGA ##

# mass
diet.mass <- setNames(d[,"diet"], rownames(d));
mass <- setNames(d[,"mass"], rownames(d))
n.mass <- length(mass)

# rlga

rlga.drop <- which(is.na(d[,"RLGA"]))
d.rlga <- d[-rlga.drop,]
td <- treedata(best.tree, d.rlga)
rlga <- setNames(td$data[,"RLGA"], rownames(td$data))
rlga.diet <- setNames(td$data[,"diet"], rownames(td$data))
n.rlga <- length(rlga)

## set working directory for output files to go to ##

## setwd("wherever you want them to go")

for(i in 1:500) {
	
	### if doing this from scratch, you'll want to ###
	### sample trees from the post-burn-in sample  ###
	### there are 50k trees and we'll discared the ###
	### first 25% as burnin    					   ###
	
	## sample a tree, scale the edge lengths to time in myr ##
	# samples <-  sample(seq(5000, length(t.file1)), size=500, replace = F);
	
	### however, if you want to replicate analyses     ###
	### as they are in the paper, use the same samples ###
	
	samples <- read.csv("~/Desktop/CanidsForDryad/posteriorSample.csv")[,1];
	
	######################################################
	
	tree <- samples[i]
	phy.tmp <- t.file1[[tree]];
	phy.tmp$edge.length <-  phy.tmp$edge.length / p.file1[tree,"Clockrate"]
	

	phy.mass <- treedata(phy.tmp, mass)$phy;
	phy.rlga <- treedata(phy.tmp, rlga)$phy;
	## trees are ready	
	
	est.div.mass <- estDiversity(phy.mass, diet.mass, method = "asr", model ="SYM"); ## estimate diversity
	est.div.rlga <- estDiversity(phy.rlga, rlga.diet, method = "asr", model ="SYM"); ## estimate diversity

	
	### now fit models
	
	## MASS fi;rst, then rlga
	
	# Bm 
	bm.mass<- fitContinuous(phy.mass, mass, model="BM", control = list(niter= 10))$opt
	bm.model.res.mass[i, ] <- c(tree, bm.mass $lnL, bm.mass $aicc, bm.mass $sigsq, "NA")

	bm.rlga<- fitContinuous(phy.rlga, rlga, model="BM", control = list(niter= 10))$opt
	bm.model.res.rlga[i, ] <- c(tree, bm.rlga $lnL, bm.rlga $aicc, bm.rlga $sigsq, "NA")
	
	# ACDC 
	acdc.mass <- fitContinuous(phy.mass, mass, model="EB", bounds = list(a = c(-0.2, 0.2)), control = list(niter= 10))$opt
	acdc.model.res.mass[i, ] <- c(tree, acdc.mass $lnL, acdc.mass $aicc, acdc.mass $sigsq, acdc.mass$a)

	acdc.rlga <- fitContinuous(phy.rlga, rlga, model="EB", bounds = list(a = c(-0.2, 0.2)), control = list(niter= 10))$opt
	acdc.model.res.rlga[i, ] <- c(tree, acdc.rlga $lnL, acdc.rlga $aicc, acdc.rlga $sigsq, acdc.rlga $a)


	# Trend
	trend.mass <- fitContinuous(phy.mass, mass, model="drift", control = list(niter= 10))$opt
	trend.model.res.mass[i, ] <- c(tree, trend.mass$lnL, trend.mass$aicc, trend.mass$sigsq, trend.mass$drift)

	trend.rlga <- fitContinuous(phy.rlga, rlga, model="drift", control = list(niter= 10))$opt
	trend.model.res.rlga[i, ] <- c(tree, trend.rlga $lnL, trend.rlga $aicc, trend.rlga $sigsq, trend.rlga $drift)

	# diversity NOT by diet
	div.mass <- fitDiversityModel(phy.mass, mass, showTree =F)
	div.mass$aicc <- aic.comp(div.mass $logL, k = 3, n  = n.mass, small.sample = T)
	div.model.res.mass[i, ] <- c(tree, div.mass $logL, div.mass $aicc, div.mass $sig0, div.mass $psi)

	div.rlga <- fitDiversityModel(phy.rlga, rlga, showTree =F)
	div.rlga$aicc <- aic.comp(div.rlga $logL, k = 3, n  = n.rlga, small.sample = T)
	div.model.res.rlga[i, ] <- c(tree, div.rlga $logL, div.rlga $aicc, div.rlga $sig0, div.rlga $psi)


	# diversity by diet
	divdiet.mass <- fitDiversityModel(phy.mass, mass, d = est.div.mass, showTree =F)
	divdiet.mass$aicc <- aic.comp(divdiet.mass $logL, k = 3, n  = n.mass, small.sample = T)
	div.diet.model.res.mass[i, ] <- c(tree, divdiet.mass $logL, divdiet.mass $aicc, divdiet.mass $sig0, divdiet.mass $psi)

	divdiet.rlga <- fitDiversityModel(phy.rlga, rlga, d = est.div.rlga, showTree =F)
	divdiet.rlga$aicc <- aic.comp(divdiet.rlga$logL, k = 3, n  = n.rlga, small.sample = T)
	div.diet.model.res.rlga[i, ] <- c(tree, divdiet.rlga $logL, divdiet.rlga $aicc, divdiet.rlga $sig0, divdiet.rlga $psi)

	## set up data for OUwie and perform ancestral state estimation ##
	
	asr.mass <- rerootingMethod(phy.mass, diet.mass, type = "discrete", model = "SYM")
	asr.rlga <- rerootingMethod(phy.rlga, rlga.diet, type = "discrete", model = "SYM")
	
	best<- function(x) return(order(x, decreasing = T)[1])
	node.states.mass <- apply(asr.mass$marginal.anc,1,best)
	names(node.states.mass) <- seq(n.mass+1, n.mass+ phy.mass$Nnode)
	
	node.states.rlga <- apply(asr.rlga$marginal.anc,1,best)
	names(node.states.rlga) <- seq(n.rlga+1, n.rlga+ phy.rlga$Nnode)
	

	ouwie.mass.phy <- phy.mass; ouwie.mass.phy$node.label <- node.states.mass
	ouwie.mass.data <-data.frame(names(diet.mass), diet.mass, mass)
	
	ouwie.rlga.phy <- phy.rlga; ouwie.rlga.phy$node.label <- node.states.rlga
	ouwie.rlga.data <-data.frame(names(rlga.diet), rlga.diet, rlga)
	
	## perform OUwie analyses
	ou.mass <- OUwie(ouwie.mass.phy, ouwie.mass.data, model = "OUM", root.station =T, lb = 1e-50, ub = 10, quiet = T)
	ou.model.res.mass[i, ] <- c(tree, ou.mass $loglik, ou.mass $AICc, ou.mass $ solution[2,1], ou.mass $ solution[1,1])
	ou.optima.mass[i,]  <- c(tree, ou.mass$theta[,1])
	
	
	ou.rlga <- OUwie(ouwie.rlga.phy, ouwie.rlga.data, model = "OUM", root.station =T, lb = 1e-50, ub = 10, quiet = T)
	ou.model.res.rlga[i, ] <- c(tree, ou.rlga $loglik, ou.rlga $AICc, ou.rlga $ solution[2,1], ou.rlga $ solution[1,1])
	ou.optima.rlga[i,]  <- c(tree, ou.rlga$theta[,1])

	
	# variable rates #
	bm.var.mass <- OUwie(ouwie.mass.phy,  ouwie.mass.data, model = "BMS", quiet = T)
	bm.var.res.mass[i,] <- c(tree, bm.var.mass $loglik, bm.var.mass $AICc, NA, bm.var.mass $theta[1])
	bm.var.rates.mass[i,]  <- c(tree, bm.var.mass $solution[2,])

	bm.var.rlga <- OUwie(ouwie.rlga.phy, ouwie.rlga.data, model = "BMS", quiet = T)
	bm.var.res.rlga[i,] <- c(tree, bm.var.rlga $loglik, bm.var.rlga $AICc, NA, bm.var.rlga $theta[1])
	bm.var.rates.rlga[i,]  <- c(tree, bm.var.rlga $solution[2,])
	
	#####################
	## End of analyses ##	
	#####################
	
	
	aic.mass <- setNames(c(bm.mass$aicc, bm.var.mass $AICc, acdc.mass$aicc, div.mass$aicc, divdiet.mass$aicc,trend.mass$aicc, ou.mass$AICc), c("bm", "bm.var",  "acdc", "diversity", "dietdiversity", "trend",  "ou"))
	weights.mass[i, ]<- aicw(aic.mass)[,3];
	
	aic.rlga <- setNames(c(bm.rlga$aicc, bm.var.rlga $AICc, acdc.rlga$aicc, div.rlga$aicc, divdiet.rlga$aicc,trend.rlga$aicc, ou.rlga$AICc), c("bm", "bm.var",  "acdc", "diversity", "dietdiversity", "trend",  "ou"))
	weights.rlga[i, ]<- aicw(aic.rlga)[,3];
	
	if(i%%10 == 0) {
		print(paste("Samples complete = ", i))
	}	
}

### processing posterior sample output ###

## mass first

res.mass <- list(bm= bm.model.res.mass, acdc= acdc.model.res.mass, trend = trend.model.res.mass, bmv= bm.var.rates.mass, div= div.model.res.mass, divdiet = div.diet.model.res.mass, oum= ou.model.res.mass)

foo <- function(x) return(apply(x, 2, median))
mass.res<- matrix(unlist(lapply(res.mass, foo)), nrow = length(res.mass), ncol = 6, byrow=T)[,-1]
colnames(mass.res) <- c("loglk", "AICc", "sigmasq", "PARAM")
rownames(mass.res) <- names(res.mass)


## compute ou and bmv diet specific characters ##


diet.spec.mass<- rbind(apply(bm.var.rates.mass,2,median), apply(ou.optima.mass, 2, median))[ , -1]
rownames(diet.spec.mass) <- c("bmv", "oum")

## note that optimal masses are in Ln(Kgs) ##
 
 exp(diet.spec.mass[2,]) ## gives optimal in Kgs ##
 
## Now we repeat for RLGA

res.rlga <- list(bm= bm.model.res.rlga, acdc= acdc.model.res.rlga, trend = trend.model.res.rlga, bmv= bm.var.rates.rlga, div= div.model.res.rlga, divdiet = div.diet.model.res.rlga, oum= ou.model.res.rlga)

foo <- function(x) return(apply(x, 2, median))
rlga.res<- matrix(unlist(lapply(res.rlga, foo)), nrow = length(res.rlga), ncol = 6, byrow=T)[,-1]
colnames(rlga.res) <- c("loglk", "AICc", "sigmasq", "PARAM")
rownames(rlga.res) <- names(res.rlga)


## compute ou and bmv diet specific characters ##


diet.spec.rlga<- rbind(apply(bm.var.rates.rlga,2,median), apply(ou.optima.rlga, 2, median))[ , -1]
rownames(diet.spec.rlga) <- c("bmv", "oum")


##### End of Analysis and Script ######


