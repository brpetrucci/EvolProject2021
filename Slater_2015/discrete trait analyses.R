#### perform analyses to determine best model for dietary strategy over posterior distribution of canid trees ####
rm(list=ls())
require(geiger);
require(phytools);
require(OUwie);
require(MASS)

load('canid_tree_data.rdata'); ## load trees 

d <-  read.csv("all.final.data.csv", stringsAsFactors =F, row.names =1) ## 

td <- treedata(best.tree, d, sort = T)

best.tree <- td$phy
d <- td$data

#### set up data #####
diet <- setNames(as.numeric(d[ ,"diet"]), rownames(d))

n <- length(diet)

diet

## run analyses on "best" tree ##

er <- fitDiscrete(best.tree, diet, model = "ER", control = list(niter=10))
sym <- fitDiscrete(best.tree, diet, model = "SYM", control = list(niter=10))
ard <- fitDiscrete(best.tree, diet, model = "ARD", control = list(niter=10))

## save AICc values ##

aic.disc <- setNames(c(er$opt$aicc, sym$opt$aicc, ard$opt$aicc), c("er", "sym", "ard"))

## compute Akaike Weights ##

aicw(aic.disc)

#### plot ASR results ####
#### Figure 1 in paper####

asr <- rerootingMethod(best.tree, diet, model = "SYM")
asr$marg

co <- c("red", "blue", "green")
co[as.factor(diet)]

quartz(width = 7.007874, height = 3.4)
par(mar=c(2,3,2,2))
tp <- plot(ladderize(best.tree), cex = 0.4, edge.width = 1, edge.col = "darkgray", direction = "upwards", tip.col = co[as.factor(diet)]);
axisPhylo(side =2, las=2, cex.axis=0.7, lwd=1, col="black", col.axis="black")
mtext("Millions of Years Ago", 2, line = 2, at = 20, cex=1, col="black")

text(24.5, 61.5, label ="Borophaginae", cex=0.6, font=2, col ="black")
segments(0, 60.5, 48, 60.5, lwd = 1, col = "black")
text(60, 61.5, label ="Hesperocyoninae", cex=0.6, font=2, col ="black")
segments(49, 60.5, 69, 60.5, lwd = 1, col = "black")
text(83, 61.5, label ="Caninae", cex=0.6, font=2, col ="black")
segments(70, 60.5, 93, 60.5, lwd =1, col = "black")

nd <- max(diag(vcv(best.tree)))

rect(94, nd, 96, nd-2.588, col = rgb(red=255,green =242,blue=174,alpha = 255 , maxColorValue = 255) , border = "black") # pleistocene
rect(94,  nd-2.588, 96, nd-5.332,col = rgb(red=255,green =255,blue=153,alpha = 255 , maxColorValue = 255) , border = "black")# pliocene
rect(94,  nd-5.332, 96, nd-23.03, col = rgb(red=255,green =255,blue=0,alpha = 255 , maxColorValue = 255) , border = "black") # miocene
rect(94,  nd-23.03, 96, nd-33.9, col = rgb(red=253,green =192,blue=122,alpha = 255 , maxColorValue = 255) , border = "black") #oligocene
rect(94,  nd-33.9, 96, 0,col =  rgb(red=253,green =180,blue=108,alpha = 255 , maxColorValue = 255) , border = "black")


nodelabels(pie = asr$marg, piecol = co, cex = 0.4)
legend("bottomleft", legend =c("hypercarnivore", "mesocarnivore", "hypocarnivore"), pch=21, pt.bg=c("red", "blue", "green"),pt.cex=1.5,cex=0.8, box.lwd=0, bty ="o",bg="white")


###### run analyses over posterior sample ######
###### then compute model averaged rates  ######

## results matrices ##
er.res <- sym.res <- ard.res<- matrix(NA, nrow = 500, ncol = 8, dimnames = list(NULL, c("q12", "q13","q21","q23","q31","q32","lnl", "aicc")))
weights <- matrix(NA, nrow = 500, ncol =3, dimnames = list(NULL, c("er", "sym", "ard")))

##### LOOP THROUGH THE TREES HERE ######
##### use only run 1 from the     ######
##### MrBayes analysis. The runs  ######
##### converged, so this is fine  ######

tree.sample <-read.csv("posteriorSample.csv")[,1]

for(i in 1:length(tree.sample)) {
	
	## sample a tree, scale the edge lengths to time in my
	p <- t.file1[[tree.sample[i]]]
	p$edge.length <-  p$edge.length / p.file1[tree.sample[i],"Clockrate"]
	phy.tmp <- treedata(p, diet)$phy ## tree is ready
	
	# fit models
	

	er <- fitDiscrete(phy.tmp, diet, model = "ER", control = list(niter=10))
	er.res[i, ] <- as.numeric(unlist(er$opt)[c(1,2,3,4,5,6,7,10)])

	sym <- fitDiscrete(phy.tmp, diet, model = "SYM", control = list(niter=10))
	sym.res[i, ] <- as.numeric(unlist(sym$opt)[c(1,2,3,4,5,6,7,10)])

	ard <- fitDiscrete(phy.tmp, diet, model = "ARD", control = list(niter=10))
	ard.res[i, ] <- as.numeric(unlist(ard$opt)[c(1,2,3,4,5,6,7,10)])
	

	
	aic <- setNames(c(er$opt$aicc,sym$opt$aicc, ard$opt$aicc ), c("er", "sym",  "ard"))
	weights[i, ]<- aicw(aic)[,3];
	
	## progress counter ##
	if(i%%50 == 0) {
		print(paste("Samples complete = ", i))
	}	
}


## create output files ##
write.csv(er.res, "equal_rates_500.csv")
write.csv(sym.res, "sym_rates_500.csv")
write.csv(ard.res, "diff_rates_500.csv")
write.csv(weights, "discrete_weights_500.csv")

## get median param estimates ##
median.params <- round(rbind(apply(er.res,2, median), apply(sym.res,2, median), apply(ard.res,2, median)),4)
write.csv(cbind(median.params, apply(weights, 2, median)), "median.res_500.csv")

## model average the rates 

	er <- er.res
	sym <- sym.res
	ard <- ard.res
	model.av.weights <- function(er, sym, ard, weights) {
	
	er.rates <-er[,1:6]
	sym.rates <-sym[,1:6]
	ard.rates <-ard[,1:6]
	
	weighted.er <- er.rates * weights[,"er"]
	weighted.sym <- sym.rates * weights[,"sym"]
	weighted.ard <- ard.rates * weights[,"ard"]
	
	mars <- matrix(NA, nrow(er), 6)
	
	for(i in 1:ncol(mars)) {
		mars[,i] <- (apply(cbind(weighted.er[,i], weighted.sym[,i], weighted.ard[,i]), 1, sum)) 
	}
	colnames(mars) <- colnames(weighted.ard)
	return(mars)
	
}

maw <- model.av.weights(er, sym, ard, weights)

## get median rates over 500 model averaged samples ##
maw.median <- apply(maw, 2, median)


### plot transition rates ###

require(diagram);


M <- c(0, "q12", "q13", "q21", 0, "q23", "q31", "q32", 0 )
M <- matrix(data=M, nrow=3, ncol=3, byrow=T)


A <- maw.median[match(M, names(maw.median))]
A[is.na(A)] <- 0
A <- matrix(data = A, nrow = 3, ncol = 3,byrow = T)
A <- A *100
#A[which(A<1 & A>0)] <- 0.01


quartz(width= 3.4, height= 3.4)
par(mar=c(0,0,0,0))


plot.line.col="black"
pp <- plotmat(M, pos = c(1, 2), name = c("hypercarnivore", "mesocarnivore", "hypocarnivore"),lwd = 1, lcol = plot.line.col, arr.lwd = A, arr.col = plot.line.col, box.lwd = 2,cex.txt = 0, box.size = 0.2, box.type = "ellipse", box.prop = 0.5, box.col = c("red",  "blue","green"), box.cex=1)


