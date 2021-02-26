### discriminant analysis multivariable ###
rm(list=ls())

## read in extant canid and procyonid data ##

ed <- read.csv("extant_ratios.csv", row.names =1)

ed[,"mass"] <- log(ed[,"mass"])


require(MASS)

## classify as many species as possible using lower teeth and moment arms only. 

disc<- lda(diet~ RBL + M1BS + M2S + MAT +p4S, data = as.data.frame(ed))
dis.pred <- predict(disc)

## look at posterior classifications
cbind(rownames(ed), ed[,"diet"], dis.pred$class)
plot(dis.pred$x, pch = 21, bg = dis.pred$class )
text(dis.pred$x, label = rownames(dis.pred$x), cex=0.4)

write.csv(coef(disc), "discriminant_1_coefficients.csv")

### fossil data

fd <- read.csv("fossil.ratios.csv", row.names =1)
na.foo  <- function (x) sum(!is.na(x)) 

fd.pred <- predict(disc, as.data.frame(fd))
fossil.classif<- cbind(fd.pred$class, fd.pred$posterior, fd.pred$x)
fossil.classif[-which(is.na(fossil.classif[,1])),1]

points(fd.pred$x, pch=22, bg = fd.pred$class)
text(fd.pred$x, rownames(fd.pred$x), cex=0.4)

fossil.classif1 <- fossil.classif[-which(is.na(fossil.classif[,1])),]


write.csv(fossil.classif1, "fossil_classification1.csv")
unclass <- which(is.na(fossil.classif[,1]))
apply(fd[unclass,], 2, na.foo)

### reduced dfa - take out low representation lever arms

disc <- lda(diet~ RBL + M1BS + M2S + p4S, data = as.data.frame(ed))

disc
write.csv(coef(disc), "discriminant_2_coefficients.csv")

dis.pred <- predict(disc)
## look at posterior classifications
cbind(rownames(ed), ed[,"diet"], dis.pred$class)

plot(dis.pred$x, pch = 21, bg = dis.pred$class )
text(dis.pred$x, label = rownames(dis.pred$x), cex=0.4)

# ok, looks good. try fossils still to be classified

fd.pred2 <- predict(disc, as.data.frame(fd[unclass,]))
fossil.classif<- cbind(fd.pred2 $class, fd.pred2$posterior, fd.pred2$x)

fossil.class.2 <- fossil.classif[-which(is.na(fossil.classif[,1])),]

write.csv(fossil.class.2, "fossil_classification2.csv")
unclass <- which(is.na(fossil.classif[,1]))
apply(fd[match(names(unclass), rownames(fd)),], 2, na.foo)


### try with further reduced function - remove M1BS, 

disc<- lda(diet~ RBL + M2S + p4S, data = as.data.frame(ed))
disc
write.csv(coef(disc), "discriminant_3_coefficients.csv")

dis.pred <- predict(disc)

## look at posterior classifications
cbind(rownames(ed), ed[,"diet"], dis.pred$class)

plot(dis.pred$x, pch = 21, bg = dis.pred$class )
text(dis.pred$x, label = rownames(dis.pred$x), cex=0.4)

# still looks ok. try fossils still to be classified

dis.pred <- predict(disc, as.data.frame(fd[match(names(unclass), rownames(fd)),]))
fossil.classif<- cbind(dis.pred$class, dis.pred$posterior, dis.pred$x)

fossil.class.3 <- fossil.classif[-which(is.na(fossil.classif[,1])),]
write.csv(fossil.class.3, "fossil_classification3.csv")

unclass <- which(is.na(fossil.classif[,1]))
apply(fd[match(names(unclass), rownames(fd)),], 2, na.foo)

### the remaining taxa lack many lower tooth characters but may be present in the tree and body mass data sets.
###

# target Urocyon galushi and Enhydrocyon stenocephalus

## Urocyon

fd["Urocyon_galushi",]
disc<- lda(diet~ M2S + MAM + P4P + P3S, data = as.data.frame(ed))

disc
write.csv(coef(disc), "fossil_classification4")
dis.pred <- predict(disc)

## look at posterior classifications
cbind(rownames(ed), ed[,"diet"], dis.pred$class)
plot(dis.pred$x, pch=21, bg=dis.pred$class)

dis.pred <- predict(disc, as.data.frame(fd[match(names(unclass), rownames(fd)),]))
fossil.classif<- cbind(dis.pred$class, dis.pred$posterior, dis.pred$x)
points(dis.pred$x)
fossil.class.4 <- fossil.classif[-which(is.na(fossil.classif[,1])),]

## Now Enhydrocyon
disc<- lda(diet~ RBL + P3S, data = as.data.frame(ed))
disc
write.csv(coef(disc), "fossil_classification5")
dis.pred <- predict(disc)

## look at posterior classifications
cbind(rownames(ed), ed[,"diet"], dis.pred$class)
plot(dis.pred$x, pch=21, bg=dis.pred$class)

dis.pred <- predict(disc, as.data.frame(fd[match(names(unclass), rownames(fd)),]))
fossil.classif<- cbind(dis.pred$class, dis.pred$posterior, dis.pred$x)
points(dis.pred$x)
fossil.class.5 <- fossil.classif[-which(is.na(fossil.classif[,1])),]

fossil.class.4 <- rbind(fossil.class.4, fossil.class.5["Enhydrocyon_stenocephalus",])
rownames(fossil.class.4) <- c("Urocyon_galushi", "Enhydrocyon_stenocephalus")


write.csv(fossil.class.4, "fossil_classification4.csv")

all.disc.results <- rbind(fossil.classif1[,(1:5)], fossil.class.2[,(1:5)], fossil.class.3[,(1:5)], fossil.class.4[,(1:5)])

colnames(all.disc.results) <- c("diet", "P.hyper", "P.meso", "P.hypo", "P.herb")

write.csv(all.disc.results, "discriminant_Results.csv")

## now pull everything together and extract North American extant canids to get a complete dataset ##

diet <- all.disc.results[match(rownames(fd), rownames(all.disc.results)),1]

final.fossil.data <- cbind(diet, fd)[-which(is.na(diet)),c("diet", "RLGA", "mass")]

rownames(ed)
NA.extant <- ed[c("Canis_latrans", "Canis_lupus", "Cuon_alpinus", "Urocyon_cinereoargenteus", "Urocyon_littoralis", "Vulpes_lagopus", "Vulpes_velox", "Vulpes_vulpes"),c("diet", "RLGA", "mass")]

rownames(NA.extant)[3] <- "Cuon_javanicus"

all.final.data <- rbind(final.fossil.data, NA.extant)
all.final.data$diet[which(all.final.data$diet==4)]<-3

write.csv(all.final.data, "all.final.data.csv")