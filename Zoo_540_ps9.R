#Janelle, Hangkai, Quinn
# Possibly useful libraries
library(ape)
library(phylolm)
library(lme4)
library(lmerTest)
library(nlme)
library(mvtnorm)
library(phytools)
library(phyr)
library(tidyverse)

library(factoextra)
library(cluster)
library(dendextend)


# this is the inverse logit function that might be useful for simulations (although I'm not suggesting simulations)
inv.logit <- function(x){ 
  
  1/(1 + exp(-x))
}


# Data structure

# •16 sites, 20 plots in each (50m2)

# Response Variables

# •Shrub Abundance
# •Shrub Coverage (not included in the data sets)

# Environmental Variables

# •Light Environment
# •Gap fraction, Leaf Area Index
# •Climate 
# •Mean Annual Temperature (daily min, max, ave)
# •Mean Annual Precipitation
# •Soil Environment
# •A1 horizon depth, % sand, silt, and clay

# Central question
# How do light and edaphic conditions affect understory shrub community composition in mesic forests?

############################################################
# load and prep data
############################################################

# load the plot-level data
d <- read.csv("Zoo_540_PS9_PlotLevelData.csv")
d$plotid <- as.factor(d$plotid)
d$site <- as.factor(d$site)
summary(d)

# make a list of the names of the species
name.list <- names(d)[13:ncol(d)]

# load the phylogeny
phy <- read.tree("output_tree.tre")

# compare the number of species in d and phy
length(name.list)
Ntip(phy)

# check to see if the names match
setdiff(phy$tip.label, name.list)
setdiff(name.list, phy$tip.label)

# since there are some different names, compare them directly
cbind(name.list, sort(phy$tip.label))

# change the names in d to match the phylogeny
names(d)[names(d) == "Viburnum_rafinesquianum"] <- "Viburnum_rafinesqueanum"
names(d)[names(d) == "Lonicera_spp"] <- "Lonicera"
names(d)[names(d) == "Euonymus_spp"] <- "Euonymus"
names(d)[names(d) == "Rubus_spp"] <- "Rubus"
names(d)[names(d) == "Crataegus_spp"] <- "Crataegus"
names(d)[names(d) == "Amelanchier_spp"] <- "Amelanchier"

# reconstruct name.list
name.list <- names(d)[13:ncol(d)]
setdiff(phy$tip.label, name.list)
setdiff(name.list, phy$tip.label)

# load the site-level data that has additional environmental variables
d.summary <- read.csv("Zoo_540_PS9_SiteLevelData.csv")

# add site-level environmental variables to d
d <- merge(d, d.summary[,c("site","Abv.Dir","MinT","MaxT","DailyP","AveT")])

# create d.agg from d that aggregates the data to site
env.list <- names(d)[c(7:12, 57:61)]
d.agg <- aggregate(as.formula(paste0(env.list[1]," ~ site")), data = d, FUN = mean)
for(i in 2:length(env.list)) {
  d.agg <- cbind(d.agg, aggregate(as.formula(paste0(env.list[i]," ~ site")), data = d, FUN = mean)[,2])
  names(d.agg)[ncol(d.agg)] <- env.list[i]
}
d.agg <- merge(d.agg, d.summary[,c("site","region")])
for(i in 1:length(name.list)) {
  d.agg <- cbind(d.agg, aggregate(as.formula(paste0(name.list[i]," ~ site")), data = d, FUN = sum)[,2])
  names(d.agg)[ncol(d.agg)] <- name.list[i]
}
d.agg$region <- as.factor(d.agg$region)
summary(d.agg)

# compare d.agg with d.summary of site-level data
cbind(d.summary[,1:8], d.agg[,1:7])
# these all check except there are smallish differences in A1depth

# another way to measure community composition is the number of plots in a site containing at least one individual of a species
d.occupancy <- aggregate(as.formula(paste0(env.list[1]," ~ site")), data = d, FUN = mean)
for(i in 2:length(env.list)) {
  d.occupancy <- cbind(d.occupancy, aggregate(as.formula(paste0(env.list[i]," ~ site")), data = d, FUN = mean)[,2])
  names(d.occupancy)[ncol(d.occupancy)] <- env.list[i]
}
d.occupancy <- merge(d.occupancy, d.summary[,c("site","region")])
FUN <- function(x) return(sum(x > 0))
for(i in 1:length(name.list)) {
  d.occupancy <- cbind(d.occupancy, aggregate(as.formula(paste0(name.list[i]," ~ site")), data = d, FUN = FUN)[,2])
  names(d.occupancy)[ncol(d.occupancy)] <- name.list[i]
}
d.occupancy$region <- as.factor(d.occupancy$region)
summary(d.occupancy)

# take a look at the occupancy and abundance of the species
sp.distribution <- t(rbind(colSums(0 < d.agg[,-(1:13)]), colSums(d.occupancy[,-(1:13)]), colSums(d.agg[,-(1:13)])))
colnames(sp.distribution) <- c("sites","plots","abundance")
sp.distribution

# sort d.agg and d.occupancy so the order of species names match phy
d.agg <- d.agg[,c(1:13, 13 + match(phy$tip.label, name.list))]
d.occupancy <- d.occupancy[,c(1:13, 13 + match(phy$tip.label, name.list))]
cbind(names(d.agg)[-(1:13)], phy$tip.label)

# plots to compare mean abundance vs. proportion accupancy
par(mfrow=c(2,1), mai = c(.3,.3,.3,.3))
m.agg <- t(log(as.matrix(d.agg[,-(1:13)])+1))
image(m.agg, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "mean abundance", col=gray(.01*(1:100)))
mtext(side=1, "site", padj = .5, cex=1)
mtext(side=2, "sp", padj = -.5, cex=1)

m.occupancy <- t(as.matrix(d.occupancy[,-(1:13)]))
image(m.occupancy, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "mean occupancy", col=gray(.01*(1:100)))
mtext(side=1, "site", padj = .5, cex=1)
mtext(side=2, "sp", padj = -.5, cex=1)

# figure of d.agg with phylogeny (after Fig. 4.3 in the textbook)
layout(matrix(c(1,2,2,3,4,4), ncol=2))
par(mai=c(0,.35,.5,0))
plot(phy, direction="downwards", cex=.01, cex.main=1.5)
mtext(side=3, "mean abundance", padj = -.5, cex=1)
par(mai=c(1,.5,0,.15))
image(m.agg, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "", col=gray(.01*(1:100)))
mtext(side=1, "sp", padj = .5, cex=1)
mtext(side=2, "site", padj = -.5, cex=1)

par(mai=c(0,.35,.5,0))
plot(phy, direction="downwards", cex=.01, cex.main=1.5)
mtext(side=3, "occupancy", padj = -.5, cex=1)
par(mai=c(1,.5,0,.15))
image(m.occupancy, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "", col=gray(.01*(1:100)))
mtext(side=1, "sp", padj = .5, cex=1)
mtext(side=2, "site", padj = -.5, cex=1)

###########################################
# Create datasets with "common" species

# For analyses of community composition at the individual species level, it is often best to remove rare species, because there isn't enough information to analyze their distributions. 

# set a minimum occupancy to be included
min.occupancy <- 3

common.species.index <- (min.occupancy <= colSums(d.occupancy[,-(1:13)]))
common.name.list <- names(common.species.index)[common.species.index]
rare.name.list <- setdiff(name.list, common.name.list)
length(common.name.list)

# create versions of d.agg and d.occupancy for common species
d.agg.common <- d.agg[,!is.element(names(d.agg), rare.name.list)]
d.occupancy.common <- d.occupancy[,!is.element(names(d.occupancy), rare.name.list)]

# take a look at the occupancy and abundance of the species
sp.distribution <- t(rbind(colSums(0 < d.agg.common[,-(1:13)]), colSums(d.occupancy.common[,-(1:13)]), colSums(d.agg.common[,-(1:13)])))
colnames(sp.distribution) <- c("sites","plots","abundance")
sp.distribution

# create a reduced phylogeny to include only common species (although for pglmm() this actually isn't necessary)
phy.common <- drop.tip(phy, rare.name.list) 

# change the species order in d.agg.common and d.occupancy.common to match phy.common
d.agg.common <- d.agg.common[,c(1:13, 13 + match(phy.common$tip.label, common.name.list))]
d.occupancy.common <- d.occupancy.common[,c(1:13, 13 + match(phy.common$tip.label, common.name.list))]
cbind(names(d.agg.common)[-(1:13)], phy.common$tip.label)

# figure of d.agg.common with phylogeny (after Fig. 4.3 in the textbook)
m.agg.common <- t(log(as.matrix(d.agg.common[,-(1:13)])+1))
m.occupancy.common <- t(log(as.matrix(d.occupancy.common[,-(1:13)])+1))

layout(matrix(c(1,2,2,3,4,4), ncol=2))
par(mai=c(0,.42,.5,0.07))
plot(phy.common, direction="downwards", cex=.01, cex.main=1.5)
mtext(side=3, "mean abundance", padj = -.5, cex=1)
par(mai=c(1,.5,0,.15))
image(m.agg.common, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "", col=gray(.01*(1:100)))
mtext(side=1, "sp", padj = .5, cex=1)
mtext(side=2, "site", padj = -.5, cex=1)

par(mai=c(0,.42,.5,0.07))
plot(phy.common, direction="downwards", cex=.01, cex.main=1.5)
mtext(side=3, "occupancy", padj = -.5, cex=1)
par(mai=c(1,.5,0,.15))
image(m.occupancy.common, xaxt="n", yaxt="n", ylab = "", xlab = "", main = "", col=gray(.01*(1:100)))
mtext(side=1, "sp", padj = .5, cex=1)
mtext(side=2, "site", padj = -.5, cex=1)

# for some analyses, having data.frames in "long" format is necessary
library(reshape2)

dat <- melt(d.agg.common[,-(2:13)], id = "site")
names(dat)[2:3] <- c("sp","abundance")

dat.occ <- melt(d.occupancy.common[,-(2:13)], id = "site")
names(dat.occ)[2:3] <- c("sp","occupancy")
dat <- merge(dat, dat.occ)

# add environmental variables
dat <- merge(dat, d.agg.common[,1:13])
summary(dat)

# double check merging of files and compare abundance with occupancy
par(mfrow = c(1,1))
plot(abundance ~ occupancy, data = dat)

# finally, it is a good idea to scale independent variables
for(i.col in env.list) dat[,i.col] <- scale(dat[,i.col])

#### Exploratory Code ####
# look at data, plot level and site level richness plus calc mean
# subset plot data by site
AW <- subset(plotlevel, plotlevel$site == "AW")
BAX <- subset(plotlevel, plotlevel$site == "BAX")
BH <- subset(plotlevel, plotlevel$site == "DH")
CH <- subset(plotlevel, plotlevel$site == "CH")
DH <- subset(plotlevel, plotlevel$site == "DH")
EH <- subset(plotlevel, plotlevel$site == "EH")
KS <- subset(plotlevel, plotlevel$site == "KS")
LL <- subset(plotlevel, plotlevel$site == "LL")
MW <- subset(plotlevel, plotlevel$site == "MW")
NOW <- subset(plotlevel, plotlevel$site == "NOW")
PL <- subset(plotlevel, plotlevel$site == "PL")
PLH <- subset(plotlevel, plotlevel$site == "PLH")
RL <- subset(plotlevel, plotlevel$site == "RL")
SPD <- subset(plotlevel, plotlevel$site == "SPD")
TI <- subset(plotlevel, plotlevel$site == "TI")
WW <- subset(plotlevel, plotlevel$site == "WW")

# barplot code collapsed
#####
# AW
AWrmean <- mean(AW$richness)
print(AWrmean)
barplot(AW$richness, names.arg = AW$plotid, las = 2, cex.names = .8, main = "AW Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", AWrmean))

# BAX
BAXrmean <- mean(BAX$richness)
print(BAXrmean)
barplot(BAX$richness, names.arg = BAX$plotid, las = 2, cex.names = .7, main = "BAX Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", BAXrmean))

# BH
BHrmean <- mean(BH$richness)
print(BHrmean)
barplot(BH$richness, names.arg = BH$plotid, las = 2, cex.names = .7, main = "BH Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", BHrmean))

# CH
CHrmean <- mean(CH$richness)
print(CHrmean)
barplot(CH$richness, names.arg = CH$plotid, las = 2, cex.names = .7, main = "CH Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", CHrmean))

# DH
DHrmean <- mean(DH$richness)
print(DHrmean)
barplot(DH$richness, names.arg = DH$plotid, las = 2, cex.names = .7, main = "DH Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", DHrmean))

# EH
EHrmean <- mean(EH$richness)
print(EHrmean)
barplot(EH$richness, names.arg = EH$plotid, las = 2, cex.names = .7, main = "EH Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", EHrmean))

# KS
KSrmean <- mean(KS$richness)
print(KSrmean)
barplot(KS$richness, names.arg = KS$plotid, las = 2, cex.names = .7, main = "KS Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", KSrmean))

# LL
LLrmean <- mean(LL$richness)
print(LLrmean)
barplot(LL$richness, names.arg = LL$plotid, las = 2, cex.names = .7, main = "LL Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", LLrmean))

# MW
MWrmean <- mean(MW$richness)
print(MWrmean)
barplot(MW$richness, names.arg = MW$plotid, las = 2, cex.names = .7, main = "MW Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", MWrmean))

# NOW
NOWrmean <- mean(NOW$richness)
print(NOWrmean)
barplot(NOW$richness, names.arg = NOW$plotid, las = 2, cex.names = .7, main = "NOW Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", NOWrmean))

# PL
PLrmean <- mean(PL$richness)
print(PLrmean)
barplot(PL$richness, names.arg = PL$plotid, las = 2, cex.names = .7, main = "PL Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", PLrmean))

# PLH
PLHrmean <- mean(PLH$richness)
print(PLHrmean)
barplot(PLH$richness, names.arg = PLH$plotid, las = 2, cex.names = .7, main = "PLH Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", PLHrmean))

# RL
RLrmean <- mean(RL$richness)
print(RLrmean)
barplot(RL$richness, names.arg = RL$plotid, las = 2, cex.names = .7, main = "RL Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", RLrmean))

# SPD
SPDrmean <- mean(SPD$richness)
print(SPDrmean)
barplot(SPD$richness, names.arg = SPD$plotid, las = 2, cex.names = .7, main = "SPD Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", SPDrmean))

# TI
TIrmean <- mean(TI$richness)
print(TIrmean)
barplot(TI$richness, names.arg = TI$plotid, las = 2, cex.names = .7, main = "TI Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", TIrmean))

# WW
WWrmean <- mean(WW$richness)
print(WWrmean)
barplot(WW$richness, names.arg = WW$ploWWd, las = 2, cex.names = .7, main = "WW Plot Richness", xlab = "plot id", ylab = "richness")
legend("topright", legend = paste("Mean =", WWrmean))
#####

# plot means by site
site <- c("AW", "BAX", "BH", "CH", "DH", "EH", "KS", "LL", "MW", "NOW", "PL", "PLH", "RL", "SPD", "TI", "WW")
nordsud <- c("s", "s", "s", "n", "n", "s", "n", "n", "s", "s", "n", "n", "n", "n", "s", "s")
meanrich <- c(AWrmean, BAXrmean, BHrmean, CHrmean, DHrmean, EHrmean, KSrmean, LLrmean, MWrmean, NOWrmean, PLrmean, PLHrmean, RLrmean, SPDrmean, TIrmean, WWrmean)

mrichsite <- data.frame(site, nordsud, meanrich)

colr <- ifelse(mrichsite$nordsud == "n", "red", "blue")

barplot(mrichsite$meanrich, names.arg = mrichsite$site, las = 2, cex.names = .8, main = "Site Richness", xlab = "stie id", ylab = "mean richness", col = colr)
legend("topright", legend = unique(mrichsite$nordsud), fill = unique(colr))


summary_file <- file("model_summaries.txt", "w")

# Traverse the list of species and construct a model for each species
for(species in name.list) {
  formula <- as.formula(paste(species, "~ A1depth + clay + silt + sand + (1 | site)", collapse = " "))
  
  model <- lmer(formula, data = d)
  
  model_summary <- summary(model)
  print(paste("Model for species:", species))
  print(model_summary)
  writeLines(paste("Model for species:", species), summary_file)
  writeLines(capture.output(print(model_summary)), summary_file)
  
  # residual plots
  plot(residuals(model), main=paste("Residuals for species:", species))
  
  # visualization: Relationship between environmental factors and species richness
  ggplot(d, aes_string(x = "A1depth", y = species)) + 
    geom_point() + 
    geom_smooth(method = "lm") +
    labs(title = paste("Relationship between A1depth and", species), x = "A1depth", y = "Richness")
}

close(summary_file)

####Clustering####

#A heavy handed batching approach 

ez_dat <- dat %>%
  select(site, sp, abundance)

ez_matrix <- spread(ez_dat, key=sp, value=abundance)
#Change value to "sp" if you want some real funk 

dend<- hclust(dist(ez_matrix), method = "ward.D2")

dend1 <- as.dendrogram(dend)
plot(dend1, main = "Hierarchical Clustering of Sites based on Species Richness")
#Labels won't work nicely - but this tree should match up with the "ez_matrix" sites

#Big cluster by site 
tp_dat <- dat %>%
  select(sp, abundance, DailyP, AveT, gap_fraction, LAI, A1depth, clay, silt, sand)

big_matrix <- spread(tp_dat, key = sp, value = abundance, fill = 0)

#Standardizing variables 
standardized_matrix <- scale(big_matrix)
#Don't think we need to do the following:
#env_variables <- tp_dat %>%
#  select(AveT, DailyP,gap_fraction, LAI, A1depth, clay, silt, sand )
#standardized_env_variables <- scale(env_variables)

#combined_data <- cbind(standardized_species_matrix, standardized_env_variables)

dend2 <- hclust(dist(standardized_matrix), method = "ward.D2")
dend3 <- as.dendrogram(dend2)
plot(dend3, main = "Hierarchical Clustering of Sites based on Species Richness and Environmental Variables")
#Again - labels are not behaving here, but match this up with "ez matrix" sites 


#Goofin 

ez_dat0 <- dat %>%
  select(site, sp, abundance)

ez_matrix0 <- spread(ez_dat0, key=sp, value=abundance)

ez_matrix0 <- ez_matrix0[,-1]

ez_matrix00 <- scale(ez_matrix0)
#Change value to "sp" if you want some real funk and use ez_matrix0

dend0<- hclust(dist(ez_matrix00), method = "ward.D2")

dend00 <- as.dendrogram(dend0)
plot(dend00, main = "Hierarchical Clustering of Sites based on Species Richness")
#Labels won't work nicely - but this tree should match up with the "ez_matrix" sites
