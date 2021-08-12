################################################################################
# 08.2021 - Script written for volume and shape analyses of turtle endocasts
# Published in Ferreira GS, Werneburg I, Lautenschlager S, Evers SW. Contrasting 
#   brains and bones: neuroanatomical evolution of turtles (Testudinata). 
#       In: Dozo MT, Paulina-Carabajal A, Macrini TE, Walsh S. "Paleoneurology
#           of Amniotes: new directions in the study of fossil endocasts"
#
#
# What this script does:
# 1. Volume analyses:
#       1a. PGLS model of endocast volume ~ box volume
#		1b. Ancestral states estimate of endocast volumen
# 		1c. PGLS model: brain volume ~ endocast volume
#
# 2. Shape Analyses with linear measurements
#
# 3. Semilandmark (outline) analyses
#		3a. PCA complete dataset
#		3b. PCA adults only dataset
#
################################################################################

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(paleotree)
library(rcompanion)
library(geomorph)
library(rr2)

getwd()

# Load image
load(file = "results/Ferreiraetal.RData")

##############################################################
# 1. Volume Analyses
##############################################################

# Load Farina et al. in press pruned tree
vol.MCCT = read.tree("Data/Farina_pruned_tree.tre")
plot(ladderize(vol.MCCT, right = F))  ## check the tree
axisPhylo()

# Load endocast volume raw data
vol.data = read.table("Data/volume_data_clades.txt", header = T, row.names = 1)

# check if names on the tree and data are the same
name.check(vol.MCCT, vol.data)  ## check if names are the same

# extract columns
endocast.vol = vol.data[, "endocast_vol"]
box.vol = vol.data[, "box_vol"]
eb.vol = vol.data[, "eb_vol"]
clades = vol.data$large_clades

# transform data
endocastLog = log(endocast.vol)
boxLog = log(box.vol)
eb.volLog = log(eb.vol)

# apply names
names(endocastLog) = names(boxLog) = names(eb.volLog) = rownames(vol.data)
names(clades) = rownames(vol.data)



################################################################################
# 1a. PGLS model: endocast volume ~ box volume
################################################################################

# using Pagel's Lambda
tip.heights = diag(vcv(phy = vol.MCCT))  ## a weights argument is necessary when
                                ## the tree include fossils as tips

pagelModel = gls(endocastLog ~ boxLog, correlation = corPagel(1, phy = vol.MCCT),
                 data = vol.data, weights = varFixed(~tip.heights))
summary(pagelModel)

# Calculate R2 statistics for the pglsModel
pglsModel.r2 = R2.pred(pagelModel, phy = vol.MCCT)
pglsModel.r2

# plot log-transformed data with pglsModel line
plot(endocastLog ~ boxLog, xlab = "Log-transformed Box volume", 
     ylab = "Log-transformed Endocast volume")
abline(coef(pagelModel))

# plot untransformed data for comparison
plot(endocast.vol ~ box.vol, xlab = "Box volume", 
     ylab = "Endocast volume")
abline(lm(endocast.vol ~ box.vol))



################################################################################
# 1b. Ancestral states estimate of endocast volume
################################################################################

fit.EBvol = fastAnc(vol.MCCT, eb.volLog, vars = TRUE, CI = TRUE)
range(eb.volLog)

# plot log-transformed endocast/box volume on the tree
obj = contMap(vol.MCCT, eb.volLog, plot = F)
obj = setMap(obj, colors=c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677"))
plot(obj, main = "Trait", legend=0.5*max(nodeHeights(vol.MCCT)), fsize=c(0.7,0.9),
     leg.txt = "log(Endocast volume)")

# plot in a pdf
pdf("results/ancestral-state-volumes.pdf")

plot(obj, main = "Trait", legend=0.5*max(nodeHeights(vol.MCCT)), fsize=c(0.7,0.9),
     leg.txt = "log(Endocast volume)")

dev.off()


################################################################################
# 1c. PGLS model: brain volume ~ endocast volume
################################################################################

BrainEndoVol = read.csv("Data/brain_endocast_vol.csv", header = T, sep = ",")

BrainEndoVol.log = BrainEndoVol
BrainEndoVol.log$endocast = log(BrainEndoVol$endocast)
BrainEndoVol.log$brain = log(BrainEndoVol$brain)

# load tree with taxa that contain brain data
dev.off()

brain.tree = read.tree("Data/brain_reduced_tree.tre")
plot(brain.tree)
axisPhylo()

# pgls using Pagel's Lambda
tip.heights.brain = diag(vcv(phy = brain.tree))  ## a weights argument is 
                    ## necessary the tree include fossils as tips

pagelModel.brain = gls(brain ~ endocast, correlation = corPagel(1, phy = brain.tree),
                 data = BrainEndoVol.log, weights = varFixed(~tip.heights.brain))
summary(pagelModel.brain)

# Calculate R2 statistics for the pglsModel
pglsModel.r2 = R2.pred(pagelModel.brain, phy = brain.tree)
pglsModel.r2

# plot brain volume ~ endocast volume
plot(BrainEndoVol$brain ~ BrainEndoVol$endocast, bty = "l", pch = 18,
     xlab = "log-transformed Endocast volume", ylab = "log-transformed Brain volume")
title("A", adj = 0)
legend("topleft", paste("y = ", round(res.BrEn$coefficients[1], digits = 3), 
                        " + ", round(res.BrEn$coefficients[2], digits = 3), "
R^2 = ", round(res.BrEn$r.squared, digits = 3)), bty = "n")
abline(reg.BrainEndo, col = "#AA3377", lwd = 1.5)



################################################################################
# 2. Shape Analyses with linear measurements
################################################################################

endo_meas = read.table("Data/endocast_measurements.txt", header = T, row.names = 1, sep = " ")
brain_meas = read.table("Data/brain_measurements.txt", header = T, row.names = 1, sep = " ")

label.names = colnames(endo_meas)

WOB = data.frame(brain_meas$WOB/endo_meas$WOB, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(WOB)
names(WOB) = c("be", "juv", "id")

WCH = data.frame(brain_meas$WCH/endo_meas$WCH, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(WCH)
names(WCH) = c("be", "juv", "id")

WOR = data.frame(brain_meas$WOR/endo_meas$WOR, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(WOR)
names(WOR) = c("be", "juv", "id")

HCH = data.frame(brain_meas$HCH/endo_meas$HCH, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(HCH)
names(HCH) = c("be", "juv", "id")

HOB = data.frame(brain_meas$HOB/endo_meas$HOB, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(HOB)
names(HOB) = c("be", "juv", "id")

HOR = data.frame(brain_meas$HOR/endo_meas$HOR, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(HOR)
names(HOR) = c("be", "juv", "id")

BE = data.frame(BrainEndoVol$brain/BrainEndoVol$endocast, endo_meas$juv,endo_meas$id, row.names = rownames(endo_meas))
colnames(BE)
names(BE) = c("be", "juv", "id")

# plot to pdf
pdf("results/brain_endocast_comparison.pdf", height = 4, width = 4)
par(mfrow = c(1, 1), cex = 0.7, cex.axis = 0.65, bty = "l")
plot(1, type="l", xlim=c(0.35,1), ylim=c(0,14), xlab = "% of endocast filled by the brain",
     ylab="", yaxt="n",axes = T, cex.lab = 0.7, cex.axis = 0.65)
axis(2, at = 1, labels = "WOB", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 3, labels = "HOB", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 5, labels = "WCH", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 7, labels = "HCH", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 9, labels = "WOR", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 11, labels = "HOR", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)
axis(2, at = 13, labels = "BE", las=2, tck = -0.03, lty = 1, cex.axis = 0.65)

# plot lines and points for WOB 
points(WOB$be[which(WOB$juv == "y")], rep(1, length(WOB$be[which(WOB$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.1, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(WOB$be[which(WOB$juv == "y")], rep(1, length(WOB$be[which(WOB$juv == "y")])), 
       cex = 0.5, lab = endo_meas$id[which(WOB$juv == "y")])
points(WOB$be[which(WOB$juv == "n")], rep(1, length(WOB$be[which(WOB$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(WOB$be[which(WOB$juv == "n")], rep(1, length(WOB$be[which(WOB$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(WOB$juv == "n")])

# plot lines and points for HOB 
points(HOB$be[which(HOB$juv == "y")], rep(3, length(HOB$be[which(HOB$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(HOB$be[which(HOB$juv == "y")], rep(3, length(HOB$be[which(HOB$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(HOB$juv == "y")])
points(HOB$be[which(HOB$juv == "n")], rep(3, length(HOB$be[which(HOB$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(HOB$be[which(HOB$juv == "n")], rep(3, length(HOB$be[which(HOB$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(HOB$juv == "n")])

# plot lines and points for WCH 
points(WCH$be[which(WCH$juv == "y")], rep(5, length(WCH$be[which(WCH$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(WCH$be[which(WCH$juv == "y")], rep(5, length(WCH$be[which(WCH$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(WCH$juv == "y")])
points(WCH$be[which(WCH$juv == "n")], rep(5, length(WCH$be[which(WCH$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(WCH$be[which(WCH$juv == "n")], rep(5, length(WCH$be[which(WCH$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(WCH$juv == "n")])

# plot lines and points for HCH 
points(HCH$be[which(HCH$juv == "y")], rep(7, length(HCH$be[which(HCH$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(HCH$be[which(HCH$juv == "y")], rep(7, length(HCH$be[which(HCH$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(HCH$juv == "y")])
points(HCH$be[which(HCH$juv == "n")], rep(7, length(HCH$be[which(HCH$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(HCH$be[which(HCH$juv == "n")], rep(7, length(HCH$be[which(HCH$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(HCH$juv == "n")])

# plot lines and points for WOR 
points(WOR$be[which(WOR$juv == "y")], rep(9, length(WOR$be[which(WOR$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(WOR$be[which(WOR$juv == "y")], rep(9, length(WOR$be[which(WOR$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(WOR$juv == "y")])
points(WOR$be[which(WOR$juv == "n")], rep(9, length(WOR$be[which(WOR$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(WOR$be[which(WOR$juv == "n")], rep(9, length(WOR$be[which(WOR$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(WOR$juv == "n")])

# plot lines and points for HOR 
points(HOR$be[which(HOR$juv == "y")], rep(11, length(HOR$be[which(HOR$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(HOR$be[which(HOR$juv == "y")], rep(11, length(HOR$be[which(HOR$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(HOR$juv == "y")])
points(HOR$be[which(HOR$juv == "n")], rep(11, length(HOR$be[which(HOR$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(HOR$be[which(HOR$juv == "n")], rep(11, length(HOR$be[which(HOR$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(HOR$juv == "n")])

# plot lines and points for BE
points(BE$be[which(BE$juv == "y")], rep(13, length(BE$be[which(BE$juv == "y")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("brown2"))/255, alpha = 0.7), col = "#A9A9A9", pch=23)
text(BE$be[which(BE$juv == "y")], rep(13, length(BE$be[which(BE$juv == "y")])), 
     cex = 0.5, lab = endo_meas$id[which(BE$juv == "y")])
points(BE$be[which(BE$juv == "n")], rep(13, length(BE$be[which(BE$juv == "n")])), 
       cex = endo_meas$ML/5, lwd = 0.15, bg = rgb(t(col2rgb("darkcyan"))/255, alpha = 0.7), col = "#A9A9A9", pch=21)
text(BE$be[which(BE$juv == "n")], rep(13, length(BE$be[which(BE$juv == "n")])), 
     cex = 0.5, lab = endo_meas$id[which(BE$juv == "n")])


dev.off()


################################################################################
# 3. Semilandmark (outline) analyses
################################################################################
# 3a. PCA complete dataset
################################################################################

# import tps file
dorsal.raw.data = readland.tps("Data/landmark_curve.TPS", specID = "ID", readcurves = T)

# import classifiers
classifiers.GM = read.csv("Data/classifiers.csv", sep = " ", header = T)

# classifier as factor
classifiers.GM$group = as.factor(classifiers.GM$group)
is.factor(classifiers.GM$group)
classifiers.GM$juvenile = as.factor(classifiers.GM$juvenile)
is.factor(classifiers.GM$juvenile)

# import curve slide
semi.land = read.csv("Data/curveslide.csv")

# procrustes fit
gpa.data = gpagen(dorsal.raw.data, curves = semi.land, surfaces = NULL, PrinAxes = T,
                     max.iter = NULL, ProcD = F, Proj = T, print.progress = T)

# plot aligned data
plot(gpa.data)

# PCA
PCA.dorsal = gm.prcomp(gpa.data$coords)

# plot extremes PC1
med = mshape(gpa.data$coords)  # medium shape

PC1 = PCA.dorsal$x[,1]
preds1 = shape.predictor(gpa.data$coords, x= PC1, Intercept = F, pred1 = min(PC1),
                         pred2 = max(PC1)) # PC 1 extremes


# plot extremes PC2
PC2 = PCA.dorsal$x[,2]
preds2 = shape.predictor(gpa.data$coords, x= PC2, Intercept = F, pred1 = min(PC2), 
                         pred2 = max(PC2)) # PC 2 extremes

# Export extremes of PCs 1 and 2
pdf("results/PC1&2_extremes_complete.pdf", height = 11, width = 8)

par(mfrow = c(2,1))
plotRefToTarget(med, preds1$pred1, method = "points")
title(main = "PC1 minimum values")
plotRefToTarget(med, preds1$pred2, method = "points")
title(main = "PC1 maximum values")

plotRefToTarget(med, preds2$pred1, method = "points")
title(main = "PC2 minimum values")
plotRefToTarget(med, preds2$pred2, method = "points")
title(main = "PC2 maximum values")

dev.off()

# check PCA plot colored by juveniles vs. adults
par(mfrow = c(1, 1))
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "yes")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "yes")], cex = 1, pch = 23, 
       col = "brown2", bg = "brown2")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "no")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "no")], cex = 1, pch = 21, 
       col = "darkcyan", bg = "darkcyan")
text(PCA.dorsal$x[,2]~PCA.dorsal$x[,1], labels = rownames(gpa.data$data), cex = 0.3, pos = 4)
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", 
       col = c("darkcyan", "brown2"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))



# PCA plot colored by phylogenetic groups & age
# preparing data for plotting
# Juvenile chelids
juv.chelid1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelidae"))]
juv.chelid2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelidae"))]
# Juvenile chelonioids
juv.chelon1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelonioidea"))]
juv.chelon2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelonioidea"))]
# Juvenile chelydrids
juv.chelyd1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelydroidea"))]
juv.chelyd2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Chelydroidea"))]
# Juvenile pelomedusoids
juv.pelo1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Pelomedusoides"))]
juv.pelo2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Pelomedusoides"))]
# Juvenile testudinoids
juv.testu1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Testudinoidea"))]
juv.testu2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Testudinoidea"))]
# Juvenile Trionychia
juv.trion1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Trionychia"))]
juv.trion2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "yes") & 
                                        (classifiers.GM$group == "Trionychia"))]


# prepare points for plot
adu.chelid1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                       (classifiers.GM$group == "Chelidae"))]
adu.chelid2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                       (classifiers.GM$group == "Chelidae"))]
adu.chelon1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Chelonioidea"))]
adu.chelon2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Chelonioidea"))]
adu.chelyd1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Chelydroidea"))]
adu.chelyd2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Chelydroidea"))]
adu.pelo1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Pelomedusoides"))]
adu.pelo2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Pelomedusoides"))]
adu.testu1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Testudinoidea"))]
adu.testu2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Testudinoidea"))]
adu.trion1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Trionychia"))]
adu.trion2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Trionychia"))]
adu.angol1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Angolachelonia"))]
adu.angol2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Angolachelonia"))]
others1 = PCA.dorsal$x[, 2][which(classifiers.GM$Testudines == "non_crown")]
others2 = PCA.dorsal$x[, 1][which(classifiers.GM$Testudines == "non_crown")]


# Plot all data before color
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)

# plot juveniles colored points 
points(juv.chelid1~juv.chelid2, cex = 2, pch = 23, col = "#A9A9A9", 
       bg = "#7AC5CD")
points(juv.chelon1~juv.chelon2, cex = 2, pch = 23, col = "#A9A9A9", 
       bg = "#1C86EE")
points(juv.chelyd1~juv.chelyd2, cex = 2, pch = 23, col = "#A9A9A9", 
       bg = "#00EE76")
points(juv.pelo1~juv.pelo2, cex = 2, pch = 23, col = "#A9A9A9", 
       bg = "#8B0000")
points(juv.testu1~juv.testu2, cex = 2, pch = 23, col = "#A9A9A9", bg = "#EE9A00")
points(juv.trion1~juv.trion2, cex = 2, pch = 23, col = "#A9A9A9", 
       bg = "#EE6A50")

# plot adults colored points
points(adu.chelid1~adu.chelid2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#7AC5CD")
points(adu.chelon1~adu.chelon2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#1C86EE")
points(adu.chelyd1~adu.chelyd2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#00EE76")
points(adu.pelo1~adu.pelo2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#8B0000")
points(adu.testu1~adu.testu2, cex = 2, pch = 21, col = "#A9A9A9", bg = "#EE9A00")
points(adu.trion1~adu.trion2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#EE6A50")
points(adu.angol1~adu.angol2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#B23AEE")
points(others1~others2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#EE2C2C")

# plot legend for each group
legend.groups = levels(classifiers.GM$group)
legend.groups = legend.groups[-c(5, 6, 7, 9, 12)]
legend.groups = c(legend.groups, "early-diverging")

legend(x = "topleft", legend = legend.groups, cex = 0.9, bty = "n", pch = 21, pt.cex = 2,
       pt.bg = c("#B23AEE", "#7AC5CD", "#1C86EE", "#00EE76", "#8B0000", "#EE9A00",
                "#EE6A50", "#EE2C2C"))



################################################################################
# 3b. PCA adults only dataset
################################################################################

# import tps file
dorsal.adu.data = readland.tps("Data/landmark_curve_adults.TPS", specID = "ID", readcurves = T)

# import classifiers
classifiers.GM.adu = read.csv("Data/classifiers_adults.csv", sep = " ", header = T)

# make classifier as factor
classifiers.GM.adu$group = as.factor(classifiers.GM.adu$group)
is.factor(classifiers.GM.adu$group)
classifiers.GM.adu$juvenile = as.factor(classifiers.GM.adu$juvenile)
is.factor(classifiers.GM.adu$juvenile)

# curve
semi.land = read.csv("Data/curveslide.csv")

# procrustes fit
gpa.data.adu = gpagen(dorsal.adu.data, curves = semi.land, surfaces = NULL, PrinAxes = T,
                  max.iter = NULL, ProcD = F, Proj = T, print.progress = T)

# plot aligned data
plot(gpa.data.adu)

# PCA
PCA.adu = gm.prcomp(gpa.data.adu$coords)

# plot extremes PC1
med = mshape(gpa.data.adu$coords)  # medium shape

PC1.adu = PCA.adu$x[,1]
preds1.adu = shape.predictor(gpa.data.adu$coords, x= PC1.adu, Intercept = F, 
                             pred1 = min(PC1.adu), pred2 = max(PC1.adu)) # PC 1 extremes


# plot extremes PC2
PC2.adu = PCA.adu$x[,2]
preds2.adu = shape.predictor(gpa.data.adu$coords, x= PC2.adu, Intercept = F, 
                             pred1 = min(PC2.adu), pred2 = max(PC2.adu)) # PC 2 extremes

# Export extremes of PCs 1 and 2
pdf("results/PC1&2_extremes_adults.pdf", height = 11, width = 8)

par(mfrow = c(2,1))
plotRefToTarget(med, preds1.adu$pred1, method = "points")
title(main = "PC1 minimum values")
plotRefToTarget(med, preds1.adu$pred2, method = "points")
title(main = "PC1 maximum values")

plotRefToTarget(med, preds2.adu$pred1, method = "points")
title(main = "PC2 minimum values")
plotRefToTarget(med, preds2.adu$pred2, method = "points")
title(main = "PC2 maximum values")

dev.off()

# PCA plot: different colors and points by groups & age

# prepare points for plot
adu2.chelid1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Chelidae")]

adu2.chelid2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Chelidae")]

adu2.chelon1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Chelonioidea")]

adu2.chelon2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Chelonioidea")]

adu2.chelyd1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Chelydroidea")]

adu2.chelyd2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Chelydroidea")]

adu2.pelo1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Pelomedusoides")]

adu2.pelo2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Pelomedusoides")]

adu2.testu1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Testudinoidea")]

adu2.testu2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Testudinoidea")]

adu2.trion1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Trionychia")]

adu2.trion2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Trionychia")]

adu2.angol1 = PCA.adu$x[, 2][which(classifiers.GM.adu$group == "Angolachelonia")]

adu2.angol2 = PCA.adu$x[, 1][which(classifiers.GM.adu$group == "Angolachelonia")]

others12 = PCA.adu$x[, 2][which(classifiers.GM.adu$Testudines == "non_crown")]
others22 = PCA.adu$x[, 1][which(classifiers.GM.adu$Testudines == "non_crown")]

# plot on PDF
pdf("results/PCA-plots.pdf", height = 4, width = 7)
par(mfrow = c(1, 2), cex = 0.7, cex.axis = 0.8, bty = "l")

# PCA plot based on complete datase colored by juveniles vs. adults
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0, bty = "n")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "yes")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "yes")], cex = 2, pch = 23, 
       col = "#A9A9A9", bg = "brown2")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "no")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "no")], cex = 2, pch = 21, 
       col = "#A9A9A9", bg = "darkcyan")
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", pt.cex = 2, 
       col = c("black", "black"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))


# PCA plot based on adults only dataset colored by clade
plot(PCA.adu, axis1 = 1, axis2 = 2, cex = 0, bty = "n")

# plot adults colored points
points(adu2.chelid1~adu2.chelid2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#7AC5CD")
points(adu2.chelon1~adu2.chelon2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#1C86EE")
points(adu2.chelyd1~adu2.chelyd2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#00EE76")
points(adu2.pelo1~adu2.pelo2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#8B0000")
points(adu2.testu1~adu2.testu2, cex = 2, pch = 21, col = "#A9A9A9", bg = "#EE9A00")
points(adu2.trion1~adu2.trion2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#EE6A50")
points(adu2.angol1~adu2.angol2, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#B23AEE")
points(others12~others22, cex = 2, pch = 21, col = "#A9A9A9", 
       bg = "#EE2C2C")

# plot legend for each group
legend.groups = levels(classifiers.GM.adu$group)
legend.groups = legend.groups[-c(5, 6, 7, 9, 12)]
legend.groups = c(legend.groups, "early-diverging")

legend(x = "topleft", legend = legend.groups, cex = 0.9, bty = "n", pch = 21, pt.cex = 2,
       pt.bg = c("#B23AEE", "#7AC5CD", "#1C86EE", "#00EE76", "#8B0000", "#EE9A00",
                 "#EE6A50", "#EE2C2C"))

dev.off()


# Save Data
save.image(file = "results/Ferreiraetal.RData")
