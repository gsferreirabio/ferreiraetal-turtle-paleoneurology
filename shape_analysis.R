################################################################################
# Script written for volume and shape analyses of turtle endocasts
# Published in Ferreira et al. in press Paleoneurology book chapter
################################################################################

################################################################################
# Linear measurements based analyses
###############################################################
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(paleotree)
library(rcompanion)

getwd()


# Load image
load(file = "results/Ferreiraetal.RData")

##############################################################
# Voume Analyses
##############################################################

# Load the MCCT tree of Farina et al. in press
MCCT = dget("Data/Evers_20M_MCCT.txt")
plot(MCCT, cex = 0.2)

# Load endocast volume raw data
vol.data = read.table("Data/volume_data_clades.txt", header = T, row.names = 1)

# Drop tips in tree not in data
vol.MCCT = dropPaleoTip(MCCT, MCCT$tip.label[!(MCCT$tip.label %in% 
                                                       rownames(vol.data))])
plot(ladderize(vol.MCCT, right = F))  ## check the tree
axisPhylo()

name.check(vol.MCCT, vol.data)  ## check if names are the same

# is there a correlation between endocast and box volume?
plot(vol.data[,c("endocast_vol", "box_vol")])

# extract columns
endocast.vol = vol.data[, "endocast_vol"]
box.vol = vol.data[, "box_vol"]
eb.vol = vol.data[, "eb_vol"]
clades = vol.data$large_clades

# check normality  & transform data
plotNormalHistogram(endocast.vol)
qqnorm(endocast.vol, ylab = "Sample Quantiles for Endocast volume values")
qqline(endocast.vol)
endocastLog = log(endocast.vol)
plotNormalHistogram(endocastLog)

plotNormalHistogram(box.vol)
qqnorm(box.vol, ylab = "Sample Quantiles for Box volume values")
qqline(box.vol)
boxLog = log(box.vol)
plotNormalHistogram(boxLog)

plotNormalHistogram(eb.vol)
qqnorm(eb.vol, ylab = "Sample Quantiles for Endocast/box volume values")
qqline(eb.vol)
eb.volLog = log(eb.vol)
plotNormalHistogram(eb.volLog)

# give names to the objects
names(endocastLog) = names(boxLog) = names(eb.volLog) = rownames(vol.data)
names(clades) = rownames(vol.data)

# make a PGLS model of endocast volume ~ box.vol
pglsModel = gls(endocastLog ~ boxLog, correlation = corBrownian(phy = vol.MCCT),
                data = vol.data, method = "ML")
summary(pglsModel)
write.table(as.character(summary(pglsModel)), "results/pglsEndoBox.txt")

# plot log-transformed data with pglsModel line
plot(endocastLog ~ boxLog, xlab = "Log-transformed Box volume", 
     ylab = "Log-transformed Endocast volume")
abline(coef(pglsModel))

# plot untransformed data for comparison
plot(endocast.vol ~ box.vol, xlab = "Box volume", 
     ylab = "Endocast volume")
abline(lm(endocast.vol ~ box.vol))

# make a PGLS model of endocast ~ box * clades
pglsCladeModel = gls(endocastLog ~ boxLog*clades, 
                     correlation = corBrownian(phy = vol.MCCT), data = vol.data,
                     method = "ML")
summary(pglsCladeModel)
write.table(as.character(summary(pglsCladeModel)), "results/pglsCladeModel.txt")

# make a PGLS model of endocast ~ clades
pglsClades = gls(endocastLog ~ clades, correlation = corBrownian(phy = vol.MCCT),
                 data = vol.data, method = "ML")
summary(pglsClades)
write.table(as.character(summary(pglsClades)), "results/pglsClades.txt")

# ancestral states estimate
# and also compute variances & 95% confidence intervals for each node
fit.EBvol = fastAnc(vol.MCCT, eb.volLog, vars = TRUE, CI = TRUE)
range(eb.volLog)

# plot log-transformed endocast/box volume on the tree
obj = contMap(vol.MCCT, eb.volLog, plot = F)
obj = setMap(obj, colors=c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677"))
plot(obj, main = "Trait", legend=0.5*max(nodeHeights(vol.MCCT)), fsize=c(0.7,0.9),
     leg.txt = "log(Endocast volume)")

# plot them to a pdf
pdf("results/ASE&traigram.pdf")

plot(obj, main = "Trait", legend=0.5*max(nodeHeights(vol.MCCT)), fsize=c(0.7,0.9),
     leg.txt = "log(Endocast volume)")

dev.off()

# brain vs. endocast volume comparison
BrainEndoVol = read.csv("Data/brain_endocast_vol.csv", header = T, sep = ",")
BrainEndoVol$endocast = log(BrainEndoVol$endocast)
BrainEndoVol$brain = log(BrainEndoVol$brain)

cor(BrainEndoVol$endocast, BrainEndoVol$brain)  ## Pearson's correlation

# fit linear regression
reg.BrainEndo = lm(BrainEndoVol$brain ~ BrainEndoVol$endocast)
res.BrEn = summary(reg.BrainEndo)

# export anova table
write.table(anova(reg.BrainEndo), "results/anova_br_endo_vol.txt", sep = ",")

# plot 
plot(BrainEndoVol$brain ~ BrainEndoVol$endocast, bty = "l", pch = 18,
     xlab = "log-transformed Endocast volume", ylab = "log-transformed Brain volume")
title("A", adj = 0)
legend("topleft", paste("y = ", round(res.BrEn$coefficients[1], digits = 3), 
                        " + ", round(res.BrEn$coefficients[2], digits = 3), "
R^2 = ", round(res.BrEn$r.squared, digits = 3)), bty = "n")
abline(reg.BrainEndo, col = "#AA3377", lwd = 1.5)


##############################################################
# Shape Analyses
##############################################################

# import data & prepare data
linear.data = read.csv("Data/linear_data.csv", header = T, sep = ";", row.names = 1)

lin.meas = linear.data[ , 5:15]
rownames(lin.meas) = rownames(linear.data)

##############################################################
# exploring data
plot(lin.meas$ML, lin.meas$WOB)
plot(lin.meas$ML, lin.meas$WCH)
plot(lin.meas$ML, lin.meas$WOR)
plot(lin.meas$ML, lin.meas$WIE)
plot(lin.meas$ML, lin.meas$VWMO)
plot(lin.meas$ML, lin.meas$HOB)
plot(lin.meas$ML, lin.meas$HCH)
plot(lin.meas$ML, lin.meas$HOR)
plot(lin.meas$ML, lin.meas$HPE)
plot(lin.meas$WOB, lin.meas$WCH)

# check normality  & transform data
plotNormalHistogram(lin.meas[1])

lin.measLog = as.data.frame(matrix(nrow = 33, ncol = 11))

for (i in 1:length(lin.meas)) {
        lin.measLog[[i]] = log(lin.meas[[i]])
}

rownames(lin.measLog) = rownames(lin.meas)
colnames(lin.measLog) = colnames(lin.meas)

# PC Analyses
lin.pca = prcomp(lin.measLog)
summary(lin.pca)

# variance plots
par(mfrow = c(1, 2))
screeplot(lin.pca, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")
abline(h = 1, col = "#AA3377", lty = 5)
legend("topright", legend = c("Eigenvalue = 1"), col = c("#AA3377"), lty = 5, cex = 0.8)

cumvar = cumsum(lin.pca$sdev^2 / sum(lin.pca$sdev^2))
plot(cumvar[0:10], xlab = "PC #", ylab = "Explained variance", 
     main = "Cumulative variance plot")
abline(v = 2, col = "#228833", lty = 5)
abline(h = 0.9327118, col = "#228833", lty = 5)
legend("bottomright", legend = c("Cut-off @ PC2"), col = c("#228833"), lty = 5, 
       cex = 0.8)

par(mfrow = c(1,1))
plot(lin.pca$x[, 1], lin.pca$x[, 2], xlab = "PC1 (75.2%)", ylab = "PC2 (18.0%)", 
     main = "PC1 / PC2 linear measurements" )
text(lin.pca$x[, 2] ~ lin.pca$x[, 1], labels = rownames(lin.measLog), 
     cex = 0.6, pos = 4)

# boxplot with brain vs. endocast measurement comparison
br_ec = read.table("Data/br_ec_measures.txt", header = T, sep = " ")
br_ec$Model = as.factor(br_ec$Model)

# transform data dividing all measures by ML
br_ec.trans = br_ec

for(i in 5:15){
        br_ec.trans[i] = br_ec[i]/br_ec[5]
}

# plot transformed data
pdf("results/boxplots_br_ec.pdf", height = 3, width = 4)

par(cex = 0.5, cex.axis = 0.8, bty = "n", lwd = 0.3)
boxplot(br_ec.trans[, 6:15], boxfill = NA, border = NA, bty = "n", lwd = 0.3)
boxplot(br_ec.trans[br_ec.trans$Model == "endocast", 6:15], xaxt = "n", yaxt = "n", add = T, boxwex = 0.35, 
        boxfill = "#4477AA", lwd = 0.6, at = 1:ncol(br_ec.trans[, 6:15]) - 0.2)
boxplot(br_ec.trans[br_ec.trans$Model == "brain", 6:15], xaxt = "n", add = T, boxwex = 0.35, 
        boxfill = "#EE6677", lwd = 0.6, at = 1:ncol(br_ec.trans[, 6:15]) + 0.2)
legend("topright", legend = c("endocast", "brain"), col = c("#4477AA", "#EE6677"),
       bty = "n", pch = 20, pt.cex = 2, cex = 0.8)

dev.off()

################################################################################
# Semilandmark (outline) based analyses
# Dorsal view data
################################################################################
library(geomorph)

# import tps file
dorsal.raw.data = readland.tps("Data/dorsal_curve.TPS", specID = "ID", readcurves = T)

# import classifiers
classifiers.GM = read.csv("Data/classifiers_GM.csv", sep = ";", header = T)

# make classifier as factor
classifiers.GM$group = as.factor(classifiers.GM$group)
is.factor(classifiers.GM$group)
classifiers.GM$juvenile = as.factor(classifiers.GM$juvenile)
is.factor(classifiers.GM$juvenile)

# curve
semi.land = read.csv("Data/curveslide.csv")

# procrustes fit
gpa.data = gpagen(dorsal.raw.data, curves = semi.land, surfaces = NULL, PrinAxes = T,
                     max.iter = NULL, ProcD = F, Proj = T, print.progress = T)

# plot aligned data
plot(gpa.data)

# visualize shape changes
PC1.neg = plot(gpa.data$coords[,, "Proganochelys_quenstedti_2"])

PC1.pos = plot(gpa.data$coords[, , "Podocnemis_erythrocephala"])

PC2.neg = plot(gpa.data$coords[,, "Yuraramirim_montealtensis"])
        
PC2.pos = plot(gpa.data$coords[,, "Bothremys_cooki"])

# PCA
PCA.dorsal = gm.prcomp(gpa.data$coords)

# PCA colored by juveniles vs. adults
par(mfrow = c(1, 1))
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "yes")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "yes")], cex = 1, pch = 23, 
       col = "brown2", bg = "brown2")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "no")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "no")], cex = 1, pch = 21, 
       col = "darkcyan", bg = "darkcyan")
text(PCA.dorsal$x[,2]~PCA.dorsal$x[,1], labels = taxa, cex = 0.3, pos = 4)
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", 
       col = c("darkcyan", "brown2"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))

################################################################################
# PCA different colors and points by groups & age
# objects for classifiers

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


#########################
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

############################
# Plot all data before color
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)

# plot juveniles colored points 
points(juv.chelid1~juv.chelid2, cex = 1, pch = 23, col = "darkgreen", 
       bg = "darkgreen")
points(juv.chelon1~juv.chelon2, cex = 1, pch = 23, col = "darkslateblue", 
       bg = "darkslateblue")
points(juv.chelyd1~juv.chelyd2, cex = 1, pch = 23, col = "darkslategray4", 
       bg = "darkslategray4")
points(juv.pelo1~juv.pelo2, cex = 1, pch = 23, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(juv.testu1~juv.testu2, cex = 1, pch = 23, col = "brown3", bg = "brown3")
points(juv.trion1~juv.trion2, cex = 1, pch = 23, col = "darkorchid", 
       bg = "darkorchid")

# plot adults colored points
points(adu.chelid1~adu.chelid2, cex = 1, pch = 21, col = "darkgreen", 
       bg = "darkgreen")
points(adu.chelon1~adu.chelon2, cex = 1, pch = 21, col = "darkslateblue", 
       bg = "darkslateblue")
points(adu.chelyd1~adu.chelyd2, cex = 1, pch = 21, col = "darkslategray4", 
       bg = "darkslategray4")
points(adu.pelo1~adu.pelo2, cex = 1, pch = 21, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(adu.testu1~adu.testu2, cex = 1, pch = 21, col = "brown3", bg = "brown3")
points(adu.trion1~adu.trion2, cex = 1, pch = 21, col = "darkorchid", 
       bg = "darkorchid")
points(adu.angol1~adu.angol2, cex = 1, pch = 21, col = "gold3", 
       bg = "gold3")
points(others1~others2, cex = 1, pch = 21, col = "gray23", 
       bg = "gray23")

# plot legend for each group
legend.groups = levels(classifiers.GM$group)
legend.groups = legend.groups[-c(5, 6, 7, 9, 12)]
legend.groups = c(legend.groups, "Others")

legend(x = "topleft", legend = legend.groups, cex = 0.9, bty = "n", pch = 20, pt.cex = 2,
       col = c("gold3", "darkgreen", "darkslateblue", "darkslategray4", 
               "darkolivegreen3", "brown3", "darkorchid", "gray23"))


# plot on PDF
pdf("results/PCA_outlines.pdf", height = 4, width = 8)
par(mfrow = c(1, 2), cex = 0.5, cex.axis = 0.8, bty = "l")

# Plot all data before color
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)

# plot juveniles colored points 
points(juv.chelid1~juv.chelid2, cex = 1, pch = 23, col = "darkgreen", 
       bg = "darkgreen")
points(juv.chelon1~juv.chelon2, cex = 1, pch = 23, col = "darkslateblue", 
       bg = "darkslateblue")
points(juv.chelyd1~juv.chelyd2, cex = 1, pch = 23, col = "darkslategray4", 
       bg = "darkslategray4")
points(juv.pelo1~juv.pelo2, cex = 1, pch = 23, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(juv.testu1~juv.testu2, cex = 1, pch = 23, col = "brown3", bg = "brown3")
points(juv.trion1~juv.trion2, cex = 1, pch = 23, col = "darkorchid", 
       bg = "darkorchid")

# plot adults colored points
points(adu.chelid1~adu.chelid2, cex = 1, pch = 21, col = "darkgreen", 
       bg = "darkgreen")
points(adu.chelon1~adu.chelon2, cex = 1, pch = 21, col = "darkslateblue", 
       bg = "darkslateblue")
points(adu.chelyd1~adu.chelyd2, cex = 1, pch = 21, col = "darkslategray4", 
       bg = "darkslategray4")
points(adu.pelo1~adu.pelo2, cex = 1, pch = 21, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(adu.testu1~adu.testu2, cex = 1, pch = 21, col = "brown3", bg = "brown3")
points(adu.trion1~adu.trion2, cex = 1, pch = 21, col = "darkorchid", 
       bg = "darkorchid")
points(adu.angol1~adu.angol2, cex = 1, pch = 21, col = "gold3", 
       bg = "gold3")
points(others1~others2, cex = 1, pch = 21, col = "gray23", 
       bg = "gray23")

legend(x = "topleft", legend = legend.groups, cex = 0.9, bty = "n", pch = 20, pt.cex = 2,
       col = c("gold3", "darkgreen", "darkslateblue", "darkslategray4", 
               "darkolivegreen3", "brown3", "darkorchid", "gray23"))


# PCA colored by juveniles vs. adults
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "yes")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "yes")], cex = 1, pch = 23, 
       col = "brown2", bg = "brown2")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "no")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "no")], cex = 1, pch = 21, 
       col = "darkcyan", bg = "darkcyan")
text(PCA.dorsal$x[,2]~PCA.dorsal$x[,1], labels = taxa, cex = 0.5, pos = 4)
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", 
       col = c("darkcyan", "brown2"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))


dev.off()

################################################################################
# allometric regression
par(mfrow = c(1, 1))
allom.reg = procD.lm(two.d.array(gpa.data$coords)~gpa.data$Csize, iter = 999)
summary(allom.reg)

plot(allom.reg)

# plot residuals
plot(allom.reg$residuals)
text(allom.reg$residuals, labels = classifiers.GM$ID, cex = 0.7, pos = 4)

# plotAllometry with RegScore
allom.plot = plotAllometry(allom.reg, size = gpa.data$Csize, method = "RegScore")
text(allom.plot$plot.args$x, allom.plot$RegScore, labels = classifiers.GM$ID, 
     cex = 0.7, pos = 4)
################################################################################
# shape~specimen age

# gpa + classifiers data frame
gpa.data.class = c(gpa.data, classifiers.GM)

# regression analysis and plot
age.reg = procD.lm(gpa.data.class$coords~gpa.data.class$juvenile, iter = 999)
summary(age.reg)

################################################################################
# ontogenetic trajectory analysis
################################################################################

# load raw coordinates data
raw.ontog = readland.tps("Data/dorsal_curve_ontog.TPS", specID = "ID", 
                         readcurves = T)

# load classifiers
classifiers.GM.ontog = read.csv("Data/classifiers_GM_ontog.csv", header = T, sep = ";")

# make classifier as factor
classifiers.GM.ontog$group = as.factor(classifiers.GM.ontog$group)
is.factor(classifiers.GM.ontog$group)


gpa.ontog = gpagen(raw.ontog, curves = semi.land, surfaces = NULL, Proj = T, ProcD = F, 
                               print.progress = T)


# PCA
PCA.ontog = gm.prcomp(gpa.ontog$coords)

# plot PCA results
plot(PCA.ontog, axis1 = 1, axis2 = 2, cex = 0)
text(PCA.ontog$x[, 2]~PCA.ontog$x[, 1], labels = classifiers.GM.ontog$ID, cex = 0.7, pos = 4)
points(PCA.ontog$x[, 2][which(classifiers.GM.ontog$juvenile == "no")]~
               PCA.ontog$x[, 1][which(classifiers.GM.ontog$juvenile == "no")], 
       cex = 1.5, pch = 21, col = "darkcyan", bg = "darkcyan")
points(PCA.ontog$x[, 2][which(classifiers.GM.ontog$juvenile == "yes")]~
               PCA.ontog$x[, 1][which(classifiers.GM.ontog$juvenile == "yes")], 
       cex = 1.5, pch = 23, col = "brown3", bg = "brown3")
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", 
       col = c("darkcyan", "brown2"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))

################################################################################
# ontogenetic regression
ontog.reg = procD.lm(two.d.array(gpa.ontog$coords)~gpa.ontog$Csize, iter = 999)
summary(ontog.reg)


# export anova table 
anovaOntogReg = anova(ontog.reg)
write.table(anovaOntogReg$table, "results/anova_ontog_regres.txt", sep = ",")

plot(ontog.reg)

# plot residuals
plot(ontog.reg$residuals)
text(ontog.reg$residuals, labels = classifiers.GM.ontog$ID, cex = 0.7, pos = 4)

# ontogenetic plotAllometry with RegScore
ontog.plot = plotAllometry(ontog.reg, size = gpa.ontog$Csize, method = "RegScore", cex = 0)
points(ontog.plot$plot.args$x[which(classifiers.GM.ontog$juvenile == "yes")], 
               ontog.plot$RegScore[which(classifiers.GM.ontog$juvenile == "yes")], 
       cex = 1.5, pch = 23, col = "brown3", bg = "brown3")
points(ontog.plot$plot.args$x[which(classifiers.GM.ontog$juvenile == "no")], 
       ontog.plot$RegScore[which(classifiers.GM.ontog$juvenile == "no")], 
       cex = 1.5, pch = 23, col = "darkcyan", bg = "darkcyan")
text(ontog.plot$plot.args$x, ontog.plot$RegScore, labels = classifiers.GM.ontog$ID, 
     cex = 0.7, pos = 4)
legend(x = "topleft", legend = c("adult", "juvenile"), bty = "n", 
       col = c("darkcyan", "brown2"), pt.bg = c("darkcyan", "brown2"), pch = c(21, 23))



################################################################################
# 3D PGA data
tridim.data = read.table("Data/turtle_endocasts_procrustes.txt", header = T)
ID_string = rownames(lin.meas)
ID_string = ID_string[-30]
rownames(tridim.data) = ID_string
class(tridim.data$coord1)
# take a look
class(tridim.data)
dim(tridim.data)

warningMessages = ""
for (i in 1:length(tridim.data)) {
        if(is.numeric(tridim.data[[i]])){
                next
        }
        else{
                warning("The following column is not numeric: ", i)
        }
}
warnings()
tridim.data$coord51 = as.numeric(tridim.data$coord51)

# PCA
PCA.3d = gm.prcomp(tridim.data)


################################################################################
# plot figure

pdf("results/figure_quantitative.pdf", height = 8, width = 8)
par(mfrow = c(2, 2), bty = "l")
# 1st plot: size evolution
plot(BrainEndoVol$brain ~ BrainEndoVol$endocast, pch = 18,
     xlab = "log-transformed Endocast volume", ylab = "log-transformed Brain volume")
title("A", adj = 0)
legend("topleft", paste("y = ", round(res.BrEn$coefficients[1], digits = 3), 
                        " + ", round(res.BrEn$coefficients[2], digits = 3), "
R^2 = ", round(res.BrEn$r.squared, digits = 3)), bty = "n")
abline(reg.BrainEndo, col = "#AA3377", lwd = 1.5)

# 2nd plot: brain vs. endocast differences
boxplot(br_ec.trans[, 6:15], bty = "l", boxfill = NA, border = NA, lwd = 0.3)
boxplot(br_ec.trans[br_ec.trans$Model == "endocast", 6:15], xaxt = "n", yaxt = "n", add = T, boxwex = 0.35, 
        boxfill = "#4477AA", lwd = 0.6, at = 1:ncol(br_ec.trans[, 6:15]) - 0.2)
boxplot(br_ec.trans[br_ec.trans$Model == "brain", 6:15], xaxt = "n", add = T, boxwex = 0.35, 
        boxfill = "#EE6677", lwd = 0.6, at = 1:ncol(br_ec.trans[, 6:15]) + 0.2)
title("B", adj = 0)
legend("topright", legend = c("endocast", "brain"), col = c("#4477AA", "#EE6677"),
       bty = "n", pch = 20, pt.cex = 2, cex = 0.8)

# 3rd plot: PCA juveniles vs. adults shape
plot(PCA.dorsal, axis1 = 1, axis2 = 2, cex = 0)  ## empty plot

points(juv.chelid1~juv.chelid2, cex = 1, pch = 21, col = "#AA3377", 
       bg = "#AA3377")
points(juv.chelon1~juv.chelon2, cex = 1, pch = 22, col = "#AA3377", 
       bg = "#AA3377")
points(juv.chelyd1~juv.chelyd2, cex = 1, pch = 22, col = "#AA3377", 
       bg = "#AA3377")
points(juv.pelo1~juv.pelo2, cex = 1, pch = 21, col = "#AA3377", 
       bg = "#AA3377")
points(juv.testu1~juv.testu2, cex = 1, pch = 22, col = "#AA3377", bg = "#AA3377")
points(juv.trion1~juv.trion2, cex = 1, pch = 23, col = "#AA3377", 
       bg = "#AA3377")

# plot adults colored points
points(adu.chelid1~adu.chelid2, cex = 1, pch = 21, col = "#228833", 
       bg = "#228833")
points(adu.chelon1~adu.chelon2, cex = 1, pch = 22, col = "#228833", 
       bg = "#228833")
points(adu.chelyd1~adu.chelyd2, cex = 1, pch = 22, col = "#228833", 
       bg = "#228833")
points(adu.pelo1~adu.pelo2, cex = 1, pch = 21, col = "#228833", 
       bg = "#228833")
points(adu.testu1~adu.testu2, cex = 1, pch = 22, col = "#228833", bg = "#228833")
points(adu.trion1~adu.trion2, cex = 1, pch = 23, col = "#228833", 
       bg = "#228833")
points(adu.angol1~adu.angol2, cex = 1, pch = 25, col = "#228833", 
       bg = "#228833")
points(others1~others2, cex = 1, pch = 25, col = "#228833", 
       bg = "#228833")

# plot legend for each group
title("C", adj = 0)
legendGroups = c("juvenile", "adult", "Pleurodira", "Durocryptodira", 
                  "Trionychia", "Others")

legend("topleft", legend = legendGroups, cex = 0.8, bty = "n", pch = c(21, 21, 1, 0, 5, 6), pt.cex = 1.5,
       col = c("#AA3377", "#228833", "black", "black", "black", "black", "black", "black"),
       pt.bg = c("#AA3377", "#228833", "black", "black", "black", "black", "black", "black"))

dev.off()

# Save Data
save.image(file = "results/Ferreiraetal.RData")
