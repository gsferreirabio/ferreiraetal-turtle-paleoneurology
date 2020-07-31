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

library(factoextra)
library(ggplot2)


library(FactoMineR)
library(ggfortify) # this lets ggplot2 know how to interpret PCA objects

getwd()

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
clades = vol.data$large_clades

# give names to the objects
names(endocast.vol) = names(box.vol) = names(clades) = rownames(vol.data)

# make a PGLS model of endocast volume ~ box.vol
pglsModel = gls(endocast.vol ~ box.vol, correlation = corBrownian(phy = vol.MCCT),
                data = vol.data, method = "ML")
summary(pglsModel)

plot(endocast.vol ~ box.vol)
abline(coef(pglsModel))

# make a PGLS model of endocast ~ box * clades
pglsCladeModel = gls(endocast.vol ~ box.vol*clades, 
                     correlation = corBrownian(phy = vol.MCCT), data = vol.data,
                     method = "ML")
summary(pglsCladeModel)

# make a PGLS model of endocast ~ clades
pglsClades = gls(endocast.vol ~ clades, correlation = corBrownian(phy = vol.MCCT),
                 data = vol.data, method = "ML")
summary(pglsClades)

# for ancestral state estimate we divide endocast by box volume
vol.index = endocast.vol/box.vol

# ancestral states estimate
# and also compute variances & 95% confidence intervals for each node
fit.endocast.vol = fastAnc(vol.MCCT, vol.index, vars = TRUE, CI = TRUE)
range(vol.index)

# plot endocast/box volume on the tree
obj = contMap(vol.MCCT, vol.index, plot=FALSE)
obj = setMap(obj, colors=c("blue", "cyan", "green", "yellow", "red"))
plot(obj, legend=0.7*max(nodeHeights(vol.MCCT)), fsize=c(0.7,0.9))

##############################################################
# Shape Analyses
##############################################################


# import data & prepare data
linear.data = read.csv("Data/linear_data.csv", header = T, sep = ";", row.names = 1)
View(linear.data)

lin.meas = linear.data[ , 4:14]
rownames(lin.meas) = rownames(linear.data)

##############################################################
#I DID THE FOLLOWING ANALYSES WITH cor_linear_measurements.csv
##############################################################
# exploring data
plot(lin.meas$cML, lin.meas$cWOB)
plot(lin.meas$cML, lin.meas$cWCH)
plot(lin.meas$cML, lin.meas$cWOR)
plot(lin.meas$cML, lin.meas$cWIE)
plot(lin.meas$cML, lin.meas$cVWMO)
plot(lin.meas$cML, lin.meas$cHOB)
plot(lin.meas$cML, lin.meas$cHCH)
plot(lin.meas$cML, lin.meas$cHOR)
plot(lin.meas$cML, lin.meas$cHPE)
plot(lin.meas$cWOB, lin.meas$cWCH)

ggplot(lin.meas, aes(x = cML, y = cWOB)) +
        geom_point(aes(), show.legend = T) +
        geom_smooth(method = "lm", se = T)

#
lin.meas.val = lin.meas[, -c(7:11)]
lin.meas = lin.meas[-11,]
lin.meas = lin.meas[-c(2,4,7,9,10,12,15,20,23),]
lin.meas = lin.meas[,-10]
View(lin.meas.val)
#



# PC Analyses
lin.pca = prcomp(lin.meas)
summary(lin.pca)

# extract eigenvalues
eig.val = get_eigenvalue(lin.pca)
eig.val

# plot variation
fviz_eig(lin.pca, addlabels = T, ylim = c(0,90))

# extract variable results
var = get_pca_var(lin.pca)
ind = get_pca_ind(lin.pca)

# plot PCA variables
fviz_pca_var(lin.pca, col.var = "blue")

# plot PCA biplot
fviz_pca_biplot(lin.pca, title = "Linear Measurements Turtle Endocasts")


# with log transformed data
log.lin.meas = log(lin.meas)

# PC Analyses
lin.log.pca = prcomp(log.lin.meas)

# extract eigenvalues
eig.val = get_eigenvalue(lin.log.pca)
eig.val

# plot variation
fviz_eig(lin.log.pca, addlabels = T, ylim = c(0,90))

# extract variable results
var = get_pca_var(lin.log.pca)
ind = get_pca_ind(lin.log.pca)

# plot PCA variables
fviz_pca_var(lin.log.pca, col.var = "blue")

# plot PCA biplot
fviz_pca_biplot(lin.log.pca, title = "Linear Measurements (Log) Turtle Endocasts")
plot(lin.log.pca$x)
text(lin.log.pca$x[, 2]~lin.log.pca$x[, 1], labels = row.names(lin.meas), cex = 0.8, pos = 4)


################################################################################
# Semilandmark (outline) based analyses
# Dorsal view data
################################################################################
library(geomorph)

# import tps file
dorsal.raw.data = readland.tps("Data/dorsal_curve.TPS", specID = "ID", readcurves = T)

# create a copy of the raw data
#raw = readland.tps("Data/dorsal_curve.TPS", specID = "ID", readcurves = T)

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
# PCA colored by groups & age
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

legend(x = "topleft", legend = legend.groups, cex = 0.9, bty = "n", 
       col = c("gold3", "darkgreen", "darkslateblue", "darkslategray4", 
               "darkolivegreen3", "brown3", "darkorchid", "gray23"), fill = 
               c("gold3", "darkgreen", "darkslateblue", "darkslategray4", 
                 "darkolivegreen3", "brown3", "darkorchid", "gray23"))

################################################################################
# allometric regression
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





# Save Data
save.image()
