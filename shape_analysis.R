################################################################################
# Script written for shape analyses of turtle endocasts
# Published in Ferreira et al. in press Paleoneurology book chapter
################################################################################

library(geomorph)
library(ggplot2)
getwd()

# import tps file
dorsal.data = readland.tps("Data/dorsal_curve.TPS", specID = "ID", readcurves = T)

# import classifiers
classifiers.GM = read.csv("Data/classifiers_GM.csv", sep = ";", header = T)

taxa = dimnames(dorsal.data)
taxa = taxa[[3]]
ordphyl = match(classifiers.GM$ID, taxa)
which(is.na(ordphyl))

# curve
semi.land = read.csv("Data/curveslide.csv")

# procrustes fit
dorsal.aligned = gpagen(dorsal.data, curves = semi.land, surfaces = NULL, PrinAxes = T,
                     max.iter = NULL, ProcD = F, Proj = T, print.progress = T)

# PCA
PCA.dorsal = gm.prcomp(dorsal.aligned$coords)

# PCA colored by juveniles vs. adults
plot(PCA.dorsal, axis1 = 1, axis2 = 2)
text(PCA.dorsal$x[,2]~PCA.dorsal$x[,1], labels = taxa, cex = 0.9, pos = 4)
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "yes")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "yes")], cex = 2, pch = 23, 
       col = "brown2", bg = "brown2")
points(PCA.dorsal$x[, 2][which(classifiers.GM$juvenile == "no")]~PCA.dorsal$x[, 1]
       [which(classifiers.GM$juvenile == "no")], cex = 2, pch = 21, 
       col = "brown2", bg = "brown2")


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

# Plot
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
adu.sand1 = PCA.dorsal$x[, 2][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Sandownidae"))]
adu.sand2 = PCA.dorsal$x[, 1][which((classifiers.GM$juvenile == "no") & 
                                        (classifiers.GM$group == "Sandownidae"))]
others1 = PCA.dorsal$x[, 2][which(classifiers.GM$Testudines == "non_crown")]
others2 = PCA.dorsal$x[, 1][which(classifiers.GM$Testudines == "non_crown")]


# Plot all data
plot(PCA.dorsal, axis1 = 1, axis2 = 2)

# plot juveniles
points(juv.chelid1~juv.chelid2, cex = 2, pch = 23, col = "darkgreen", 
       bg = "darkgreen")
points(juv.chelon1~juv.chelon2, cex = 2, pch = 23, col = "darkslateblue", 
       bg = "darkslateblue")
points(juv.chelyd1~juv.chelyd2, cex = 2, pch = 23, col = "darkslategray4", 
       bg = "darkslategray4")
points(juv.pelo1~juv.pelo2, cex = 2, pch = 23, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(juv.testu1~juv.testu2, cex = 2, pch = 23, col = "brown3", bg = "brown3")
points(juv.trion1~juv.trion2, cex = 2, pch = 23, col = "darkorchid", 
       bg = "darkorchid")

# plot adults
points(adu.chelid1~adu.chelid2, cex = 2, pch = 21, col = "darkgreen", 
       bg = "darkgreen")
points(adu.chelon1~adu.chelon2, cex = 2, pch = 21, col = "darkslateblue", 
       bg = "darkslateblue")
points(adu.chelyd1~adu.chelyd2, cex = 2, pch = 21, col = "darkslategray4", 
       bg = "darkslategray4")
points(adu.pelo1~adu.pelo2, cex = 2, pch = 21, col = "darkolivegreen3", 
       bg = "darkolivegreen3")
points(adu.testu1~adu.testu2, cex = 2, pch = 21, col = "brown3", bg = "brown3")
points(adu.trion1~adu.trion2, cex = 2, pch = 21, col = "darkorchid", 
       bg = "darkorchid")
points(adu.sand1~adu.sand2, cex = 2, pch = 21, col = "gold3", 
       bg = "gold3")
points(others1~others2, cex = 2, pch = 21, col = "gray23", 
       bg = "gray23")





# Save Data
save.image()
