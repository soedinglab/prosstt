set.seed(42)

library(merlot)
library(limSolve)
# Read the example Guo dataset that is distributed with the package
tree <- readRDS("~/Documents/data/axolotl/ElasticTree_top50ClustBl-topMatMarker_OnlyBlastema_NEnd4.RDS")
emb <- readRDS("~/Documents/data/axolotl/EmbeddedTree_NEnd4_top50ClustBl-topMatMarker_OnlyBlastema.RDS")
pt <- readRDS("~/Documents/data/axolotl/Pseudotime_top50ClustBl-topMatMarker_OnlyBlastema_NEnd4.RDS")

write.csv(emb$Cells2Branches, "~/Documents/repos/prosstt/dev/learn_dataset/axolotl_merlot.txt")

branch_names <- c("progenitor", "nonskeletal", "intermediate", "bone", "cartilage")
time <- unlist(lapply(emb$Branches, length))
time

# now get the change in average expression:
for (i in 1:5) {
  b <- emb$Branches[[i]]
  emb$Nodes[b, ][emb$Nodes[b, ] <= 0] <- 1e-7
  write.csv(x = emb$Nodes[b, ],
            file = paste("~/Documents/repos/prosstt/dev/learn_dataset/", branch_names[i], ".csv", sep=""))
}

# calculate tree density (how much % of cells each node has):
tree_density <- rep(0, 100)
tabled <- table(emb$Cells2TreeNodes[,2])
tree_density[as.integer(names(tabled))] <- tabled
# tree_density <- tree_density / sum(tree_density)

for (i in 1:5) {
  b <- emb$Branches[[i]]
  dens <- tree_density[b]
  write.csv(x = dens,
            file = paste("~/Documents/repos/prosstt/dev/learn_dataset/", branch_names[i], "_density.csv", sep=""))
}

# calculate variance parameters
progenitor_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] %in% emb$Branches[[1]][1]]
nonskeletal_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] %in% emb$Branches[[2]][1]]
intermediate_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] %in% emb$Branches[[3]][1:2]]
bone_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] %in% emb$Branches[[4]][1]]
cartilage_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] %in% emb$Branches[[5]][1]]

CellCoordinates <- tree$CellCoords[,1:2]
plot(CellCoordinates, pch=16, cex=1.5)
points(CellCoordinates[progenitor_cells, ], col="red", pch=16, cex=1.2)
points(CellCoordinates[nonskeletal_cells, ], col="cyan", pch=16, cex=1.2)
points(CellCoordinates[intermediate_cells, ], col="green", pch=16, cex=1.2)
points(CellCoordinates[bone_cells, ], col="orange", pch=16, cex=1.2)
points(CellCoordinates[cartilage_cells, ], col="yellow", pch=16, cex=1.2)

progenitor_means <- apply(emb$CellCoords[progenitor_cells, ], 2, mean)
nonskeletal_means <- apply(emb$CellCoords[nonskeletal_cells, ], 2, mean)
intermediate_means <- apply(emb$CellCoords[intermediate_cells, ], 2, mean)
bone_means <- apply(emb$CellCoords[bone_cells, ], 2, mean)
cartilage_means <- apply(emb$CellCoords[cartilage_cells, ], 2, mean)

means <- rbind(progenitor_means, nonskeletal_means, intermediate_means, bone_means, cartilage_means)
rownames(means) <- NULL

progenitor_var <- apply(emb$CellCoords[progenitor_cells, ], 2, var)
nonskeletal_var <- apply(emb$CellCoords[nonskeletal_cells, ], 2, var)
intermediate_var <- apply(emb$CellCoords[intermediate_cells, ], 2, var)
bone_var <- apply(emb$CellCoords[bone_cells, ], 2, var)
cartilage_var <- apply(emb$CellCoords[cartilage_cells, ], 2, var)

vars <- rbind(progenitor_var, nonskeletal_var, intermediate_var, bone_var, cartilage_var)
rownames(vars) <- NULL

G <- dim(emb$CellCoords)[2]
alphas <- rep(0, G)
betas <- rep(0, G)

for (g in 1:G) {
  m <- means[, g]
  s <- vars[, g]
  A <- cbind(rep(0, length(m)), m, m^2)
  b <- s
  G <- matrix(nrow=2, ncol=3, byrow = TRUE, data = c(0, 0, 1, 0, 1, 0))
  h <- c(1e-7, 1+1e-7)
  constrained_model <- lsei(A = A, B = b, G = G, H = h, type=2)
  alphas[g] <- constrained_model$X[3]
  betas[g] <- constrained_model$X[2]
}

g <- 291
plot(means[,g], vars[,g])
x <- seq(min(means[,g]), max(means[,g]), 0.1)
points(x, alphas[g]*x^2+betas[g]*x, type="l", col="blue")

write.csv(alphas, "~/Documents/repos/prosstt/dev/learn_dataset/axolotl_alphas.txt")
write.csv(betas, "~/Documents/repos/prosstt/dev/learn_dataset/axolotl_betas.txt")
