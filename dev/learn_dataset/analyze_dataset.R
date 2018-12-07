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

mean_exp <- list()
for (i in 1:5) {
  b <- emb$Branches[[i]]
  emb$Nodes[b, ][emb$Nodes[b, ] <= 0] <- 1e-7
  mean_exp[[i]] <- emb$Nodes[b, ]
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
means <- array(0, c(100, dim(emb$CellCoords)[2]))
vars <- array(0, c(100, dim(emb$CellCoords)[2]))

for(n in 1:100) {
  node_cells <- emb$Cells2TreeNodes[,1][emb$Cells2TreeNodes[,2] == n]
  means[n,] <- apply(emb$CellCoords[node_cells, ], 2, mean)
  vars[n,] <- apply(emb$CellCoords[node_cells, ], 2, var)
}

genes <- dim(emb$CellCoords)[2]
alphas <- rep(0, genes)
betas <- rep(0, genes)

for (g in 1:genes) {
  m <- means[, g]
  s <- vars[, g]
  A <- cbind(rep(0, length(m)), m, m^2)
  b <- s #rep(1, length(m))
  G <- matrix(nrow=2, ncol=3, byrow = TRUE, data = c(0, 0, 1, 0, 1, 0))
  h <- c(1e-7, 1+1e-7)
  constrained_model <- tryCatch(
    {
      lsei(A = A, B = b, G = G, H = h, type=2)
    },
    error=function(cond) {
      message(cond)
      message("\nrescaling...")
      m <- means[, g] / 1000
      s <- vars[, g] / 1000
      A <- cbind(rep(0, length(m)), m, m^2)
      b <- s #rep(1, length(m))
      lsei(A = A, B = b, G = G, H = h, type=2)
    },
    finally={}
  )
  
  alphas[g] <- constrained_model$X[3]
  betas[g] <- constrained_model$X[2]
}

g <- 121
plot(means[,g], vars[,g])
x <- seq(min(means[,g]), max(means[,g]), 0.1)
points(x, alphas[g]*x^2+betas[g]*x, type="l", col="blue")

write.csv(alphas, "~/Documents/repos/prosstt/dev/learn_dataset/axolotl_alphas.txt")
write.csv(betas, "~/Documents/repos/prosstt/dev/learn_dataset/axolotl_betas.txt")
