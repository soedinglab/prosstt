library(viridis)
library(RColorBrewer)
library(destiny)

# this script assumes that the working directory is .../prosstt/scripts/
# if this is not the case please use the following lines so that all relative paths work properly
setwd("~/Documents/repos/prosstt/scripts/")

rotation_matrix <- function(comp1, comp2, angle) {
  rot <- diag(3)
  theta <- angle / 180 * pi
  rot[comp1, comp1] <- cos(theta)
  rot[comp1, comp2] <- -sin(theta)
  rot[comp2, comp1] <- sin(theta)
  rot[comp2, comp2] <- cos(theta)
  return(rot)
}

##########################
# guo plot
##########################

guo <- read.table(file="../data/guo/GuoDescription.txt",
                  sep="\t", header = T, stringsAsFactors = F, check.names = F, row.names = 1)
guoft <- read.table("../data/guo/GuoFeatures.txt")

data = guo
stage = as.factor(guoft$V4)
dm <- DiffusionMap(data, sigma = "global")
d = 3
trcoords <- sweep(eigenvectors(dm)[,1:d], eigenvalues(dm)[1:d], MARGIN=2, `*`)
plot(trcoords[,1:2])

rot <- rotation_matrix(comp1 = 1, comp2 = 3, angle = -30)
# an angle of +30 or -30 might be needed depending on the 
coords = trcoords[,1:3]
tcoords = coords %*% rot
toplot <- cbind(dm@eigenvectors[,2], tcoords[,1], dm@eigenvectors[,3])
plot(toplot[,1], toplot[,2], pch=16, col="black",
     xlab="DC2", ylab="transformed DC1", cex=1, cex.lab=1.5, cex.axis=1.5)

# svg(filename="guo.svg", height=4.8, width=6)
par(mar=c(4.1, 5.1, 2.1, 5.1))
pal_labels = viridis(length(unique(stage)), option = "viridis")
colorcells = pal_labels[stage]

plot(toplot[,1], toplot[,2], pch=16, col="black",
     xlab="DC2", ylab="transformed DC1", cex=2, cex.lab=1.5, cex.axis=1.5)
points(toplot[,1], toplot[,2], pch=16, col=colorcells,
     xlab="DC2", ylab="transformed DC1", cex=1.5, cex.lab=1.5, cex.axis=1.5)
legend("right", legend=paste(levels(stage), "C", sep=""), pt.bg = pal_labels, pch=21, title="Stage",
       xpd=TRUE, inset=c(-0.25,0), bty="n", cex=1.5)
dev.off()


##########################
# paul plot
##########################
data <- read.table(file="~/Documents/repos/merlot/inst/example/Paul2015.txt",
                   sep="\t", header = T, stringsAsFactors = F, check.names = F, row.names = 1)
annotation <- read.csv("~/Documents/data/paul/PaulCellsMapping.csv", header = F)
ddata <- data

# # three coloring schemes
# # one: according to clusters
# all_colors <- c("cyan1", "cyan4", "darkcyan",  "blue", "darkblue", "blue4", "navy", "darkgoldenrod", "darkgoldenrod1", "darkgoldenrod1", "gold", "bisque", "palegreen", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", "darkolivegreen", "chartreuse4", "darkgreen")
# colorcells <- all_colors[annotation$V2]
# 
# # two: according to lineage/cell type
# labels <- c("Ery", "Mk", "DC", "Mono", "Neu", "Eos", "Baso")
# longnames <- c("erythrocyte", "megakaryocyte", "dendritic cells", "monocyte", "neutrophil", "eosinophil", "basophil")
# cell_cols <- c("firebrick3", "deepskyblue3", "darkorchid1", "forestgreen", "chartreuse", "darkolivegreen3", "darkolivegreen1")
# 
# colorcells <- rep("black", dim(ddata)[1])
# 
# ery <- (annotation$V2 %in% c(1,2,3,4,5,6,7))
# mk <- (annotation$V2 %in% c(8))
# dc <- (annotation$V2 %in% c(11))
# mono <- (annotation$V2 %in% c(10, 14, 15))
# neu <- (annotation$V2 %in% c(9, 16, 17))
# eos <- (annotation$V2 %in% c(18))
# baso <- (annotation$V2 %in% c(12 ,13))
# 
# colorcells[ery] <- cell_cols[1]
# colorcells[mk] <- cell_cols[2]
# colorcells[dc] <- cell_cols[3]
# colorcells[mono] <- cell_cols[4]
# colorcells[neu] <- cell_cols[5]
# colorcells[eos] <- cell_cols[6]
# colorcells[baso] <- cell_cols[7]

# three: according to progenitor type
prog_cols <- c("orangered3", "gold2", "green3")
labels <- c("MEP", "CMP", "GMP")
longnames <- c("megakaryocyte/erythrocyte progenitors", "common myeloid progenitors", "granulocyte/macrophage progenitors")
cmp <- (annotation$V2 %in% 1:6)
gmp <- (annotation$V2 %in% 7:11)
mep <- (annotation$V2 %in% 12:19)
colorcells <- rep("black", dim(ddata)[1])
colorcells[cmp] <- prog_cols[1]
colorcells[gmp] <- prog_cols[2]
colorcells[mep] <- prog_cols[3]



dif <- DiffusionMap(ddata)
trcoords <- sweep(eigenvectors(dif)[,1:20], eigenvalues(dif)[1:20], MARGIN=2, `*`)
plot(trcoords[,1:2])

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/paul.svg", height=4.8, width=6)
par(mar=c(4.1, 5.1, 2.1, 5.1))

plot(trcoords[,1], trcoords[,2], pch=16, col="black",
     xlab="DC2", ylab="transformed DC1", cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(trcoords[,1], trcoords[,2], pch=16, col=colorcells,
      xlab="DC2", ylab="transformed DC1", cex=1., cex.lab=1.5, cex.axis=1.5)
legend("right", legend=labels, pt.bg=prog_cols, pch=21, title="Group", cex=1.5,
       xpd=TRUE, inset=c(-0.25,0), bty="n")
dev.off()


##########################
# treutlein plot
##########################
raw = read.delim("~/Documents/repos/hhtree/ProcessedDatasets/TreutleinDescription.txt", stringsAsFactors=FALSE, header = T)
labels = read.table(file="~/Documents/repos/hhtree/data/BarbraCellTimes.txt", sep="", header=T, stringsAsFactors = F)
tp = as.factor(labels$assignment)

dm <- DiffusionMap(raw, density.norm = T, verbose = T, sigma = "local")

toplot = eigenvectors(dm)
svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/treutlein.svg", height=4.8, width=6)
par(mar=c(4.1, 5.1, 2.1, 5.1))
pal_labels = brewer.pal(length(unique(tp)), "Set2")
colorcells = pal_labels[tp]
plot(toplot[,1], toplot[,2], pch=16, col="black",
     xlab="DC2", ylab="transformed DC1", cex=2., cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
points(toplot[,1], toplot[,2], pch=16, col=colorcells,
       xlab="DC2", ylab="transformed DC1", cex=1.5, cex.lab=1.5, cex.axis=1.5, cex.names=1.5)
# plot.new()
abridged <- c("induced", "early intrm", "early iN", "early Myo", "late intrm", "MEF", "Myo", "Neuron")
legend("right", legend=abridged, pt.bg = pal_labels, pch=21, title="Stage",
       xpd=TRUE, inset=c(-0.25,0), bty="n", cex=1.5)
dev.off()



##########################
# sim single bif plot
##########################
JobFolder <- "~/Documents/repos/prosstt/data/single/"
JobName <- "single_"
job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file=paste(job, "simulation.txt", sep=""), sep="\t", header=T, row.names=1)
data <- as.matrix(raw)
library_size = apply(data, 1, sum)
data <- 1/library_size * data

labels <- read.table(paste(job, "cellparams.txt", sep=""))$branches + 1
# calculate diffusion map ----
dif <- DiffusionMap(data)
plot3d(eigenvectors(dif)[,1:3], col=labels, size=7)

d=4
trcoords <- sweep(eigenvectors(dif)[,1:d], eigenvalues(dif)[1:d], MARGIN=2, `*`)
plot(trcoords[,1:2])

t <- 50
theta = t / 180 * pi
rot = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2, ncol=2)
coords = trcoords[,c(1,3)]
tcoords = coords %*% rot
toplot <- cbind(dif@eigenvectors[,2], tcoords[,1], dif@eigenvectors[,3])
plot(toplot[,1], toplot[,2], pch=16, col=labels,
     xlab="DC2", ylab="transformed DC1", cex=1, cex.lab=1.5, cex.axis=1.5)


CellCoordinates <- toplot[,1:2]

pal_labels <- brewer.pal(length(unique(labels)), "Set2")
colorcells=pal_labels[labels]

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/single_bif.svg", height=4.8, width=6)
par(mar=c(4.1, 5.1, 2.1, 5.1))
plot(CellCoordinates[,1], CellCoordinates[,2], pch=16, col="black",
     xlab="DC1", ylab="DC2", cex=2, cex.lab=1.5, cex.axis=1.5)
points(CellCoordinates[,1], CellCoordinates[,2], pch=16, col=colorcells,
       xlab="DC1", ylab="DC2", cex=1.5, cex.lab=1.5, cex.axis=1.5)
legend("right", legend=levels(ds$branch), pt.bg=pal_labels, pch=21, title="Branch",
       xpd=TRUE, inset=c(-0.25,0), bty="n", cex=1.5)
dev.off()


##########################
# sim double bif plot
##########################
JobFolder <- "~/Documents/repos/prosstt/data/double/"
JobName <- "double_"
job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file=paste(job, "simulation.txt", sep=""), sep="\t", header=T, row.names=1)
data <- as.matrix(raw)
library_size = apply(data, 1, sum)
data <- 1/library_size * data

labels <- read.table(paste(job, "cellparams.txt", sep=""))$branches + 1

# calculate diffusion map ----
dif <- DiffusionMap(data)
plot3d(eigenvectors(dif)[,1:3], pch=16, col=labels, main=test, size=7)

ds <- data.frame(DC1=eigenvectors(dif)[,1], DC2=eigenvectors(dif)[,2], DC3=eigenvectors(dif)[,3],
                 DC4=eigenvectors(dif)[,4], branch=labels)

t1 <- 10
t2 <- 0
t3 <- 100
rotx = rotation_matrix(1, 2, t1)
roty = rotation_matrix(1, 3, t2)
rotz = rotation_matrix(2, 3, t3)
rot <- rotx %*% roty %*% rotz
coords = dif@eigenvectors[,1:3]
tcoords = coords %*% rot
toplot <- cbind(tcoords[,1], tcoords[,2])
plot(toplot, col=labels, pch=16, cex=2)

# calculate tree
CellCoordinates <- toplot

pal_labels <- brewer.pal(length(unique(labels)), "Set2")
colorcells=pal_labels[labels]

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/double_bif.svg", height=4.8, width=6)
par(mar=c(4.1, 5.1, 2.1, 5.1))
plot(CellCoordinates[,1], CellCoordinates[,2], pch=16, col="black",
     xlab="DC1", ylab="DC2", cex=2, cex.lab=1.5, cex.axis=1.5)
points(CellCoordinates[,1], CellCoordinates[,2], pch=16, col=colorcells,
       xlab="DC1", ylab="DC2", cex=1.5, cex.lab=1.5, cex.axis=1.5)
legend("right", legend=unique(ds$branch), pt.bg=pal_labels, pch=21, title="Branch",
       xpd=TRUE, inset=c(-0.25,0), bty="n", cex=1.5)
dev.off()