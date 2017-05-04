library(rgl)
suppressPackageStartupMessages(library(rpgraph))
suppressPackageStartupMessages(library(destiny))
hhtree <- "~/Documents/repos/hhtree"
various <- paste(hhtree, "/scripts/various.R", sep="")
# evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")
treefuncs <- paste(hhtree, "/TreeTopologyFunctions.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(treefuncs))


dimensions=4
cur <- "test8"
JobFolder <- paste("~/Documents/repos/prosstt/tests/", cur, "/", sep="")
JobName <- paste(cur, "_", sep="")

job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file=paste(job, "simulation.txt", sep=""), sep="\t", header=T, row.names=1, stringsAsFactors = T)
data <- as.matrix(raw)
time <- read.table(paste(job, "pseudotimes.txt", sep=""))
time <- time - min(time) + 1
labels <- read.table(paste(job, "branches.txt", sep="")) + 1

dm <- DiffusionMap(data)
plot3d(eigenvectors(dm)[,c(1,2,3)], pch=16, col=labels$branch, size=8)


CellCoordinates <- eigenvectors(dm)[, 1:dimensions]
JobName <- paste(JobName, "CellCoord", sep="")

DataFile <- paste(JobFolder, JobName, "_TreeTopology.dat", sep="")

Topology = read_topology(DataFile)
N_yk = 100
mu_0 = 0.00625
lambda_0 = 2.03e-09

# scaling according to N_yk
mu=(N_yk-1)*mu_0
lambda=((N_yk-2)**3)*lambda_0

ElasticTree3DTopology= CalculateElasticTree(CellCoordinates = CellCoordinates, 
                                            Topology = Topology, N_yk = N_yk,
                                            input = "topology", plot=F, EP=lambda, RP=mu)

# plot_elastic_tree(CellCoordinates, ElasticTree3DTopology, 3, colorcells=labels$branch)
EmbeddedTree = GenesSpaceEmbedding(CellCoordinates, data, ElasticTree3DTopology)

dif_genes = differentially_expressed_genes_elasticpath(EmbeddedTree, 1, 2,  test_name="dcov", Replicates= 30000, mode="cells")
