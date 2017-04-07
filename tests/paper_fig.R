library(rgl)
library(SLICER)
library(monocle)

dimensions=4
cur <- "test8"
JobFolder <- paste("~/Desktop/benchmark3/", cur, "/", sep="")
JobName <- paste(cur, "_", sep="")

job <- paste(JobFolder, JobName, sep="")
raw <- read.table(file=paste(job, "simulation.txt", sep=""), sep="\t", header=T, row.names=1, stringsAsFactors = T)
data <- as.matrix(raw)
time <- read.table(paste(job, "pseudotimes.txt", sep=""))
time <- time - min(time) + 1
labels <- read.table(paste(job, "branches.txt", sep="")) + 1

dm <- DiffusionMap(data)
plot3d(eigenvectors(dm)[,c(1,2,3)], pch=16, col=labels$branch, size=8)
# now set up the files that we need for the rotation
write.table(eigenvectors(dm)[,1:4], "~/Desktop/paper_hhtree/dm_coords", sep="\t", col.names=F, row.names=F)


dim(data)

genes = select_genes(data)
# for branch selection
k = select_k(data[,genes], kmin=5)
data_lle = lle(data[,genes], m=dimensions, k)$Y
data_graph = conn_knn_graph(data_lle, 10)
ends = find_extreme_cells(data_graph, data_lle)
start = 1
cells_ordered = cell_order(data_graph, start)
stimes <- cells_ordered
sbranches = tryCatch( {assign_branches(data_graph, start)},
                      error=function(cond) {
                        print("SLICER failed with message:")
                        print(cond)
                        return(NA)
                      } )

# for plotting
k = 25
data_lle = lle(data[,genes], m=5, k)$Y
plot3d(data_lle[,c(1,2,5)], pch=16, col=labels$branch, size = 8)
write.table(data_lle, "~/Desktop/paper_slicer/lle_coords", sep="\t", col.names=F, row.names=F)




# set up monocle ----
exprs <- t(data)
mdata <- newCellDataSet(exprs, expressionFamily = negbinomial.size(), lowerDetectionLimit = 1)

# run monocle ----
mdata <- estimateSizeFactors(mdata)
mdata <- tryCatch( {estimateDispersions(mdata)}, 
                   error = function(cond) {
                     print(cond)
                     return(NA)
                   })

mbranches <- NA
mtimes <- NA
if (!is.na(mdata)) {
  mdata <- detectGenes(mdata, min_expr=0.1)
  
  disp_table <- dispersionTable(mdata)
  ordering_genes <- subset(disp_table, mean_expression >= 0.5 & dispersion_empirical >= 2*dispersion_fit)$gene_id
  
  mdata <- setOrderingFilter(mdata, ordering_genes)
  mdata <- reduceDimension(mdata, max_components=dimensions)
  mdata <- orderCells(mdata, reverse=FALSE)
  
  mbranches <- pData(mdata)$State
  mtimes <- pData(mdata)$Pseudotime
  
  # help monocle out with determining pseudotime:
  # ideal: make it start pseudotime at cell 1
  # but monocle sometimes will put cell_1 in an inner branch.
  new_root <- pData(mdata)$State[1]
  tryCatch( {
    mdata2 <- orderCells(mdata, root_state=new_root)
    mbranches <- pData(mdata2)$State
    mtimes <- pData(mdata2)$Pseudotime}, 
    error = function(cond) {
      print("heuristic found inner branch D:")
      print(cond)
    })
}
plot3d(t(reducedDimS(mdata)), pch=16, col=labels$branch, size=8)
write.table(t(reducedDimS(mdata)), "~/Desktop/paper_mon/ddr_coords", sep="\t", col.names=F, row.names=F)

# process the saved files with automate_rotation.py and generate the publication-friendly 2D coordinates
# now plot them
diffcoords <- t(read.table("~/Desktop/paper_hhtree/rotated_coords"))
slicercoords <- t(read.table("~/Desktop/paper_slicer/rotated_coords"))
moncoords <- t(read.table("~/Desktop/paper_mon/rotated_coords"))

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/slicer_emb.svg", height=2, width=2.08)
par(mar = rep(0, 4))
plot(slicercoords, pch=21, bg=labels$branch)
box(which = "plot")
dev.off()

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/mon_emb.svg", height=2, width=2.08)
par(mar = rep(0, 4))
plot(moncoords, pch=21, bg=labels$branch)
box(which = "plot")
dev.off()

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/dest_emb.svg", height=2, width=2.08)
par(mar = rep(0, 4))
plot(diffcoords, pch=21, bg=labels$branch)
box(which = "plot")
dev.off()

svg(filename="~/Documents/presentations/2017-04_PROSSTT_paper/pics/legend.svg")
plot.new()
legend("right", legend = 1:7, fill = 1:7)
dev.off()


