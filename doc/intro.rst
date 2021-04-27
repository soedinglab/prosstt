What is PROSSTT?
################

Introduction & Motivation
*************************

Single-cell RNA sequencing (scRNA-seq) is revolutionizing the study of molecular biology by making comprehensive transcriptomic profiles available for single cells. This is challenging the traditional definition and understanding of cell types, as scRNA-seq has demonstrated that transcriptional heterogeneity is present within cell populations previously thought homogeneous. Even more excitingly, scRNA-seq is providing unprecedented insights into cellular differentiation. Cellular lineage trees can be reconstructed from snapshots of single-cell transcriptomes of cells caught at various stages of their differentiation process. The study of gene expression change along the reconstructed trees can uncover the intricate regulatory processes that govern development and differentiation.

Dedicated tree inference algorithms are indispensable for the study of scRNA-seq differentiation data. While various algorithms have been developed and successfully used, there are two overlapping concerns going forward:

- while the experimental field is clearly trending towards bigger datasets with more complex lineage trees, most published data consists of tree topologies with one or at most two bifurcations. Hence most algorithms have never reconstructed complex lineage trees, and it is not clear how they would perform.

- even if such datasets were available, it is currently very difficult to quantify method performance, since there exist no datasets with known ground truths, i.e. data with known intrinsic developmental time and cell identities that can be compared to algorithm predictions.

PROSSTT addresses both needs by producing realistic scRNA-seq datasets of differentiating cells.


How does it work?
*****************

PROSSTT generates simulated scRNA-seq datasets in four steps:

1. **Generate tree:** the topology of the lineage tree and the length of each branch are either sampled or set by the user.

2. **Simulate average gene expression along tree:** gene expression levels are a linear mixture of a small number K (default K=10) of functional expression programs. The time evolution of each expression program is simulated by random walks with a momentum term in each tree segment. The contribution weight of each expression program is drawn randomly for each gene.

3. **Sample cells from tree:** Three different sampling strategies are available: (1) sampling every point on the tree, (2) sampling via a density function or (3) sampling diffusely around selected tree points.

4. **Simulate UMI counts:** We simulate unique molecular identifier (UMI) counts using a negative binomial distribution where the variance of gene expression depends on average gene expression.

The end result of a PROSSTT simulation is a counts matrix of *N* cells x *G* genes, as well as a vector that contains the branch assignments and pseudotime values for each simulated cell.

Target users
************

This tool should be ideal for bioinformaticians currently developing tree inference methods, or comparing multiple methods to each other. PROSSTT can create datasets of varying difficulty, complexity and size and help test the limits of a method. In time PROSSTT could be extended to produce simulations from the reconstructed lineage trees of real datasets or include alternative noise models. We hope that PROSSTT will, directly or indirectly, help give biological insights into how to model and interpret scRNA-seq data.

Additional information can be found in the paper and the Supplemental Material.