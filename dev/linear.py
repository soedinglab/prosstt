import numpy as np
from numpy import random
import scipy as sp

# import scanpy for dimensionality reduction and plotting
import anndata as ad
from scanpy.api.tl import tsne
from scanpy.api.tl import umap
from scanpy.api import pp

# set viridis as the default color map
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
mpl.rcParams['image.cmap'] = 'viridis'

import tree
import tree_utils as tu
import simulation as sim
import sim_utils as sut
import count_model as cmod

def bla():
    print("actually running")
    rseed = 42
    np.random.seed(rseed)
    top = [["A", 'B'], ['B', 'C'], ['C', 'D'], ['D', 'E'], ['E', 'F']]
    branches = np.unique(np.array(top).flatten())
    time = {b: 30 for b in branches}
    G = 500
    t = tree.Tree(topology=top, G=G, time=time, num_branches=len(branches), branch_points=0, modules=40)
    uMs, Ws, Hs = sim.simulate_lineage(t, intra_branch_tol=-1, inter_branch_tol=0)
    gene_scale = sut.simulate_base_gene_exp(t, uMs)
    t.add_genes(uMs, gene_scale)
    alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=t.G))
    beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=t.G)) + 1
    X1, labs1, brns1, scalings1 = sim.sample_whole_tree(t, 5, alpha=alpha, beta=beta)
    print(labs1[0:5])    
    # normalize gene expression by library size
    X1 = (X1.transpose() / scalings1).transpose()
    data1 = ad.AnnData(np.log(X1+1))
    pp.neighbors(data1, use_rep='X', n_neighbors=700)
    umap(data1)
    dm1 = data1.obsm["X_umap"]
    job_id = "linear"
    save_dir = "/home/npapado/Desktop"
    tu.save_matrices(job_id, save_dir, X1, uMs, Hs)
    tu.save_params(job_id, save_dir, t, rseed)
    tu.save_cell_params(job_id, save_dir, labs1, brns1, scalings1)
    np.savetxt("/home/npapado/Desktop/linear_umap.csv", dm1)


if __name__ == "__main__":
    bla()
