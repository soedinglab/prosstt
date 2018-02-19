"""
A script that will produce a simulation with N bifurcations. Every branch
will have the same pseudotime length (50), and hyperparameters for all sanity
checks and count sampling take default values.
"""

import warnings
import argparse

import numpy as np
from numpy import random
import pandas as pd
import scipy as sp

import matplotlib.pyplot as plt
import pylab

from prosstt import simulation as sim
from prosstt import tree
from prosstt import sim_utils as sut

with warnings.catch_warnings():
    warnings.filterwarnings(message='.*Conversion of the second.*',
                            action='ignore',
                            category=FutureWarning,
                            module='h5py')
    import anndata as ad
    from scanpy.api.tl import diffmap


def save_files(job_id, save_dir, X, labs, brns, scalings, uMs, H, gene_scale, alpha, beta):
    # make the data more presentable by adding gene and cell names
    cell_names = ["cell_" + str(i) for i in range(X.shape[0])]
    gene_names = ["gene_" + str(i) for i in range(X.shape[1])]

    expr_mat = pd.DataFrame(X, columns=gene_names, index=cell_names).astype(int)
    cell_params = pd.DataFrame({"pseudotime": labs,
                                "branches": brns,
                                "scalings": scalings},
                               index=cell_names,
                               columns=["pseudotime", "branches", "scalings"])
    gene_params = pd.DataFrame({"alpha": alpha,
                                "beta": beta,
                                "genescale": gene_scale},
                               index=gene_names,
                               columns=["alpha", "beta", "genescale"])

    expr_mat.to_csv(save_dir + "/" + job_id + "_simulation.txt", sep="\t")
    cell_params.to_csv(save_dir + "/" + job_id + "_cellparams.txt", sep="\t")
    gene_params.to_csv(save_dir + "/" + job_id + "_geneparams.txt", sep="\t")

    np.savetxt(fname=save_dir + "/" + job_id + "_h.txt", X=H)

    for branch in uMs.keys():
        np.savetxt(fname=save_dir +
                         "/" +
                         job_id +
                         "_ums" +
                         str(branch) +
                         ".txt",
                   X=uMs[branch])


def save_params(job_id, save_dir, lineage_tree, rseed):
    paramfile = save_dir + "/" + job_id + "_params.txt"
    with open(paramfile, 'w') as out:
        out.write("Genes: " + str(lineage_tree.G) + "\n")
        out.write("pseudotimes: " + str(lineage_tree.time) + "\n")
        out.write("topology: " + str(lineage_tree.topology) + "\n")
        out.write("#modules: " + str(lineage_tree.modules) + "\n")
        out.write("random seed: " + str(rseed))


def plot_diff_map(X, pseudotime, brns):
    data = ad.AnnData(X)
    diffmap(adata=data)
    diff_map = data.obsm["X_diffmap"]
    cols = np.array(list(pylab.cm.Set1.colors))

    branch_names, indices = np.unique(brns, return_inverse=True)

    fig, ax = plt.subplots(ncols=2)
    fig.set_size_inches(w=9, h=4)
    ax[0].scatter(diff_map[:, 0], diff_map[:, 1], c=cols[indices])
    ax[0].set_title("branches")
    ax[1].scatter(diff_map[:, 0], diff_map[:, 1], c=pseudotime, cmap="viridis")
    ax[1].set_title("pseudotime")
    plt.show()


def main(job_id, save_dir, num_brpoints, plot):
    rseed = np.random.randint(1000)
    random.seed(rseed)

    # sample the parameters randomly:
    G = random.randint(100, 1001)
    gene_scale = np.exp(sp.stats.norm.rvs(loc=0.7, scale=1, size=G))

    alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=G))
    beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=G)) + 1

    num_branches = 2 * num_brpoints + 1
    top = tree.Tree.gen_random_topology(num_brpoints)

    branches = np.unique(np.array(top).flatten())
    time = {b: 50 for b in branches}
    # pseudotime smaller than 20 makes little sense as that usually does
    # not give the random walk enough time to create useful variation
    # br_lengths = random.random_integers(20, 100, num_branches)

    br_compl = random.randint(10, 50+1)
    # more than 50 components is too much and needlessly increases
    # runtime

    t = tree.Tree(topology=top, time=time, num_branches=num_branches,
                  modules=br_compl, G=G)

    Ms = {}
    uMs, Ws, H = sim.simulate_lineage(t, a=0.05, intra_branch_tol=0)
    gene_scale = sut.simulate_base_gene_exp(t, uMs)
    Ms = {}
    for branch in t.branches:
        Ms[branch] = np.exp(uMs[branch]) * gene_scale
    t.add_genes(Ms)

    X, pseudotime, brns, scalings = sim.sample_density(t, t.get_max_time(),
                                                 alpha=alpha, beta=beta)

    if plot:
        plot_diff_map(X, pseudotime, brns)
    # job_id = "test"
    # save_dir = "/home/npapado/Desktop"
    save_params(job_id, save_dir, t, rseed)
    save_files(job_id, save_dir, X, pseudotime, brns, scalings, uMs, H,
               gene_scale, alpha, beta)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Generate a simulated scRNAseq \
                                     dataset with an n-fold bifurcation.')
    PARSER.add_argument("-j", "--job", dest="job", help="Job ID (prepended to \
                        all generated files)", metavar="FILE")
    PARSER.add_argument("-o", "--out", dest="outdir", help="Directory where  \
                        output files are saved", metavar="FILE")
    PARSER.add_argument("-n", "--num_brpoints", dest="n",
                        help="How many branching points the simulation contains",
                        metavar="FILE")
    PARSER.add_argument("-p", "--plot", dest="plot", action='store_true',
                        help="Plot a diffusion map.")

    args = PARSER.parse_args()

    main(args.job, args.outdir, int(args.n), args.plot)
