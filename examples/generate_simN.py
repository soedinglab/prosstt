"""A script that will produce a simulation with N bifurcations."""

import numpy as np
from numpy import random
import pandas as pd
import scipy as sp
import argparse

from prosstt import simulation as sim
from prosstt import tree


def save_files(job_id, save_dir, X, labs, brns, uMs, Hs, gene_scale):
    # make the data more presentable by adding gene and cell names
    cell_names = ["cell_" + str(i) for i in range(X.shape[0])]
    gene_names = ["gene_" + str(i) for i in range(X.shape[1])]

    pdX = pd.DataFrame(X, columns=gene_names, index=cell_names).astype(int)
    pdLabs = pd.DataFrame(labs, index=cell_names, columns=["pseudotime"])
    pdLabs = pdLabs.astype(int)
    pdBrns = pd.DataFrame(brns, index=cell_names, columns=["branch"])
    pdBrns = pdBrns.astype(int)

    pdX.to_csv(save_dir + "/" + job_id + "_simulation.txt", sep="\t")
    pdLabs.to_csv(save_dir + "/" + job_id + "_pseudotimes.txt", sep="\t")
    pdBrns.to_csv(save_dir + "/" + job_id + "_branches.txt", sep="\t")

    num_branches = len(uMs)
    for i in range(num_branches):
        np.savetxt(fname=save_dir + "/" + job_id + "_ums" + str(i) + ".txt",
                   X=uMs[i])
        np.savetxt(fname=save_dir + "/" + job_id + "_hs" + str(i) + ".txt",
                   X=Hs[i])

    np.savetxt(fname=save_dir + "/" + job_id + "_gene_scale.txt", X=gene_scale)


def save_params(job_id, save_dir, G, alpha, beta, br_lengths, br_compl, rseed):
    paramfile = save_dir + "/" + job_id + "_params.txt"
    with open(paramfile, 'w') as out:
        out.write("Genes: " + str(G) + "\n")
        out.write("alpha: " + str(alpha) + "\n")
        out.write("beta: " + str(beta) + "\n")
        out.write("pseudotimes: " + str(br_lengths) + "\n")
        out.write("#modules: " + str(br_compl) + "\n")
        out.write("random seed: " + str(rseed))


def main(job_id, save_dir, num_brpoints):
    rseed = np.random.randint(1000)
    random.seed(rseed)

    # sample the parameters randomly:
    G = random.random_integers(100, 1000)

    alpha = np.exp(random.normal(loc=np.log(0.2), scale=0.3, size=G))

    beta = np.exp(random.normal(loc=np.log(2), scale=0.3, size=G))
    while not all(beta >= 1):
        beta[beta < 1] = np.exp(random.normal(loc=np.log(2),
                                scale=0.3, size=len(beta[beta < 1])))

    num_branches = 2 * num_brpoints + 1
    top = tree.Tree.gen_random_topology(num_brpoints)

    br_lengths = [50] * num_branches
    # pseudotime smaller than 20 makes little sense as that usually does
    # not give the random walk enough time to create useful variation
    # br_lengths = random.random_integers(20, 100, num_branches)

    br_compl = random.random_integers(10, 50)
    # more than 50 components is too much and needlessly increases
    # runtime

    t = tree.Tree(topology=top, time=br_lengths, branches=num_branches,
                  branch_points=num_brpoints, modules=br_compl, G=G)

    sample_time = np.arange(0, t.get_max_time())

    Ms = None
    while not sim.are_lengths_ok(Ms, abs_max=1000, rel_dif=0.05):
        uMs, Ws, Hs = sim.simulate_branching_data(t, tol=0.2)
        gene_scale = np.exp(sp.stats.norm.rvs(loc=0.8, scale=1, size=G))
        Ms = [np.zeros((t.time[i], G)) for i in range(num_branches)]
        for i in range(num_branches):
            Ms[i] = np.exp(uMs[i]) * gene_scale

    t.add_genes(Ms)

    X, labs, brns = sim.sample_data_balanced(1, G, t, sample_time, alpha, beta)

    # job_id = "test"
    # save_dir = "/home/npapado/Desktop"
    save_params(job_id, save_dir, G, alpha, beta, br_lengths, br_compl, rseed)
    save_files(job_id, save_dir, X, labs, brns, uMs, Hs, gene_scale)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a simulated scRNAseq \
                                     dataset with a triple bifurcation.')
    parser.add_argument("-j", "--job", dest="job", help="Job ID (prepended to \
                        all generated files)", metavar="FILE")
    parser.add_argument("-o", "--out", dest="outdir", help="Directory where  \
                        output files are saved", metavar="FILE")
    parser.add_argument("-n", "--num_brpoints", dest="n",
                        help="How many branching points the simulation contains",
                        metavar="FILE")

    args = parser.parse_args()

    main(args.job, args.outdir, int(args.n))
