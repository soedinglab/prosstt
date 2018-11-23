import numpy as np
import matplotlib.pyplot as plt

def plot(real_name, real, sim, alpha=0.5):
    sim_means, sim_vars, sim_zeros_gene, sim_zeros_cell, sim_totals = sim
    real_cells, real_genes = real
    fig, ax = plt.subplots(ncols=3, nrows=3)
    fig.set_size_inches(12, 12)

    # mean ~ variance
    ax[0][0].set_title("mean ~ variance")
    ax[0][0].set_xlabel("log2(x+1) of avg. gene expr.")
    ax[0][0].set_ylabel("log2(x+1) of gene expr. variance")
    s1 = ax[0][0].scatter(np.log2(real_genes.loc['means']+1),
                        np.log2(real_genes.loc['var']+1),
                        label=real_name, alpha=alpha)
    s2 = ax[0][0].scatter(np.log2(sim_means+1),
                        np.log2(sim_vars+1),
                        label="prosstt", alpha=alpha)

    # mean gene expression
    ax[0][1].set_title("avg. gene expr.")
    ax[0][1].set_ylabel("log2(x+1) of mean gene expr.")
    data = [np.log2(real_genes.loc['means']+1),
            np.log2(sim_means+1)]
    b1 = ax[0][1].boxplot(data, labels=[real_name, "simulated"])

    # variance
    ax[0][2].set_title("gene expr. variance")
    ax[0][2].set_ylabel("log2(x+1) of gene variance")
    data = [np.log2(real_genes.loc['var']+1),
            np.log2(sim_vars+1)]
    b1 = ax[0][2].boxplot(data, labels=[real_name, "simulated"])

    # percentage zeros per gene
    ax[1][0].set_title("%zeros per gene")
    ax[1][0].set_ylabel("percentage zeros per gene")
    data = [(real_genes.loc['zeros'] / real_cells.shape[1]),
            (sim_zeros_gene / len(sim_zeros_cell))]
    b1 = ax[1][0].boxplot(data, labels=[real_name, "simulated"])

    # percentage zeros per cell
    ax[1][1].set_title("%zeros per cell")
    ax[1][1].set_ylabel("percentage zeros per cell")
    data = [(real_cells.loc['zeros'] / real_genes.shape[1]),
            (sim_zeros_cell / len(sim_zeros_gene))]
    b1 = ax[1][1].boxplot(data, labels=[real_name, "simulated"])

    # mean-zeros relationship
    ax[1][2].set_title("mean ~ %zeros (gene)")
    ax[1][2].set_xlabel("log2(x+1) avg. gene expression")
    ax[1][2].set_ylabel("percentage zeros per gene")
    s1 = ax[1][2].scatter(np.log2(real_genes.loc['means']+1),
                        real_genes.loc['zeros'] / real_cells.shape[1],
                        label=real_name, alpha=alpha)
    s2 = ax[1][2].scatter(np.log2(sim_means+1),
                        sim_zeros_gene / len(sim_zeros_cell),
                        label="simulated", alpha=alpha)

    # library size
    ax[2][0].set_title("library size")
    ax[2][0].set_ylabel("log2(x) total UMIs per cell\nby #genes")
    data = [np.log2(real_cells.loc['total'] / real_genes.shape[1]),
            np.log2(sim_totals /  len(sim_zeros_gene))]
    b1 = ax[2][0].boxplot(data, labels=[real_name, "simulated"])

    # size/zero relationship
    ax[2][1].set_title("lib. size ~ %zeros (cell)")
    ax[2][1].set_xlabel("log2(x+1) avg. gene expression")
    ax[2][1].set_ylabel("percentage zeros per gene")
    s1 = ax[2][1].scatter(np.log2(real_cells.loc['total'] + 1),
                          real_cells.loc['zeros'] / real_genes.shape[1],
                          label=real_name, alpha=alpha)
    s2 = ax[2][1].scatter(np.log2(sim_totals + 1),
                          sim_zeros_cell / len(sim_zeros_gene),
                          label="simulated", alpha=alpha)
