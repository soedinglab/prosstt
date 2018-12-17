import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches


def plot(real_name, real, sim, alpha=0.5):
    sim_means, sim_vars, sim_zeros_gene, sim_zeros_cell, sim_totals = sim
    real_cells, real_genes = real
    fig, ax = plt.subplots(ncols=4, nrows=2)
    fig.set_size_inches(20, 10)

    # mean ~ variance
    ax[1][1].set_title("mean ~ variance")
    ax[1][1].set_xlabel("log2(x+1) of avg. gene expr.")
    ax[1][1].set_ylabel("log2(x+1) of gene expr. variance")
    s2 = ax[1][1].scatter(np.log2(sim_means+1),
                        np.log2(sim_vars+1),
                        label="prosstt", alpha=alpha)
    s1 = ax[1][1].scatter(np.log2(real_genes.loc['means']+1),
                        np.log2(real_genes.loc['var']+1),
                        label=real_name, alpha=alpha)

    # mean gene expression
    ax[1][0].set_title("avg. gene expr.")
#     ax[1][0].set_xlabel("log2(x+1) avg. expr. data")
    ax[1][0].set_ylabel("log2(x+1) avg. expr. simulation")
#     ax[1][0].scatter(np.log2(real_genes.loc['means'] + 1),
#                      np.log2(sim_means+1), alpha=alpha, c="black")
#     lims = [
#         np.min([ax[1][0].get_xlim(), ax[1][0].get_ylim()]),
#         np.max([ax[1][0].get_xlim(), ax[1][0].get_ylim()]),
#     ]
#     ax[1][0].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    data = [np.log2(real_genes.loc['means'] + 1), np.log2(sim_means+1)]
    b1 = ax[1][0].boxplot(data, labels=[real_name, "simulated"])

    # variance
    ax[0][2].set_title("gene expr. variance")
    ax[0][2].set_ylabel("log2(x+1) of gene variance")
    data = [np.log2(real_genes.loc['var']+1),
            np.log2(sim_vars+1)]
    b1 = ax[0][2].boxplot(data, labels=[real_name, "simulated"])

    # percentage zeros per gene
    ax[0][0].set_title("%zeros per gene")
    ax[0][0].set_ylabel("percentage zeros per gene")
    data = [(real_genes.loc['zeros'] / real_cells.shape[1]),
            (sim_zeros_gene / len(sim_zeros_cell))]
    b1 = ax[0][0].boxplot(data, labels=[real_name, "simulated"])

    # percentage zeros per cell
    ax[0][1].set_title("%zeros per cell")
    ax[0][1].set_ylabel("percentage zeros per cell")
    data = [(real_cells.loc['zeros'] / real_genes.shape[1]),
            (sim_zeros_cell / len(sim_zeros_gene))]
    b1 = ax[0][1].boxplot(data, labels=[real_name, "simulated"])

    # mean-zeros relationship
    ax[1][2].set_title("mean ~ %zeros (gene)")
    ax[1][2].set_xlabel("log2(x+1) avg. gene expression")
    ax[1][2].set_ylabel("percentage zeros per gene")
    s2 = ax[1][2].scatter(np.log2(sim_means+1),
                        sim_zeros_gene / len(sim_zeros_cell),
                        label="simulated", alpha=alpha)
    s1 = ax[1][2].scatter(np.log2(real_genes.loc['means']+1),
                        real_genes.loc['zeros'] / real_cells.shape[1],
                        label=real_name, alpha=alpha)

    # library size
    ax[0][3].set_title("library size")
    ax[0][3].set_ylabel("log2(x) total UMIs per cell\nby #genes")
    data = [np.log2(real_cells.loc['total'] / real_genes.shape[1]),
            np.log2(sim_totals /  len(sim_zeros_gene))]
    b1 = ax[0][3].boxplot(data, labels=[real_name, "simulated"])

    ax[1][3].scatter(1, 1, label="simulated")
    ax[1][3].scatter(1, 1, label=real_name)
    ax[1][3].set_xlim(2, 3)
    ax[1][3].axis('off')
    ax[1][3].legend(loc=6, prop={'size': 20}, markerscale=3)

#     # size/zero relationship
#     ax[2][1].set_title("lib. size ~ %zeros (cell)")
#     ax[2][1].set_xlabel("log2(x+1) avg. gene expression")
#     ax[2][1].set_ylabel("percentage zeros per gene")
#     s1 = ax[2][1].scatter(np.log2(real_cells.loc['total'] + 1),
#                           real_cells.loc['zeros'] / real_genes.shape[1],
#                           label=real_name, alpha=alpha)
#     s2 = ax[2][1].scatter(np.log2(sim_totals + 1),
#                           sim_zeros_cell / len(sim_zeros_gene),
#                           label="simulated", alpha=alpha)


def plot_overlap(dm, origin, indices1, merlot_branches):
    fig, ax = plt.subplots(ncols=3, nrows=2)
    fig.set_size_inches(w=12, h=6)
    
    # plot the data 
    ax[0][0].scatter(dm[origin=="original", 3], dm[origin=="original", 4], c=cm.Set1(merlot_branches-1))
    ax[0][0].set_title("original data")
    ax[0][1].scatter(dm[origin=="simulated", 3], dm[origin=="simulated", 4], c=cm.Set1(indices1))
    ax[0][1].set_title("simulated data")
    ax[0][2].axis('off')
    
    progenitor = mpatches.Patch(color=cm.Set1(0), label='progenitor')
    nonskeletal = mpatches.Patch(color=cm.Set1(1), label='non-skeletal')
    intermediate = mpatches.Patch(color=cm.Set1(2), label='intermediate')
    cartilage = mpatches.Patch(color=cm.Set1(3), label='cartilage')
    bone = mpatches.Patch(color=cm.Set1(4), label='bone')
    branch_legend = [progenitor, nonskeletal, intermediate, cartilage, bone]
    ax[0][2].legend(loc=6, prop={'size': 10}, markerscale=2, handles=branch_legend)
    
    ax[1][0].scatter(dm[origin=="original", 3], dm[origin=="original", 4], label="original", alpha=0.05)
    ax[1][0].scatter(dm[origin=="simulated", 3], dm[origin=="simulated", 4], label="simulated", alpha=0.02)
    ax[1][0].set_title("original/simulated together")
    
    ax[1][1].axis('off')
    
    original = mpatches.Patch(color=cm.Set1(1), label='original')
    simulated = mpatches.Patch(color=cm.Set1(4), label='simulated')
    mix_legend = [original, simulated]
    ax[1][1].legend(loc=6, prop={'size': 10}, markerscale=2, handles=mix_legend)
    
    ax[1][2].axis('off')