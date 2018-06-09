#!/usr/bin/env python
# coding: utf-8
"""
This module contains all the functions that produce simulations. This includes
the simulation of expression programs, coefficients that map expr. programs to
genes, and different sampling strategies for (pseudotime, branch) pairs.
"""

import sys
import warnings

import numpy as np
from numpy import random
import pandas as pd
import scipy as sp

from prosstt import sim_utils as sut
from prosstt import count_model as cm


def simulate_branch(tree, relative_means, coefficients, branch, tol):
    topology = np.array(tree.topology)
    program = sim_expr_branch(tree.time[branch], tree.modules, cutoff=tol)
    relative_means[branch] = np.dot(program, coefficients)
    relative_means[branch] = sut.adjust_to_parent(relative_means, branch, topology)
    return program, relative_means


def simulate_expression_programs(tree, tol):
    """
    Simulate the relative expression of the lineage tree expression programs
    for each branch.

    Parameters
    ----------
    tol: float
        Correlation cut-off between expression programs

    Returns
    -------
    programs: dict
        Relative expression for all expression programs on every branch of the
        lineage tree
    """
    if tol > 1 or tol < 0:
        raise ValueError("value of 'tol' parameter should be between 0 and 1")
    programs = {}
    for branch in tree.branches:
        programs[branch] = sim_expr_branch(tree.time[branch], tree.modules,
                                           cutoff=tol)
    return programs


def sim_expr_branch(branch_length, expr_progr, cutoff=0.2, max_loops=100):
    """
    Return expr_progr diffusion processes of length T as a matrix W. The output of
    sim_expr_branch is complementary to _sim_coeff_beta.

    W encodes how a group of expr_progr coexpressed genes will behave through
    differentiation time. This matrix describes one branch of a
    differentiation tree (between two branch points or between a branch point
    and an endpoint). W describes the module in terms of relative expression
    (from 0 to a small positive float, so from a gene not being expressed to a
    gene being expressed at 2x, 3x of its "normal" level).

    After each new diffusion process is added the function checks whether the
    new diffusion correlates with any of the older ones. If the correlation is
    too high (above 0.5 per default), the last diffusion process will be
    replaced with a new one until one is found that does not correlate with any
    other columns of W or a suitable replacement hasn't been found after 100
    tries.

    Obviously this gets more probable the higher the number of components is -
    it might be advisable to change the number of maximum loops allowed or
    the cutoff in order to reduce runtime for a high number of components.

    Parameters
    ----------
    branch_length: int
        The length of the branch of the differentiation tree in pseudotime units
    expr_progr: int
        The number of components/modules of coexpression that describe the
        differentiation in this branch of the tree
    cutoff: float, optional
        Correlation above the cut-off will be considered too much. Should be
        between 0 and 1 but is not explicitly tested
    max_loops: int, optional
        The maximum number of times the method will try simulating a new
        diffusion process that doesn't correlate with all previous ones in W
        before resetting the matrix and starting over

    Returns
    -------
    W: ndarray
        Output array
    """
    programs = np.zeros((expr_progr, branch_length))
    k = 0
    loops = 0
    while k < expr_progr:
        programs[k] = diffusion(branch_length)

        correlates = sut.test_correlation(programs, k, cutoff)
        if correlates:
            # repeat and hope it works better this time
            loops += 1
            continue
        else:
            loops = 0
            k += 1

        if loops > max_loops:
            # we tried so hard
            # and came so far
            # but in the end
            # it doesn't even matter
            return sim_expr_branch(branch_length, expr_progr, cutoff=cutoff)

    return np.transpose(programs)


def diffusion(steps):
    """
    Diffusion process with momentum term. Returns a random walk with values
    usually between 0 and 1.

    Parameters
    ----------
    steps: int
        The length of the diffusion process.

    Returns
    -------
    walk: float array
        A diffusion process with a specified number of steps.
    """
    velocity = np.zeros(steps)
    walk = np.zeros(steps)

    # walk[0] = sp.stats.uniform.rvs()
    walk[0] = 0
    velocity[0] = sp.stats.norm.rvs(loc=0, scale=0.2)

    s_eps = 2 / steps
    eta = sp.stats.uniform.rvs()

    for t in range(0, steps - 1):
        walk[t + 1] = walk[t] + velocity[t]

        epsilon = sp.stats.norm.rvs(loc=0, scale=s_eps)
        # amortize the update
        velocity[t + 1] = 0.95 * velocity[t] + epsilon - eta * velocity[t]

    return walk


def simulate_coefficients(tree, fallback_a=0.04, **kwargs):
    """
    H encodes how G genes are expressed by defining their membership to K
    expression modules (coded in a matrix W). H could be told to encode
    metagenes, as it contains the information about which genes are coexpressed
    (genes that belong to/are influenced by the same modules). The influence of
    a module on a gene is measured by a number between 0 and 1, drawn from a
    (symmetric, if used with default values) beta distribution.

    The result of simulate_H is complementary to sim_expr_branch.

    Parameters
    ----------
    tree: Tree

    a: float, optional
        Shape parameter of Gamma distribution or first shape parameter of Beta
        distribution
    **kwargs: float
        Additional parameter (float b) if Beta distribution is to be used

    Returns
    -------
    A sparse matrix of the contribution of K expression programs to G genes.
    """
    if "a" not in kwargs.keys():
        warnings.warn(
            "No argument 'a' specified in kwargs: using gamma and a=0.04", UserWarning)
        return _sim_coeff_gamma(tree, fallback_a)
    # if a, b are present: beta distribution
    if "b" in kwargs.keys():
        groups = sut.create_groups(tree.modules, tree.G)
        return _sim_coeff_beta(tree, groups)
    else:
        return _sim_coeff_gamma(tree, a=kwargs['a'])


def _sim_coeff_beta(tree, groups, a=2, b=2):
    """
    Draw weights for the contribution of tree expression programs to gene
    expression from a Beta distribution.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    groups: list of ints
        A list of the two modules to which each gene belongs
    a: float, optional
        First shape parameter of the Beta distribution
    b: float, optional
        Second shape parameter of the Beta distribution

    Returns
    -------
    H: ndarray
        Output array
    """
    H = np.zeros((tree.modules, tree.G))
    for k in range(tree.modules):
        for gene in groups[k]:
            H[k][gene] += sp.stats.beta.rvs(a, b)
    return H


def _sim_coeff_gamma(tree, a=0.05):
    """
    Draw weights for the contribution of tree expression programs to gene
    expression from a Gamma distribution.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    a: float, optional
        Shape parameter of the Gamma distribution

    Returns
    -------
    H: ndarray
        Output array
    """
    K = tree.modules
    G = tree.G
    coefficients = np.reshape(sp.stats.gamma.rvs(a, size=K * G), (K, G))
    return coefficients


def simulate_lineage(tree, rel_exp_cutoff=8, intra_branch_tol=0.5,
                     inter_branch_tol=0, **kwargs):
    """
    Simulate gene expression for each point of the lineage tree (each
    possible pseudotime/branch combination). The simulation will try to make
    sure that a) gene expression programs within the same branch don't correlate
    too heavily and b) gene expression programs in parallel branches diverge
    enough.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    rel_exp_cutoff: float, optional
        The log threshold for the maximum average expression before scaling.
        Recommended values are below 9 (exp(9) is about 8100).
    intra_branch_tol: float, optional
        The threshold for correlation between expression programs in the same
        branch
    inter_branch_tol: float, optional
        The threshold for anticorrelation between relative gene expression in
        parallel branches
    **kwargs: various, optional
        Accepts parameters for coefficient simulation; float a if coefficients
        are generated by a Gamma distribution or floats a, b if the coefficients
        are generated by a Beta distribution

    Returns
    -------
    rel_means: Series
        Relative mean expression for all genes on every lineage tree branch
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree
    coefficients: ndarray
        Array that contains the contribution weight of each expr. program for
        each gene
    """
    if not len(tree.time) == tree.num_branches:
        print("the parameters are not enough for %i branches" %
              tree.num_branches)
        sys.exit(1)

    coefficients = simulate_coefficients(tree, **kwargs)
    bfs = sut.breadth_first_branches(tree)
    programs = {}
    rel_means = {}

    for branch in bfs:
        programs[branch], rel_means = simulate_branch(tree, rel_means,
                                                      coefficients, branch,
                                                      intra_branch_tol)
        above_cutoff = (np.max(rel_means[branch]) > rel_exp_cutoff)
        parallels = sut.find_parallel(tree, programs, branch)
        diverges = sut.diverging_parallel(parallels, rel_means, tree.G, tol=inter_branch_tol)
        while above_cutoff or not all(diverges):
            programs[branch], rel_means = simulate_branch(tree, rel_means,
                                                          coefficients, branch,
                                                          intra_branch_tol)
            above_cutoff = (np.max(rel_means[branch]) > rel_exp_cutoff)
            parallels = sut.find_parallel(tree, programs, branch)
            diverges = sut.diverging_parallel(parallels, rel_means, tree.G, tol=inter_branch_tol)

    return (pd.Series(rel_means),
            pd.Series(programs),
            coefficients)


def sample_whole_tree_restricted(tree, alpha=0.2, beta=3, gene_loc=0.8, gene_s=1):
    """
    Bare-bones simulation where the lineage tree is simulated using default
    parameters. Branches are assigned randomly if multiple are possible.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    alpha: float, optional
        Average alpha value
    beta: float, optional
        Average beta value
    gene_loc: float, optional
        Mean of the log-normal distribution of base gene expression values
    gene_s: float, optional
        Standard deviation of base gene expression value distribution
        (log-normal)

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    sample_pt: ndarray
        Pseudotime values of the sampled cells
    scalings: ndarray
        Library size scaling factor for each cell
    """
    sample_time = np.arange(0, tree.get_max_time())
    tree.default_gene_expression()
    alphas, betas = cm.generate_negbin_params(tree, mean_alpha=alpha, mean_beta=beta)

    return _sample_data_at_times(tree, sample_time, alpha=alphas, beta=betas)


def sample_pseudotime_series(tree, cells, series_points, point_std, alpha=0.3,
                             beta=2, scale=True, scale_v=0.7):
    """
    Simulate the expression matrix of a differentiation if the data came from
    a time series experimentree.

    Taking a sample from a culture of differentiating cells returns a mixture of
    cells at different stages of progress through differentiation (pseudotime).
    A time series experiment consists of sampling at multiple time points. This
    is simulated by drawing normally distributed pseudotime values around
    pseudotime sample points.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    cells: list or int
        If a list, then the number of cells to be sampled from each sample
        pointree. If an integer, then the total number of cells to be sampled
        (will be divided equally among all sample points)
    series_points: list
        A list of the pseudotime sample points
    point_std: list or float
        The standard deviation with which to sample around every sample pointree.
        Use a list for differing std at each time pointree.
    alpha: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    beta: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    scale: True, optional
        Apply cell-specific library size factor to average gene expression
    scale_v: float, optional
        Variance for the drawing of scaling factors (library size) for each cell

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    sample_pt: ndarray
        Pseudotime values of the sampled cells
    branches: ndarray
        The branch to which each simulated cell belongs
    scalings: ndarray
        Library size scaling factor for each cell
    """
    series_points, cells, point_std = sut.process_timeseries_input(
        series_points, cells, point_std)
    pseudotimes = []

    max_time = tree.get_max_time()
    for t, n, var in zip(series_points, cells, point_std):
        times_around_t = draw_times(t, n, max_time, var)
        pseudotimes.extend(times_around_t)
    return _sample_data_at_times(tree, pseudotimes, alpha=alpha, beta=beta,
                                 scale=scale, scale_v=scale_v)


def draw_times(timepoint, no_cells, max_time, var=4):
    """
    Draw cell pseudotimes around a certain sample time point under the
    assumption that in an asynchronously differentiating population cells are
    normally distributed around tree. The variance of the normal distribution
    controls the speed of differentiation (high spread: transient state/fast
    differentiation, low spread: bottleneck/slow differentiation).

    Parameters
    ----------
    timepoint: int
        The pseudotime point that represents the current mean differentiation
        stage of the population.
    no_cells: int
        How many cells to sample.
    max_time: int
        All time points that exceed the differentiation duration will be
        mapped to the end of the differentiation.
    var: float, optional
        Variance of the normal distribution we use to draw pseudotime points.
        In the experiment metaphor this parameter controls synchronicity.

    Returns
    -------
    sample_pt: int array
        Pseudotime points around <timepoint>.
    """
    sample_pt = sp.stats.norm.rvs(loc=timepoint, scale=var, size=no_cells)
    sample_pt = sample_pt.astype(int)
    sample_pt[sample_pt < 0] = 0
    sample_pt[sample_pt >= max_time] = max_time - 1
    return sample_pt


def sample_density(tree, no_cells, alpha=0.3, beta=2, scale=True, scale_v=0.7):
    """
    Use cell density along the lineage tree to sample pseudotime/branch pairs
    for the expression matrix.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    no_cells: int
        no_cellsumber of cells to sample
    alpha: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    beta: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    scale: True, optional
        Apply cell-specific library size factor to average gene expression
    scale_v: float, optional
        Variance for the drawing of scaling factors (library size) for each cell

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    sample_pt: ndarray
        Pseudotime values of the sampled cells
    branches: ndarray
        The branch to which each simulated cell belongs
    scalings: ndarray
        Library size scaling factor for each cell
    """
    bt = tree.branch_times()

    possible_pt = [np.arange(bt[b][0], bt[b][1] + 1) for b in tree.branches]
    possible_branches = [[b] * tree.time[b] for b in tree.branches]
    probabilities = [tree.density[b] for b in tree.branches]

    # make numpy arrays and flatten lists
    probabilities = np.concatenate(probabilities)
    possible_pt = np.concatenate(possible_pt)
    possible_branches = np.concatenate(possible_branches)

    # select according to density and take the selected elements
    sample = random.choice(np.arange(len(probabilities)),
                           size=no_cells, p=probabilities)
    sample_time = possible_pt[sample]
    sample_branches = possible_branches[sample]

    return _sample_data_at_times(tree, sample_time, alpha=alpha, beta=beta,
                                 branches=sample_branches, scale=scale,
                                 scale_v=scale_v)


def sample_whole_tree(tree, n_factor, alpha=0.3, beta=2, scale=True, scale_v=0.7):
    """
    Every possible pseudotime/branch pair on the lineage tree is sampled a
    number of times.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    n_factor: int
        How many times each pseudotime/branch combination can be present
    alpha: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    beta: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    scale: True, optional
        Apply cell-specific library size factor to average gene expression
    scale_v: float, optional
        Variance for the drawing of scaling factors (library size) for each cell

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    sample_pt: ndarray
        Pseudotime values of the sampled cells
    branches: ndarray
        The branch to which each simulated cell belongs
    scalings: ndarray
        Library size scaling factor for each cell
    """
    pseudotime, branches = cover_whole_tree(tree)

    branches = np.repeat(branches, n_factor)
    pseudotime = np.repeat(pseudotime, n_factor)

    return _sample_data_at_times(tree, pseudotime, alpha=alpha, beta=beta,
                                 branches=branches, scale=scale,
                                 scale_v=scale_v)


def cover_whole_tree(tree):
    """
    Get all the pseudotime/branch pairs that are possible in the lineage tree.

    Parameters
    ----------
    tree: Tree
        A lineage tree

    Returns
    -------
    pseudotime: ndarray
        Pseudotime values of all positions in the lineage tree
    branches: ndarray
        Branch assignments of all positions in the lineage tree
    """
    timezone = tree.populate_timezone()
    assignments = sut.assign_branches(tree.branch_times(), timezone)
    pseudotime = list()
    branches = list()

    for i, branch_timezone in enumerate(timezone):
        start = branch_timezone[0]
        end = branch_timezone[1] + 1
        length = end - start
        for branch in assignments[i]:  # for all possible branches in timezone a
            pseudotime.extend(np.arange(start, end))
            branches.extend([branch] * length)
    return pseudotime, branches


def _sample_data_at_times(tree, sample_pt, branches=None, alpha=0.3, beta=2,
                          scale=True, scale_v=0.7):
    """
    Sample cells from the lineage tree for given pseudotimes. If branch
    assignments are not specified, cells will be randomly assigned to one of the
    possible branches.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    sample_pt: ndarray
        Pseudotime values for the cells to be sampled
    branches: ndarray, optional
        Branch assignment of the cells to be sampled
    alpha: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    beta: float or ndarray, optional
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    scale: True, optional
        Apply cell-specific library size factor to average gene expression
    scale_v: float, optional
        Variance for the drawing of scaling factors (library size) for each cell

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    sample_pt: ndarray
        Pseudotime values of the sampled cells
    branches: ndarray
        The branch to which each simulated cell belongs
    scalings: ndarray
        Library size scaling factor for each cell
    """
    no_cells = len(sample_pt)
    if np.shape(alpha) == ():
        alpha = [alpha] * tree.G
    if np.shape(beta) == ():
        beta = [beta] * tree.G
    if branches is None:
        branches = sut.pick_branches(tree, sample_pt)
    scalings = sut.calc_scalings(no_cells, scale, scale_v)
    expr_matrix = draw_counts(tree, sample_pt, branches, scalings, alpha, beta)
    return expr_matrix, sample_pt, branches, scalings


def draw_counts(tree, pseudotime, branches, scalings, alpha, beta):
    """
    For all the cells in the lineage tree described by a given pseudotime and
    branch assignment, sample UMI count values for all genes. Each cell is an
    expression vector; the combination of all cell vectors builds the expression
    matrix.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    pseudotime: ndarray
        Pseudotime values for all cells to be sampled
    branches: ndarray
        Branch assignments for all cells to be sampled
    scalings: ndarray
        Library size scaling factor for all cells to be sampled
    alpha: float or ndarray
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray
    beta: float or ndarray
        Parameter for the count-drawing distribution. Float if it is the same
        for all genes, else an ndarray

    Returns
    -------
    expr_matrix: ndarray
        Expression matrix of the differentiation
    """
    no_cells = len(branches)

    cell_avg_exp = np.zeros((no_cells, tree.G))
    offsets = [tree.branch_times()[branch][0] for branch in branches]
    cell_times = pseudotime - offsets
    p_total = np.zeros(no_cells * tree.G)
    r_total = np.zeros(no_cells * tree.G)

    for n, time, branch in zip(np.arange(no_cells), cell_times, branches):
        cell_avg_exp[n] = tree.means[branch][time] * scalings[n]

    for cell in range(no_cells):
        p, r = cm.get_pr_umi(a=alpha, b=beta, m=cell_avg_exp[cell])
        p_total[cell * tree.G:(cell + 1) * tree.G] = p
        r_total[cell * tree.G:(cell + 1) * tree.G] = r

    nbinom = sp.stats.nbinom(n=r_total, p=(1 - p_total))
    expr_matrix = nbinom.rvs()
    # custm = cm.my_negbin()
    # expr_matrix = custm.rvs(p_total, r_total)
    return expr_matrix.reshape((no_cells, tree.G))
