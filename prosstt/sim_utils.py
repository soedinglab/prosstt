#!/usr/bin/env python
# coding: utf-8
"""
This module contains utility functions for the simulations, such as a printed
progress bar, functions to perform quality checks, functions to create and
assign groups, and functions to pick between alternatives when multiple options
are possible.
"""

import collections
from collections import defaultdict
from collections import deque
import operator
import numbers
import sys

import numpy as np
from numpy import random
import scipy as sp


def print_progress(iteration, total, prefix='', suffix='', decimals=1):
    """
    Call in a loop to create a terminal-friendly text progress bar. Contributed
    by Greenstick on stackoverflow.com/questions/3173320.

    Parameters
    ----------
        iteration: int
            Current iteration.
        total: int
            Total number of iterations.
        prefix: str, optional
            Prefix string before the progress bar.
        suffix: str, optional
            Suffix string after the progress bar.
        decimals: int, optional
            Positive number of decimals in percent complete.
    """
    bar_length = 80
    format_str = "{0:." + str(decimals) + "f}"
    percent = format_str.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    progress_bar = '█' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s |%s| %s%s %s' %
                     (prefix, progress_bar, percent, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def random_partition(k, iterable):
    """
    Random partition in almost equisized groups.

    Parameters
    ----------
    k: int
        How many partitions to create.
    iterable: array
        The iterable to be partitioned.

    Returns
    -------
    results: list of int lists.

    contributed by kennytm on stackoverflow.com/questions/3760752
    """
    results = [[] for i in range(k)]
    for value in iterable:
        x = random.randint(k)
        results[x].append(value)
    return results


def test_correlation(W, k, cutoff):
    """
    For a column of a matrix, test if previous columns correlate with it.

    Parameters
    ----------
    W: numpy array
        The matrix to test.
    k: int
        Compare columns from 0 to k-1 with column k.
    cutoff: float
        Correlation above the cut-off will be considered too much. Should be
        between 0 and 1 but is not explicitly tested.
    """
    for i in range(k - 1, 0):
        pearson_r = sp.stats.pearsonr(W[k], W[i])
        if pearson_r[0] > cutoff:
            return True
    return False


def create_groups(no_programs, no_genes):
    """
    Returns a list of the groups to which each gene belongs.

    Each gene g is assigned one of no_programs possible groups twice (random
    draw with replacement).

    Parameters
    ----------
    K: int
        Number of modules.
    G: int
        Number of genes.

    Returns
    -------
    groups: list of ints
        A list of the two modules to which each gene belongs.
    """
    genes = random.permutation(no_genes)
    # we want each gene to appear in two groups, in average.
    # If we draw twice it will happen that some genes will take the same
    # group twice, but it should not happen too often.
    groups1 = random_partition(no_programs, genes)
    # performing the permutation a second time is necessary, else most genes
    # will be in the same modules and we want to mix more
    genes = random.permutation(no_genes)
    groups2 = random_partition(no_programs, genes)
    groups = [[i for subz in z for i in subz] for z in zip(groups1, groups2)]
    return groups


def bifurc_adjust(child, parent):
    """
    Adjust two matrices so that the last line of one equals the first of the
    other.

    Parameters
    ----------
    child: matrix to be adjusted

    parent: matrix to adjust to
    """
    dif = child[0] - parent[-1]
    child = child - dif
    return child


def pearson_between_programs(genes, prog1, prog2):
    """
    Calculate the pearson correlation coefficient between two expression
    programs for all genes.

    Parameters
    ----------
    genes: int
        The number of genes in the lineage tree
    prog1: numpy.ndarray
        The first expression program
    prog2: numpy.ndarray
        The second expression program

    Returns
    -------
    pearson: numpy.ndarray
        The pearson correlation coefficient for all genes in the two programs
    """
    pearson = np.zeros(genes)
    common = min(prog1.shape[0], prog2.shape[0])
    for gene in range(genes):
        pearson[gene] = sp.stats.pearsonr(prog1[:common, gene], prog2[:common, gene])[0]
    return pearson


def flat_order(n):
    """
    Map from indices of flat array of size n(n-1)/2 to an upper triangular
    matrix of size nxn

    Parameters
    ----------
    n: int
        number of options to combine
    """
    size = int(n * (n - 1) / 2)
    res = np.zeros((size, 3), dtype=int)
    for i in range(n - 1):
        for j in range(i + 1, n):
            index = int((i * (2 * n - i - 3) / 2 + j - 1))
            res[index] = np.array([index, i, j])
    return res


def calc_relat_means(tree, programs, coefficients):
    """
    Calculate relative mean expression for a lineage tree given the expression
    programs and the coefficient matrix that contains the contribution of each
    expression program to each gene.

    Parameters
    ----------
    tree: Tree
        A lineage tree object.
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree.
    coefficients: numpy.ndarray
        Array that contains the contribution weight of each expr. program for
        each gene
    
    Returns
    -------
    """
    relative_means = {}
    for branch in tree.branches:
        relative_means[branch] = np.dot(programs[branch], coefficients)
    return relative_means


def diverging_parallel(branches, programs, genes, tol=0.5):
    """
    Calculate if the expression programs in all pairs of parallel branches are
    diverging enough to make the branches distinguishable from each other.

    Parameters
    ----------
    branches: list
        A list of pairs of parallel branches.
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree.
    genes: int
        The number of genes included in the lineage tree.
    tol: float, optional
        The percentage of genes that must have anticorrelated expression
        patterns over pseudotime in order for the branches to be considered
        diverging.

    Returns
    -------
    diverging: numpy.ndarray
        A list of the boolean values: whether each pair of parallel branches
        diverges or not.
    """
    branches = [b for b in branches if b is not None]
    if len(branches) == 1:
        return [True]
    indices = flat_order(len(branches))
    diverging = np.zeros(len(indices), dtype=bool)
    for index, i, j in indices:
        br1 = branches[i]
        br2 = branches[j]
        pearson = pearson_between_programs(genes, programs[br1], programs[br2])
        percent_anticorrelated = sum(pearson < 0) / (genes * 1.0)
        diverging[index] = percent_anticorrelated > tol
    return diverging


def assign_branches(branch_times, timezone):
    """
    Assigns a branch to every timezone::

                -- T[1]------
        -T[0]--|          -- T[3]------
                -- T[2]--|
                          -- T[4]-
        timezones:
        ---0---|----1----|-2-|-3--|-4--

    ========  ========
    timezone  branch
    ========  ========
    0         0
    1         1,2
    2         1,3,4
    3         3,4
    4         3
    ========  ========

    A time point in timezone i can belong to one of k possible branches

    Parameters
    ----------
    branch_times: list of int lists
        The pseudotime at which branches start and end.
    timezone: int array
        Array that contains the timezone information for each pseudotime point.

    Returns
    -------
    res: list of int lists
        A list of the possible branches for each timezone.
    """
    res = defaultdict(list)
    for i, zone in enumerate(timezone):
        for k in branch_times:
            branch_time = branch_times[k]
            if belongs_to(zone, branch_time):
                res[i].append(k)
    return res


def belongs_to(timezone, branch):
    """
    Checks whether a timezone start and end are contained within the
    pseudotime of a branch.

    Timezones are constructed such that they don't go over branch boundaries.
    This method is used to determine which branches are possible for a
    timezone.

    Parameters
    ----------
    timezone: int array
        The pseudotime at which the timezone starts and ends.
    branch: int array
        The pseudotime at which the branch starts and ends.

    Returns
    -------
    bool
        Whether the timezone is contained within the branch.
    """
    return (timezone[0] >= branch[0]) and (timezone[1] <= branch[1])


def pick_branches(tree, pseudotime):
    """
    Lorem ipsum dolor et ames.

    Parameters
    ----------
    a: int
        Lorem ipsum

    Returns
    -------
    b: int
        Lorem ipsum as well
    """
    timezone = tree.populate_timezone()
    assignments = assign_branches(tree.branch_times(), timezone)
    branches = np.zeros(len(pseudotime), dtype=str)
    for n, t in enumerate(pseudotime):
        branches[n] = pick_branch(tree, t, timezone, assignments)
    return branches


def pick_branch(tree, pseudotime, timezones, assignments):
    """
    Picks one of the possible branches for a cell at a given time point.

    Parameters
    ----------
    pseudotime: int
        A pseudotime point.
    timezones: int array
        The pseudotimes at which the timezones start and end.
    assignments: int array
        A list of the possible branches for each timezone.

    Returns
    -------
    branch: int
        The branch to which the cell belongs.
    """
    branch = -1
    for i, zone in enumerate(timezones):
        if pseudotime >= zone[0] and pseudotime <= zone[1]:
            branch = i
            break
    possibilities = assignments[branch]
    where_in_branch = pseudotime - timezones[branch][0]
    densities = np.zeros(len(possibilities))
    for i, b in enumerate(possibilities):
        densities[i] = tree.density[b][where_in_branch]
    probabilities = densities / densities.sum()
    try:
        return random.choice(possibilities, p=probabilities)
    except IndexError:
        print(pseudotime)
        print(timezones)
        print(assignments)


def max_relat_exp(tree, relative_means):
    """
    Finds maximum relative gene expression for each gene along the lineage tree.

    Parameters
    ----------
    tree: Tree
        A lineage tree object.
    relative_means: Series
        Relative mean expression for all genes on every lineage tree branch.

    Returns
    -------
    maxes: numpy.ndarray
        An array with the maximum relative expression of each gene along the
        lineage tree.
    """
    maxes = np.zeros((tree.G, len(tree.branches)))
    for i, branch in enumerate(tree.branches):
        maxes[:, i] = np.max(np.exp(relative_means[branch]), axis=0)
    return maxes


def simulate_base_gene_exp(tree, relative_means, abs_max=5000, gene_mean=0.8, gene_std=1):
    """
    Samples appropriate base expression values for each gene. The criterion
    applied is that the absolute average gene expression does not surpass a
    certain threshold.

    Parameters
    ----------
    tree: Tree
        A lineage tree object.
    relative_means: Series
        Relative mean expression for all genes on every lineage tree branch
    abs_max: int, optional
        Highest allowed value for the absolute average expression of a gene
        along the lineage tree
    gene_mean: float, optional
        Average of the log-normal distribution from which the base gene
        expression values are sampled
    gene_std: float, optional
        Standard deviation of the log-normal distribution from which the base
        gene expression values are sampled

    Returns
    -------
    base_gene_exp: numpy.ndarray
        An array that contains base expression values for each gene
    """
    base_gene_exp = np.zeros(tree.G)

    log_generator = sp.stats.norm(loc=gene_mean, scale=gene_std)

    max_gene_per_branch = max_relat_exp(tree, relative_means)
    max_per_gene = np.max(max_gene_per_branch, axis=1)

    for gene in range(tree.G):
        tmp = np.exp(log_generator.rvs())
        counter = 0
        while tmp * max_per_gene[gene] > abs_max:
            counter = counter + 1
            tmp = np.exp(log_generator.rvs())
        base_gene_exp[gene] = tmp
    return base_gene_exp


def calc_scalings(cells, scale=True, scale_v=0.7):
    """
    Obtain library size factors for each cell.

    Parameters
    ----------
    cells: int
        The number of cells
    scale: bool, optional
        Whether to simulate different scaling factors for each cell
    scale_v: float, optional
        The standard deviation of the library size distribution (log-normal
        distribution around 0)

    Returns
    -------
    scalings: numpy.ndarray
        A library size factor for each cell
    """
    if scale:
        scalings = np.exp(sp.stats.norm.rvs(loc=0., scale=scale_v, size=cells))
    else:
        scalings = np.ones(cells)
    return scalings


def process_timeseries_input(series_points, cells, point_std):
    """
    Process the input of sample_pseudotime_series to make everything the same
    shape.

    Parameters
    ----------
    series_points: list
        The pseudotime sample points for the time series experiment
    cells: int or list
        Either the total number of cells to be sampled (in which case it is
        split equally among all sample points) or the number of cells to be
        sampled at each sample point
    point_std: float, list
        Standard deviation of cell density around each sample point. If it is a
        float, then it is the same for every sample point

    Returns
    -------
    series_points: numpy.ndarray
        The pseudotime sample points for the time series experiment
    cells: numpy.ndarray
        The cells to be sampled at each sample point of the time series
        experiment
    point_std: numpy.ndarray
        The cell density at each sample point of the time series experiment
    """
    no_samples = len(series_points)
    if isinstance(cells, collections.Iterable):
        cells = np.array(cells, dtype=int)
    elif isinstance(cells, numbers.Number):
        cells = np.array([cells / no_samples] * no_samples, dtype=int)

    if isinstance(point_std, collections.Iterable):
        point_std = np.array(point_std, dtype=float)
    elif isinstance(point_std, numbers.Number):
        point_std = np.array([point_std / no_samples] * no_samples, dtype=float)

    if not isinstance(series_points, np.ndarray):
        series_points = np.array(series_points, dtype=int)

    return series_points, cells, point_std


def breadth_first_branches(tree):
    """
    Performs a breadth-first traversal of the tree topology.

    Parameters
    ----------
    tree: Tree
        A lineage tree object.

    Returns
    -------
    bfs: list
        The tree branches in the order of traversal (breadth-first).
    """
    start = tree.root
    levels = {branch: -1 for branch in tree.branches}
    levels[start] = 0
    bfs_topology = bfs_finder(np.array(tree.topology), start)
    for pair in bfs_topology:
        levels[pair[1]] = levels[pair[0]] + 1
    sorted_levels = sorted(levels.items(), key=operator.itemgetter(1))
    bfs = np.array(sorted_levels)[:, 0]
    return bfs


def bfs_finder(graph, start):
    """
    Perform a breadth-first search where the graph is a list of connections.

    Parameters
    ----------
    graph: numpy.ndarray
        A numpy array of shape (N, 2). Every row [a, b] describes a connection from
        branch a to branch b.
    start: int
        The root node from which to start the traversal.

    Returns
    -------
    output: numpy.ndarray
        The input graph sorted by breadth-first traversal order.

    Originally answered by StackOverflow user https://stackoverflow.com/users/2988730
    for question https://stackoverflow.com/questions/50589804.
    """
    sorter = np.argsort(graph[:, 0])
    done = set()
    todo = deque([start])
    output = np.empty_like(graph)
    pos = 0
    while todo:
        key = todo.popleft()
        if key in done:
            continue
        done.add(key)
        left = np.searchsorted(graph[:, 0], key, 'left', sorter)
        if left >= graph.shape[0] or graph[sorter[left], 0] != key:
            continue
        right = np.searchsorted(graph[:, 0], key, 'right', sorter)
        next_pos = pos + right - left
        output[pos:next_pos, :] = graph[sorter[left:right], :]
        todo.extend(output[pos:next_pos, 1])
        pos = next_pos
    return output


def adjust_to_parent(relative_means, current, topology):
    """
    Adds a vector to the relative means of the current branch such that its first
    row is equal to the last row of the relative means of the branch that precedes
    it.

    Parameters
    ----------
    relative_means: Series
        Relative mean expression for all genes on every (currently available) lineage
        tree branch.
    current: int
        The current branch to be adjusted.
    topology: numpy.ndarray
        The tree topology in array format.

    Returns
    -------
    res: numpy.ndarray
        Adjusted relative expression matrix for the current branch.
    """
    parent_loc = (topology[:, 1] == current)
    if not np.any(parent_loc):
        return relative_means[current]
    branch_from = topology[parent_loc][0][0]
    branch_to = current
    adjust_from = relative_means[branch_to]
    adjust_to = relative_means[branch_from]
    res = bifurc_adjust(adjust_from, adjust_to)
    return res


def find_parallel(tree, programs, branch):
    """
    Find all branches that are parallel to the input branch (have same parent branch).

    Parameters
    ----------
    tree: Tree
        A lineage tree object.
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree.
    branch: int
        The branch to examine.

    Returns
    -------
    list
        A list of branches that are parallel to the input branch, including the 
        branch itself.
    """
    parallel = tree.get_parallel_branches()
    for parallels in parallel.values():
        if branch in parallels:
            return np.intersect1d(parallels, list(programs.keys()))
    return [branch, None]
