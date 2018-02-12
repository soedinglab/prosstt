#!/usr/bin/env python
# coding: utf-8

import collections
from collections import defaultdict
import numbers
import sys
import shutil

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
    bar_length = shutil.get_terminal_size((80, 20))[0]
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
    results: list of int lists


    contributed by kennytm on stackoverflow.com/questions/3760752
    """
    results = [[] for i in range(k)]
    for value in iterable:
        x = random.randint(k+1)
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
        p = sp.stats.pearsonr(W[k], W[i])
        if p[0] > cutoff:
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
    

def pearson_between_programs(genes, a, b):
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
    pearson = np.zeros(genes)
    for g in range(genes):
        pearson[g] = sp.stats.pearsonr(a[:, g], b[:, g])[0]
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
        A lineage tree
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree
    coefficients: ndarray
        Array that contains the contribution weight of each expr. program for
        each gene
    """
    relative_means = {}
    for branch in tree.branches:
        relative_means[branch] = np.dot(programs[branch], coefficients[branch])
    return relative_means


def diverging_parallel(branches, programs, genes, tol=0.5):
    """
    Calculate if the expression programs in all pairs of parallel branches are
    diverging enough to make the branches distinguishable from each other.

    Parameters
    ----------
    branches: list
        A list of pairs of parallel branches
    programs: Series
        Relative expression for all expression programs on every branch of the
        lineage tree
    genes: int
        The number of genes included in the lineage tree
    tol: float, optional
        The percentage of genes that must have anticorrelated expression
        patterns over pseudotime in order for the branches to be considered
        diverging
    
    Returns
    -------
    diverging: ndarray
        A list of the boolean values: whether each pair of parallel branches
        diverges or not
    """
    indices = flat_order(len(branches))
    diverging = np.zeros(len(indices), dtype=bool)
    for index, i, j in indices:
        a = branches[i]
        b = branches[j]
        pearson = pearson_between_programs(genes, programs[a], programs[b])
        percent_anticorrelated = sum(pearson < 0) / (genes * 1.0)
        diverging[index] = percent_anticorrelated > tol
    return diverging


def assign_branches(branch_times, timezone):
    """
    Assigns a branch to every timezone.

            -- T[1]------
    -T[0]--|          -- T[3]------
            -- T[2]--|
                      -- T[4]-
    timezones:
    ---0---|----1----|-2-|-3--|-4--

    tz   br
    0    0
    1    1,2
    2    1,3,4
    3    3,4
    4    3

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
    for i in range(len(timezone)):
        z = timezone[i]
        for k in branch_times:
            b = branch_times[k]
            if belongs_to(z, b):
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
        branches[n] = pick_branch(t, timezone, assignments)
    return branches
    

def pick_branch(t, timezones, assignments):
    """
    Picks one of the possible branches for a cell at a given time point.

    Parameters
    ----------
    t: int
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
    b = -1
    for i in range(len(timezones)):
        branch = timezones[i]
        if t >= branch[0] and t <= branch[1]:
            b = i
            break
    possibilities = assignments[b]
    try:
        return random.choice(possibilities)
    except IndexError:
        print(t)
        print(timezones)
        print(assignments)


def are_lengths_ok(average_expression, abs_max=800, rel_dif=0.3):
    """
    Performs a quality check for the matrix of average gene expressions for all
    branches.

    Specifically, it is required that no gene has a mean expression value above
    a cutoff, and that no branch has a maximum average expression that is a
    (specified) order of magnitude above the other branches.

    Parameters
    ----------
    average_expression: Series
        The average expression of each gene over pseudotime for each branch
    abs_max: Numbers.number, optional
        Maximum average expression allowed for a gene
    rel_dif: float, optional
        Maximum value of the ratio between the smallest and greatest maximum
        average gene expression value of all genes in all branches.

    Returns
    -------
    (): bool
        Whether both conditions are fulfilled
    """
    # worked out kinda ok
    if not average_expression:
        return False

    max_crit = False
    rel_crit = False

    maxes = np.array([np.max(np.ptp(x, axis=0)) for x in average_expression.values()])
    if np.all(maxes < abs_max):
        max_crit = True

    if np.min(maxes) / np.max(maxes) > rel_dif:
        rel_crit = True

    return max_crit and rel_crit


def calc_scalings(N, scale=True, scale_v=0.7):
    """
    Obtain library size factors for each cell.

    Parameters
    ----------
    N: int
        The number of cells
    scale: bool, optional
        Whether to simulate different scaling factors for each cell
    scale_v: float, optional
        The standard deviation of the library size distribution (log-normal
        distribution around 0)
    
    Returns
    -------
    scalings: ndarray
        A library size factor for each cell
    """
    if scale:
        scalings = np.exp(sp.stats.norm.rvs(loc=0., scale=scale_v, size=N))
    else:
        scalings = np.ones(N)
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
    series_points: ndarray
        The pseudotime sample points for the time series experiment
    cells: ndarray
        The cells to be sampled at each sample point of the time series
        experiment
    point_std: ndarray
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
