#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy import random
import scipy as sp
from scipy.special import gamma as Gamma
from scipy.special import loggamma
from scipy import stats
import pandas as pd
import random as rng
import sys
import warnings
from collections import defaultdict


def generate_negbin_params(tree, alpha=0.2, beta=1, a_scale=1.5, b_scale=1.5):
    mu_a = np.log(alpha)
    s2_a = np.log(a_scale)
    mu_b = np.log(beta)
    s2_b = np.log(b_scale)
    alphas = np.exp(sp.stats.norm.rvs(loc=mu_a, scale=s2_a, size=tree.G))
    betas = np.exp(sp.stats.norm.rvs(loc=mu_b, scale=s2_b, size=tree.G)) + 1
    return alphas, betas


def printProgress(iteration, total, prefix='', suffix='', decimals=1,
                  barLength=100):
    """
    Call in a loop to create a terminal-friendly text progress bar. Contributed
    by Greenstick on stackoverflow.com/questions/3173320.

    Parameters
    ----------
        iteration: int
            Current iteration.
        total: int
            Total iterations.
        prefix: str, optional
            Prefix string.
        suffix: str, optional
            Suffix string.
        decimals: int, optional
            Positive number of decimals in percent complete.
        barLength: int, optional
            Character length of bar.
    """
    formatStr = "{0:." + str(decimals) + "f}"
    percent = formatStr.format(100 * (iteration / float(total)))
    filledLength = int(round(barLength * iteration / float(total)))
    bar = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)),
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
        x = rng.randrange(k)
        results[x].append(value)
    return results


def lognegbin(x, theta):
    """
    Alternative formulation of the negative binomial distribution pmf.

    scipy does not support the extended definition so we have to supply it
    ourselves.

    Parameters
    ----------
    x: int
        The random variable.
    theta: real array [p, r]
        p is the probability of success of the Bernoulli test and r the number
        of "failures".

    Returns
    -------
    Probability that a discrete random variable is exactly equal to some value.
    """
    p, r = theta
    if p == 0 and r == 0:
        return 0
    else:
        return (loggamma(r + x) + np.log(1 - p) * r + np.log(p) * x -
                (loggamma(r) + loggamma(x + 1)))


def negbin(x, theta):
    """
    Alternative formulation of the negative binomial distribution pmf.

    scipy does not support the extended definition so we have to supply it
    ourselves.

    Parameters
    ----------
    x: int
        The random variable.
    theta: real array [p, r]
        p is the probability of success of the Bernoulli test and r the number
        of "failures".

    Returns
    -------
    Probability that a discrete random variable is exactly equal to some value.
    """
    p, r = theta
    if p == 0 and r == 0:
        return 1 if x == 0 else 0
    else:
        return (Gamma(r + x) * (1 - p)**r * p**x /
                (Gamma(r) * sp.special.factorial(x)))


def get_pr_amp(mu_amp, s2_amp, ksi):
    """
    Calculate parameters the negative binomial that describes the distribution
    of the original transcripts before amplification. We make the (strong)
    assumption that the amplification has no sequence bias and is the same for
    all transcripts.

    Parameters
    ----------
    mu_amp: float
        Mean expression of the amplification.
    s2_amp: float
        Variance of the amplification.
    ksi: int
        Number of initial transcripts present.

    Returns
    -------
    p_amp: float
        The probability of success of the Bernoulli test.
    r_amp: float
        The number of "failures" of the Bernoulli test.
    """
    s2 = ksi * s2_amp
    m = ksi * mu_amp
    p_amp = (s2 - m) / s2 if s2 > 0 else 0
    r_amp = (m**2) / (s2 - m) if s2 > 0 else 0
    return p_amp, r_amp


def get_pr_umi(a, b, m):
    """
    Calculate parameters for my_negbin from the mean and variance of the
    distribution.

    For single cell RNA sequencing data we assume that the distribution of the
    transcripts is described by a negative binomial where the variance s^2
    depends on the mean mu by a relation s^2 = a*mu^2 + b*mu.

    Parameters
    ----------
    a: float
        Coefficient for the quardratic term. Dominates for high mean expression.
    b: float
        Coefficient for the linear term. Dominates for low mean expression.
    m: float
        Mean expression of a gene.

    Returns
    -------
    p: float
        The probability of success of the Bernoulli test.
    r: float
        The number of "failures" of the Bernoulli test.
    """
    s2 = (a * m**2 + b * m)
    p = (s2 - m) / s2 if s2>0 else 0
    r = (m**2) / (s2 - m) if s2>0 else 0
    return p, r


def diffusion(T):
    """
    amortized diffusion process, usually between 0 and 1

    Parameters
    ----------
    T: int
        The length of the diffusion process (T steps).

    Returns
    -------
    W: float array
        A diffusion process of length T.
    """
    V = np.zeros(T)
    W = np.zeros(T)

    W[0] = sp.stats.uniform.rvs()
    V[0] = sp.stats.norm.rvs(loc=0, scale=0.2)

    s_eps = 1 / T
    eta = sp.stats.uniform.rvs()

    for t in range(0, T - 1):
        W[t + 1] = W[t] + V[t]

        epsilon = sp.stats.norm.rvs(loc=0, scale=s_eps)
        # amortize the update
        V[t + 1] = 0.95 * V[t] + epsilon - eta * V[t]

        # quality control: we are not (?) allowed to go below 0. If it happens,
        # reverse and dampen velocity
        # if W[t+1] <= 0:
        #     W[t+1] = W[t]
        #     V[t+1] = -0.2 * V[t]
    return W


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


def simulate_W(T, K, cutoff=0.2, maxloops=100):
    """
    Return K diffusion processes of length T as a matrix.

    The result of simulate_W is complementary to simulate_H.

    W encodes how K modules of coexpressed genes will behave through
    differentiation time. Each matrix W describes one branch of a
    differentiation tree (between two branch points or between a branch point
    and an endpoint). W describes the module in terms of relative expression
    (from 0 to a small positive float, so from a gene not being expressed to a
    gene being expressed at 2x, 3x of its "normal" level).

    After each new diffusion process is added the function checks whether the
    new diffusion correlates with any of the older ones. If the correlation is
    too high (above 0.3 per default), the last diffusion process will be
    replaced with a new one until one is found that does not correlate with any
    other columns of W or a suitable replacement hasn't been found after 100
    tries.

    Obviously this gets more probable the higher the number of components is -
    it might be advisable to change the number of maximum loops allowed or
    the cutoff in order to reduce runtime for a high number of components.

    Parameters
    ----------
    T: int
        The length of the branch of the differentiation tree.
    K: int
        The number of components/modules of coexpression that describe the
        differentiation in this branch of the tree.
    cutoff: float, optional
        Correlation above the cut-off will be considered too much. Should be
        between 0 and 1 but is not explicitly tested.
    maxloops: int, optional
        The maximum number of times the method will try simulating a new
        diffusion process that doesn't correlate with all previous ones in W
        before resetting the matrix and starting over.
    """
    W = np.zeros((K, T))
    k = 0
    loops = 0
    while k < K:
        W[k] = diffusion(T)

        correlates = test_correlation(W, k, cutoff)
        if correlates:
            # repeat and hope it works better this time
            loops += 1
            continue
        else:
            loops = 0
            k += 1

        if loops > maxloops:
            # we tried so hard
            # and came so far
            # but in the end
            # it doesn't even matter
            return simulate_W(T, K, cutoff=cutoff)

    return np.transpose(W)


def _simulate_H_beta(tree, groups, a=2, b=2):
    """
    Return sparse membership matrix for G genes in K groups.

    The result of simulate_H is complementary to simulate_W.

    H encodes how G genes are expressed by defining their membership to K
    expression modules (coded in a matrix W). H could be told to encode
    metagenes, as it contains the information about which genes are coexpressed
    (genes that belong to/are influenced by the same modules). The influence of
    a module on a gene is measured by a number between 0 and 1, drawn from a
    (symmetric, if used with default values) beta distribution.

    Parameters
    ----------
    K: int
        Number of modules.
    G: int
        Number of genes.
    groups: array of int arrays
        groups[i] contains all genes whose expression is influenced by module i
    a: float
        First shape parameter of beta distribution, should be positive.
    b: float
        Second shape parameter of beta distribution, should be positive.
    """
    H = np.zeros((tree.modules, tree.G))
    for k in range(tree.modules):
        for g in groups[k]:
            H[k][g] += sp.stats.beta.rvs(a, b)
    return H


def _simulate_H_gamma(tree, a=0.05):
    coefficients = {}
    K = tree.modules
    G = tree.G
    for branch in tree.branches:
        coefficients[branch] = np.reshape(
            sp.stats.gamma.rvs(a, size=K * G), (K, G))
    return coefficients


def simulate_coefficients(tree, a=0.05, **kwargs):
    coefficients = {}
    if "a" not in kwargs.keys():
        warnings.warn( "No argument 'a' specified in kwargs: using gamma and a=0.05", UserWarning)
        return _simulate_H_gamma(tree, a)
    # if a, b are present: beta distribution
    if "b" in kwargs.keys():
        groups = create_groups(tree.modules, tree.G)
        for branch in tree.time.keys():
            coefficients[branch] = _simulate_H_beta(tree, groups)
    else:
        return _simulate_H_gamma(tree, a=kwargs['a'])
    return coefficients


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


def simulate_branching_data(tree, intra_branch_tol=0.2, inter_branch_tol=0.5, **kwargs):
    if not len(tree.time) == tree.num_branches:
        print("the parameters are not enough for %i branches" %
              tree.num_branches)
        sys.exit(1)

    coefficients = simulate_coefficients(tree, **kwargs)
    programs = simulate_expression_programs(tree, intra_branch_tol)

    # check that parallel branches don't overlap too much
    programs, relative_means = correct_parallel(tree, programs, coefficients, intra_branch_tol, inter_branch_tol)

    # adjust the ends of the relative mean expression matrices
    for b in tree.topology:
        relative_means[b[1]] = bifurc_adjust(relative_means[b[1]], relative_means[b[0]])

    return (pd.Series(relative_means),
            pd.Series(programs),
            pd.Series(coefficients))


def pearson_between_programs(genes, a, b):
    pearson = np.zeros(genes)
    for g in range(genes):
        pearson[g] = sp.stats.pearsonr(a[:, g], b[:, g])[0]
    return pearson


def flat_order(n):
    size = int(n * (n - 1) / 2)
    res = np.zeros((size, 3), dtype=int)
    for i in range(n - 1):
        for j in range(i + 1, n):
            index = int((i * (2 * n - i - 3) / 2 + j - 1))
            res[index] = np.array([index, i, j])
    return res


def simulate_expression_programs(tree, tol):
    programs = {}
    for branch in tree.branches:
        programs[branch] = simulate_W(
            tree.time[branch], tree.modules, cutoff=tol)
    return programs


def calc_relat_means(tree, programs, coefficients):
    relative_means = {}
    for branch in tree.branches:
        relative_means[branch] = np.dot(programs[branch], coefficients[branch])
    return relative_means


def diverging_parallel(branches, programs, genes, tol=0.5):
    indices = flat_order(len(branches))
    res = np.zeros(len(indices), dtype=bool)
    for index, i, j in indices:
        a = branches[i]
        b = branches[j]
        pearson = pearson_between_programs(genes, programs[a], programs[b])
        percent_anticorrelated = sum(pearson < 0) / (genes * 1.0)
        res[index] = percent_anticorrelated > tol
    return res


def correct_parallel(tree, programs, coefficients, intra_branch_tol=0.2, inter_branch_tol=0.5):
    # calculate relative means over tree
    relative_means = calc_relat_means(tree, programs, coefficients)
    # find parallel branches
    parallel = tree.get_parallel_branches()
    # for all parallel branches check if they are diverging.
    # if not, fix.
    for key in parallel:
        diverges = diverging_parallel(
            parallel[key], relative_means, tree.G, tol=inter_branch_tol)
        while not diverges:
            for branch in parallel[key]:
                programs[branch] = simulate_W(tree.time[branch], tree.modules, cutoff=intra_branch_tol)
                relative_means[branch] = np.dot(
                    programs[branch], coefficients[branch])
            diverges = diverging_parallel(
                parallel[key], relative_means, tree.G, tol=inter_branch_tol)
    return programs, relative_means


def sample_data(N, G, tree, sample_times, sample_spreads, c_spread=10,
                alpha=0.3, beta=2, verbose=True, scale_v=0.7):
    # how many cells do we get from each time point?
    # we want to have ~N cells in the end, so divide N by the sample point size
    n_cells = (N * 1.) / len(sample_times)
    cell_spread = np.sqrt(c_spread)
    # we sample n_cells in each sample point, but we lose some for technical
    # reasons
    cells = sp.stats.norm.rvs(loc=n_cells, scale=cell_spread,
                              size=len(sample_times))
    cells = np.abs(cells)
    cells = cells.astype(int)

    timestamps = []

    # sample from the total time. It starts at 0 and ends at total_time
    total_time = tree.get_max_time()
    for tp, n in zip(sample_times, cells):
        # print(tp, n, total_time)
        timestamps.extend(collect_timestamps(tp, n, total_time,
                          t_s=sample_spreads))
    return sample_data_at_times(cells, tree, timestamps, alpha=alpha,
                                beta=beta, verbose=verbose, scale_v=scale_v)


def sample_data_at_times(tree, pseudotime, alpha=0.3, beta=2,
                         scale_v=0.7, branches=None, verbose=True, scale=True):
    N = len(pseudotime)
    if np.shape(alpha) == ():
        alpha = [alpha] * tree.G
    if np.shape(beta) == ():
        beta = [beta] * tree.G
    scalings = calc_scalings(N, scale, scale_v)
    branches = pick_branches(tree, pseudotime)

    X = draw_counts(tree, pseudotime, branches, scalings, alpha, beta, verbose)

    return X, pseudotime, branches, scalings


def restricted_simulation(t, alpha=0.2, beta=3, mult=2, gene_loc=0.8, gene_s=1):
    sample_time = np.arange(0, t.get_max_time())
    gene_scale = np.exp(sp.stats.norm.rvs(loc=0.8, scale=1, size=t.G))
    Ms = {}
    while not are_lengths_ok(Ms):
        uMs, Ws, Hs = simulate_branching_data(t, a=0.05)
        for i in t.branches:
            Ms[i] = np.exp(uMs[i]) * gene_scale

    t.add_genes(Ms)
    alphas, betas = generate_negbin_params(t, alpha=alpha, beta=beta)

    X, labs, brns, scalings = sample_data_at_times(t, sample_time, alphas, betas)
    return X, labs, brns, scalings


def sample_density(tree, N, alpha=0.3, beta=2, scale_v=0.7, verbose=True, scale=True):
    bt = tree.branch_times()

    possible_pt = [np.arange(bt[b][0], bt[b][1] + 1) for b in tree.branches]
    possible_branches = [[b] * tree.time[b] for b in tree.branches]
    probabilities = [tree.density[b] for b in tree.branches]

    # make numpy arrays and flatten lists
    probabilities = np.array(probabilities).flatten()
    possible_pt = np.array(possible_pt).flatten()
    possible_branches = np.array(possible_branches).flatten()

    # select according to density and take the selected elements
    sample = random.choice(np.arange(len(probabilities)), size=N, p=probabilities)
    sample_time = possible_pt[sample]
    sample_branches = possible_branches[sample]

    X, labels, branches, scalings = sample_data_at_times(tree, sample_time,
                            alpha=alpha, beta=beta, branches=sample_branches)
    return X, labels, branches, scalings


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


def collect_timestamps(timepoint, n, Tmax, t_s=4):
    """
    Samples around a timepoint.

    This function simulates how a differentiating cell population would look
    at a time point when the mean differentiation stage would be at pseudotime
    t, with cells be spread around t. Since we are sampling a discrete space
    we need to round the samples. We can control the asynchronicity via the
    variance of the normal distribution.

    Parameters
    ----------
    timepoint: int
        The pseudotime point that represents the current mean differentiation
        stage of the population.
    n: int
        How many cells to sample.
    Tmax: int
        All time points that exceed the differentiation duration will be
        mapped to the end of the differentiation.
    t_s: float, optional
        Variance of the normal distribution we use to draw pseudotime points.
        In the experiment metaphor this parameter controls synchronicity.

    Returns
    -------
    timestamps: int array
        Pseudotime points around <timepoint>.
    """
    timestamps = sp.stats.norm.rvs(loc=timepoint, scale=t_s, size=n)
    timestamps = timestamps.astype(int)
    timestamps[timestamps < 0] = 0
    timestamps[timestamps >= Tmax] = Tmax - 1
    return timestamps


class my_negbin(sp.stats.rv_discrete):
    """
    Class definition for the alternative negative binomial pmf so that we can
    sample it using rvs().
    """
    def _pmf(self, x, p, r):
        theta = [p, r]
        res = np.exp(lognegbin(x, theta))
        res = np.real(res)
        return res.astype("float")


class sum_negbin(sp.stats.rv_discrete):
    """
    Class definition for the convoluted negative binomial pmf that describes
    non-UMI data.
    """
    def _pmf(self, x, mu_amp, s_amp, p, r):
        theta = [p, r]
        ksis = np.arange(2 * int(x) + 3)
        res = 0

        for ksi in ksis:
            p_amp, r_amp = get_pr_amp(mu_amp, s_amp, ksi)
            theta_amp = [p_amp, r_amp]
            # print(theta_amp)
            # print(negbin(x, theta_amp), negbin(ksi, theta))
            tmp = lognegbin(x, theta_amp) + lognegbin(ksi, theta)
            res += np.real(np.exp(tmp))
        return res.astype("float")


def are_lengths_ok(Ms=None, abs_max=800, rel_dif=0.3):
    """
    Performs a quality check for the matrix of average gene expressions for all
    branches.

    Specifically, it is required that no gene has a mean expression value above
    a certain cutoff and that no branch has a maximum average expression that
    is an order of magnitude above the other branches.

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
    # worked out kinda ok
    if not Ms:
        return False

    max_crit = False
    rel_crit = False

    maxes = np.array([np.max(np.ptp(x, axis=0)) for x in Ms.values()])
    if np.all(maxes < abs_max):
        max_crit = True

    if np.min(maxes) / np.max(maxes) > rel_dif:
        rel_crit = True
    # print(max_crit, rel_crit)

    # means = [np.mean(np.ptp(x, axis=0)) for x in Ms]
    return (max_crit and rel_crit)


def sample_tree_points(tree):
    timezone = tree.populate_timezone()
    assignments = assign_branches(tree.branch_times(), timezone)
    pseudotime = list()
    branches = list()

    for n, a in enumerate(timezone):
        start = a[0]
        end = a[1] + 1
        length = end - start
        for b in assignments[n]:  # for all possible branches in timezone a
            pseudotime.extend(np.arange(start, end))
            branches.extend([b] * length)
    return pseudotime, branches


def calc_scalings(N, scale=True, scale_v=0.7):
    if scale:
        scalings = np.exp(sp.stats.norm.rvs(loc=0., scale=scale_v, size=N))
    else:
        scalings = np.ones(N)
    return scalings


def draw_counts(tree, pseudotime, branches, scalings, alpha, beta, verbose):
    N = len(branches)
    X = np.zeros((N, tree.G))
    custm = my_negbin()

    for n, t, b in zip(np.arange(N), pseudotime, branches):
        T_off = tree.branch_times()[b][0]
        M = tree.means[b]

        for g in range(tree.G):
            try:
                mu = M[t - T_off][g] * scalings[n]
            except IndexError:
                print("IndexError for g=%d, t=%d, T_off=%d in branch %s" % (g, t, T_off, b))
                mu = M[-1][g] * scalings[n]
            p, r = get_pr_umi(a=alpha[g], b=beta[g], m=mu)
            X[n][g] = custm.rvs(p, r)

        if verbose:
            printProgress(n, len(branches))
    return X


def sample_whole_tree(n_factor, tree, alpha=0.3, beta=2, scale_v=0.7,
                      verbose=True, scale=True):
    if np.shape(alpha) == ():
        alpha = [alpha] * tree.G
    if np.shape(beta) == ():
        beta = [beta] * tree.G

    pseudotime, branches = sample_tree_points(tree)
    N = len(branches) * n_factor

    scalings = calc_scalings(N, scale, scale_v)
    branches = np.repeat(branches, n_factor)
    pseudotime = np.repeat(pseudotime, n_factor)

    X = draw_counts(tree, pseudotime, branches, scalings, alpha, beta, verbose)

    return X, pseudotime, branches, scalings
