#!/usr/bin/env python
# coding: utf-8
"""
This module contains all functions that pertain to the model that generates
UMI count data for the simulated cells.
"""

import numpy as np
import scipy as sp
from scipy.special import gamma as Gamma
from scipy.special import loggamma

def generate_negbin_params(tree, mean_alpha=0.2, mean_beta=1, a_scale=1.5, b_scale=1.5):
    """
    Generate default hyperparameters for the negative binomial distributions
    that are used to simulate UMI count data.

    For motivation of the choice of distribution and default values, please
    refer to the paper.

    Parameters
    ----------
    tree: Tree
        A lineage tree
    mean_alpha: float, optional
        The average alpha value
    mean_beta: float, optional
        The average beta value
    a_scale: float, optional
        The standard deviation for alpha
    b_scale: float, optional
        The standard deviation for beta

    Returns
    -------
    alphas: ndarray
        Alpha values for each gene
    betas: ndarray
        Beta values for each gene
    """
    mu_a = np.log(mean_alpha)
    s2_a = np.log(a_scale)
    mu_b = np.log(mean_beta)
    s2_b = np.log(b_scale)
    alphas = np.exp(sp.stats.norm.rvs(loc=mu_a, scale=s2_a, size=tree.G))
    betas = np.exp(sp.stats.norm.rvs(loc=mu_b, scale=s2_b, size=tree.G)) + 1
    return alphas, betas


def lognegbin(x, theta):
    """
    Alternative formulation of the log negative binomial distribution pmf. scipy
    does not support the extended definition so we have to supply it ourselves.

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
    Alternative formulation of the negative binomial distribution pmf. scipy
    does not support the extended definition so we have to supply it ourselves.

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
    of the original transcripts before amplification. We make the assumption
    that the amplification has no sequence bias and is the same for all
    transcripts.

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
    p = (s2 - m) / s2 if s2 > 0 else 0
    r = (m**2) / (s2 - m) if s2 > 0 else 0
    return p, r


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
