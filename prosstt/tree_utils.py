#!/usr/bin/env python
# coding: utf-8
"""
This module contains utility functions for the Tree class.
"""

import numpy as np

def parse_newick(newick_tree, def_time):
    """
    Function that translates a Newick tree to a prosstt Tree object.

    Parameters
    ----------
    newick_tree: Newick tree
        A Newick tree object as returned by the loads function of the newick
        package
    def_time: int
        The default pseudotime length of each tree branch

    Returns
    -------
    topology: list
        A list of connected branch names
    time: dict
        Contains the pseudotime length of each branch
    branches: list
        The name of each branch
    branch_points: int
        The number of branch points in the lineage tree
    root: str
        The name of the root of the lineage tree
    """
    topology = list()
    time = {}
    branches = 0
    branch_points = 0
    root = None
    for node in newick_tree[0].walk():
        branches += 1
        if node.length == 0:
            time.update({node.name: def_time})
        else:
            time.update({node.name: np.int(node.length)})

        if not node.descendants:
            continue
        else:
            branch_points += 1
            for descendant in node.descendants:
                topology.append([node.name, descendant.name])

        if node.ancestor is None:
            root = node.name
    return topology, time, branches, branch_points, root
