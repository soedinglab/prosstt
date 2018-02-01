#!/usr/bin/env python
# coding: utf-8

import numpy as np

def parse_newick(newick_tree, def_time):
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