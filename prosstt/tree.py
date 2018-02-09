#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import numpy as np
import pandas as pd
import newick
from prosstt import tree_utils as tu

class Tree(object):
    'topology of the differentiation tree'

    # default values for when the user is not decided
    def_time = 40
    def_comp = 15
    def_genes = 500

    def __init__(self, topology=[["A", "B"], ["A", "C"]],
                 time={"A":def_time, "B":def_time, "C":def_time},
                 num_branches=3,
                 branch_points=1,
                 modules=def_comp,
                 G=def_genes,
                 density=None,
                 root=None):
        self.topology = topology
        self.time = pd.Series(time, name="time")
        self.num_branches = num_branches
        self.branch_points = branch_points
        self.modules = modules
        self.G = G
        self.means = None
        self.branches = list(time.keys())

        if root is None:
            self.root = self.branches[0]
        else:
            self.root = root

        if density is None:
            self.density = self.default_density()
        else:
            self.density = density

    @staticmethod
    # assert branch_points>0
    def gen_random_topology(branch_points):
        n = branch_points
        b = 2 * n + 1
        seeds = [0]
        avail = list(reversed(range(1, b)))
        res = []
        while avail:
            root = np.random.choice(seeds)
            a = avail.pop()
            b = avail.pop()
            res.append([root, a])
            res.append([root, b])
            seeds.append(a)
            seeds.append(b)
            seeds.remove(root)
        return res

    @classmethod
    def from_newick(cls, newick_tree,
                    modules=10,
                    genes=200,
                    density=None):
        """
        Generate a lineage tree from a Newick-formatted string.
        """
        tree = newick.loads(newick_tree)
        top, time, branches, br_points, root = tu.parse_newick(tree, cls.def_time)
        tree = Tree(top, time, branches, br_points, modules, genes, density, root)
        return tree

    @classmethod
    def from_random_topology(cls, branch_points, time, modules, genes):
        """
        Generate a random binary tree topology given a number of branch points.
        """
        topology = Tree.gen_random_topology(branch_points)
        branches = len(np.unique(topology))
        return cls(topology, time, branches, branch_points, modules, genes)

    def default_density(self):
        """
        Initializes the density with a uniform distribution (every cell has the
        same probability of being picked. This is in case the users want to use
        the density sampling function.
        """
        total_time = 0
        density = {}
        for v in self.time.values:
            total_time += v

        for k in self.time.keys():
            density[k] = np.array([1. / total_time] * np.int(self.time[k]))
        return density

    def add_genes(self, Ms):
        """
        Sets the average gene expression trajectories of genes for all branches
        after performing a sanity check.

        Parameters
        ----------
        Ms: float array
            A list of arrays. List length is #branches and array Ms[b] has the
            dimensions time[b], G.
        """
        # sanity check of dimensions so that in case a user messes up there is
        # no cryptic IndexOutOfBounds exception they have to trace.
        if not len(Ms) == self.num_branches:
            msg = "The number of arrays in Ms must be equal to the number of \
                   branches in the topology"
            raise ValueError(msg)

        for branch in Ms:
            mean = Ms[branch]
            if not mean.shape == (self.time[branch], self.G):
                msg = "Branch " + branch + " was expected to have a shape " \
                        + str((self.time[branch], self.G)) + " and instead is " \
                        + str(mean.shape)
                raise ValueError(msg)

        self.means = Ms

    def set_density(self, dens):
        """
        Sets the density as a function of the pseudotime and the branching. If
        N points from the tree were picked randomly, then the density is the
        probability of a pseudotime point in a certain branch being picked.

        Parameters
        ----------
        dens: float array
            A list of arrays. List length is #branches and len(dens[b]) equals
            time[b].
        """
        if not len(dens) == self.branches:
            msg = "The number of arrays in dens must be equal to the number \
                  of branches in the topology"
            raise ValueError(msg)
        for i, d in enumerate(dens):
            if not len(d) == self.time[i]:
                msg = "Branch " + str(i) + " was expected to have a length " \
                      + str((self.time[i], self.G)) + " and instead is " \
                      + str(d.shape)
                raise ValueError(msg)
        self.density = dens

    def get_max_time(self):
        """
        Calculate the maximum pseudotime duration possible for the tree.

        Returns
        -------
        start: str
            Name of the starting node.
        """
        # find paths to leaves in dict:
        tree_paths = self.paths(self.root)

        total_lengths = np.zeros(len(tree_paths))

        for i, path in enumerate(tree_paths):
            path_length = [self.time[branch] for branch in path]
            total_lengths[i] = np.sum(path_length)

        return int(np.max(total_lengths))

    def as_dictionary(self):
        """
        Converts the tree topology to a dictionary where the ID of every branch
        points to the branches that bifurcate from it.

        Returns
        -------
        dict
            The topology of the tree in dictionary form.
        """
        treedict = defaultdict(list)
        for b in self.topology:
            treedict[b[0]].append(b[1])
        return(treedict)

    def paths(self, start):
        """
        Finds all paths from a given start point to the leaves.

        Parameters
        ----------
        start: str
            The starting point.

        Returns
        -------
        rooted_paths: int array
            An array that contains all paths from the starting point to all
            tree leaves.
        """
        # make a list of all possible paths through the tree
        # and calculate the length of those paths, then keep max
        # first take topology and make it a dictionary:
        treedict = self.as_dictionary()
        if not treedict[start]:
            return [[start]]
        else:
            rooted_paths = []
            root = [start]
            for v in treedict[start]:
                usable = self.paths(v)
                for path in usable:
                    rooted_paths.append(root + path)
            return rooted_paths

    def populate_timezone(self):
        """
        Returns an array that assigns pseudotime to time zones.

        This function first determines the timezones by considering the length
        of the branches and then assigns a timezone to each pseudotime range.
        E.g. for Ts = [25, 25, 30] we would have timezone[0:24] = 0,
        timezone[25:49] = 1, timezone[50:54] = 2.

        Returns
        -------
        timezone: int array
            Array of length total_time, contains the timezone information for
            each pseudotime point.
        updated_Ts: int array
            Converts from relative time to absolute time: given
            Ts=[25,25,25,25,25] branch 0 starts at pseudotime 0, but branches 1
            and 2 start at pseudotime 25 and branches 3,4 at pseudotime 50.
        """
        res = []
        tpaths = self.paths(self.root)
        stacks = [self.morph_stack(self.time[x].tolist()) for x in tpaths]

        while stacks:
            lpaths = len(stacks)
            curr = [stacks[i][0] for i in range(lpaths)]
            starts = np.array([curr[i][0] for i in range(lpaths)])
            ends = np.array([curr[i][1] for i in range(lpaths)])

            if all(ends == np.max(ends)):
                res.append([np.max(starts), np.max(ends) - 1])
                for s in stacks:
                    s.pop(0)
            else:
                res.append([np.max(starts), np.min(ends) - 1])
                for s in stacks:
                    if s[0][1] != np.min(ends):
                        s.insert(1, [np.min(ends), s[0][1]])
                    s.pop(0)

            newstacks = [x for x in stacks if x]
            stacks = newstacks
        return res

    def branch_times(self):
        """
        Calculates the pseudotimes at which branches start and end.

        Returns
        -------
        branch_time: dict
            Dictionary that contains the start and end time for every branch.

        Examples
        --------
        >>> from prosstt.tree import Tree
        >>> t = Tree.from_topology([[0,1], [0,2]])
        >>> t.branch_times()
        defaultdict(<class 'list'>, {0: [0, 39], 1: [40, 79], 2: [40, 79]})
        """
        branch_time = defaultdict(list)

        branch_time[self.root] = [0, self.time[self.root] - 1]
        for b in self.topology:
            # the start time of b[1] is the end time of b[0]
            b0_end = branch_time[b[0]][1]
            branch_time[b[1]] = [b0_end + 1, b0_end + self.time[b[1]]]
        return branch_time

    # stacks = [self.morph_stack(ntime[np.array(x)].tolist()) for x in tpaths]
    def morph_stack(self, stack):
        """
        The pseudotime start and end of every branch in a path. Very similar to
        branch_times().

        Parameters
        ----------
        stack: int array
            The pseudotime length of all branches that make up a path in the
            tree (from the origin to a leaf).

        Returns
        -------
        stack: list of 2D arrays
            The pseudotime start and end of every branch in the path.
        """
        curr = 0
        prev = 0
        for i in range(len(stack)):
            curr = stack[i]
            stack[i] = [prev, prev + stack[i]]
            prev = curr + prev
        return stack

    def get_parallel_branches(self):
        top_array = np.array(self.topology)
        parallel = {}
        for branch in np.unique(top_array[:, 0]):
            matches = top_array[:, 0] == branch
            parallel[branch] = top_array[matches, 1]
        return parallel
