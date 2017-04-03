#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import numpy as np


class Tree(object):
    'topology of the differentiation tree'

    # default values for when the user is not decided
    def_time = 40
    def_comp = 15

    def __init__(self, topology, time, branches, branch_points, modules, G):
        self.topology = topology
        self.time = time
        self.branches = branches
        self.branch_points = branch_points
        self.modules = modules
        self.G = G
        self.means = None

    @classmethod
    def from_topology(cls, topology):
        """
        Alternative constructor that creates a default tree object from a
        given topology.

        Parameters
        ----------
        topology: int array
            An array that describes which branches are connected to each other.
            [[0, 1], [0, 2]] describes a single bifurcation where branch 0
            is connected with branches 1 and 2.

        Returns
        -------
        Tree
            An object of the Tree class with default branch lengths and branch
            times.
        """
        # information about branches/branchpoints is in the topology
        branches, branch_points = cls.analyze_topology(cls, topology)

        # now we can create everything else:
        time = [cls.def_time] * branches
        modules = [cls.def_comp] * branches

        return cls(topology, time, branches, branch_points, modules)

    def analyze_topology(self, topology):
        """
        Identifies the number of branch points and branches described by the
        input topology.

        Parameters
        ----------
        topology: int array
            The topology of a tree. Each entry in the array is a pair of IDs for
            two different branches. The branch described by the second ID
            follows the branch described by the first.

        Returns
        -------
        branches: int
            The number of different branch IDs in the input topology.
        branch_points: int
            The number of branch points in the topology.
        """
        # for [[0,1], [0,2], [2,3], [2,4], [2,5]]
        # #diff numbers at position 1 = #branchpoints
        # #diff numbers in topology = #branches
        all_branches = [branch for pair in topology for branch in pair]
        unique_branches = set(all_branches)

        starts = [pair[0] for pair in topology]
        unique_starts = set(starts)

        return len(unique_branches), len(unique_starts)

    def add_genes(self, Ms):
        """
        Sets the average gene expression trajectories of genes for all branches
        after performing a sanity check.

        Parameters
        ----------
        Ms: float array
            A list of arrays. List length is #branches and array Ms[b] has the
            dimensions time_length[b], G.
        """
        # sanity check of dimensions so that in case a user messes up there is
        # no cryptic IndexOutOfBounds exception they have to trace.
        if not len(Ms) == self.branches:
            msg = "The number of arrays in Ms must be equal to the number of \
                   branches in the topology"
            raise ValueError(msg)
        for i, m in enumerate(Ms):
            if not m.shape == (self.time[i], self.G):
                msg = "Branch " + str(i) + " was expected to have a shape " \
                      + str((self.time[i], self.G)) + " and instead is " \
                      + str(m.shape)
                raise ValueError(msg)

        self.means = Ms

    def get_max_time(self):
        """
        Given the lengths of the branches calculate the maximum duration
        described by the tree.

        Returns
        -------
        total: int
            Total pseudotime duration of the differentiation tree.
        """
        # find paths to leaves in dict:
        tree_paths = np.array(self.paths(0))

        total_lengths = np.zeros(len(tree_paths))

        for i in range(len(tree_paths)):
            total_lengths[i] = np.sum(np.array(self.time)[tree_paths[i]])

        return np.max(total_lengths)

    def dictionify(self):
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

    def paths(self, key):
        """
        Finds all paths from a given start point to the leaves.

        Parameters
        ----------
        key: int
            The starting point.

        Returns
        -------
        rooted_paths: int array
            An array that contains all paths from the starting branch to all
            branches that don't bifurcate.
        """
        # make a list of all possible paths through the tree
        # and calculate the length of those paths, then keep max
        # first take topology and make it a dictionary:
        treedict = self.dictionify()
        if not treedict[key]:
            return [[key]]
        else:
            rooted_paths = []
            root = [key]
            for v in treedict[key]:
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
        ntime = np.array(self.time)
        tpaths = self.paths(0)
        stacks = [self.morph_stack(ntime[np.array(x)].tolist()) for x in tpaths]

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

        branch_time[0] = [0, self.time[0] - 1]
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
