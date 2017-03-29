#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import numpy as np

class Tree(object):
	'topology of the differentiation tree'

	# default values for when the user is not decided
	def_time = 40
	def_comp = 15	

	def __init__(self, topology, time, branches, branch_points, modules):
		self.topology = topology
		self.time = time
		self.branches = branches
		self.branch_points = branch_points
		self.modules = modules
		self.means = None

	@classmethod
	def from_topology(cls, topology):
		# information about branches/branchpoints is in the topology
		branches, branch_points = analyze_topology(topology)

		# now we can create everything else:
		time = [def_time] * branches
		modules = [def_comp] * branches

		return cls(topology, time, branches, branch_points, modules)		


	def analyze_topology(topology):
		# for [[0,1], [0,2], [2,3], [2,4], [2,5]]
		# #diff numbers at position 1 = #branchpoints
		# #diff numbers in topology = #branches

		all_branches = [branch for pair in topology for branch in pair]
		unique_branches = set(all_branches)

		starts = [pair[0] for pair in topology]
		unique_starts = set(starts)

		return len(unique_branches), len(unique_starts)


	def add_genes(self, Ms):
		self.means = Ms


	def get_max_time(self):
		"""
		Given the lengths of the branches calculate the maximum duration described by the tree.

		Parameters
		----------
		Ts: int array
			The duration of the branches of the differentiation tree. It is assumed that branches in the same timezone have identical length.
		branch_points: int
			The number of branch points of the differentiation tree.

		Returns
		-------
		total: int
			Total pseudotime duration of the differentiation tree.
		"""

		# make a list of all possible paths through the tree
		# and calculate the length of those paths, then keep max
		# first take topology and make it a dictionary:
		treedict = self.dictionify()

		# find paths to leaves in dict:
		tree_paths = np.array(self.paths(0, treedict))

		total_lengths = np.zeros(len(tree_paths))

		for i in range(len(tree_paths)):
			total_lengths[i] = np.sum(np.array(self.time)[tree_paths[i]])

		return np.max(total_lengths)


	def dictionify(self):
		treedict = defaultdict(list)
		for b in self.topology:
			treedict[b[0]].append(b[1])
		return(treedict)


	# find all paths to leaves from given start point
	def paths(self, key, tree):
		if not tree[key]:
			return [[key]]
		else:
			rooted_paths = []
			root = [key]
			for v in tree[key]:
				usable = self.paths(v, tree)
				for path in usable:
					rooted_paths.append(root + path)
			return rooted_paths

	
	def populate_timezone(self):
		"""
		Returns an array that assigns pseudotime to time zones.

		This function first determines the timezones by considering the length of the branches and then assigns a timezone to each pseudotime range. E.g. for Ts = [25, 25, 25] we would have timezone[0:24] = 0, timezone[25:49] = 1.

		Parameters
		----------
		total_time: int
			The total pseudotime duration of the differentiation tree.
		Ts: int array
			The duration of the branches of the differentiation tree.
		branch_points: int
			The number of branch points of the differentiation tree.

		Returns
		-------
		timezone: int array
			Array of length total_time, contains the timezone information for each pseudotime point.
		updated_Ts: int array
			Converts from relative time to absolute time: given Ts=[25,25,25,25,25] branch 0 starts at pseudotime 0, but branches 1 and 2 start at pseudotime 25 and branches 3,4 at pseudotime 50.
		"""

		res = []
		dtree = self.dictionify()
		ntime = np.array(self.time)
		tpaths = self.paths(0, dtree)
		stacks = [self.morph_stack(ntime[np.array(x)].tolist()) for x in tpaths]


		while stacks:
			lpaths = len(stacks)
			curr = [stacks[i][0] for i in range(lpaths)]
			starts = np.array([curr[i][0] for i in range(lpaths)])
			ends = np.array([curr[i][1] for i in range(lpaths)])

			if all(ends == np.max(ends)):
				res.append([np.max(starts), np.max(ends)-1])
				for s in stacks:
					s.pop(0)
			else:
				res.append([np.max(starts), np.min(ends)-1])
				for s in stacks:
					if s[0][1] != np.min(ends):
						s.insert(1, [np.min(ends), s[0][1]])
					s.pop(0)

			newstacks = [x for x in stacks if x]
			stacks = newstacks
		return res


	def branch_times(self):
		branch_time = defaultdict(list)
		
		branch_time[0] = [0, self.time[0]-1]
		for b in self.topology:
			# the start time of b[1] is the end time of b[0]
			b0_end = branch_time[b[0]][1]
			branch_time[b[1]] = [b0_end+1, b0_end+self.time[b[1]]]
		return branch_time


	def morph_stack(self, stack):
		curr = 0
		prev = 0
		for i in range(len(stack)):
			curr = stack[i]
			stack[i] = [prev, prev + stack[i]]
			prev = curr + prev
		return stack