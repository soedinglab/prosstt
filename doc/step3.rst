Sample cells from tree
======================

PROSSTT offers a number of ways to sample cells from the lineage tree:

* sample all points on lineage tree

* sample a pseudotime series

* sample density

* sample particular points of the lineage tree

Sample whole tree
-----------------

In this mode, all possible pairs of pseudotime and branch :math:`t, b` on the tree are sampled a number of times. This mode creates simulations that are easier for trajectory inference algorithms to reconstruct, since all points of the lineage tree are available, i.e. there are no gaps that the algorithm must fill on its own.

To sample every point on the tree 3 times:

```
from prosstt import simulation as sim
sim.sample_whole_tree(tree, 3)
```

In the same vein, PROSSTT includes a function to speed up the simulation process; it will generate expression programs, average gene expression and sample the whole tree in one step using default parameters. This makes it easy for users to get a first taste of PROSSTT without having to decide on parameter values or familiarizing themselves with the software first (see the `minimal example <https://github.com/soedinglab/prosstt/blob/master/examples/minimal_example.ipynb>` in the PROSSTT github repo).

Sample pseudotime series
------------------------

In this mode, PROSSTT simulates a time series experiment given the lineage tree. In a time series experiment a cell population is sampled in different time points along its developmental trajectory. Owing to the asynchrony of cellular differentiation, this means that in every sample (at every time point) there is a mixture of cells from different developmental stages.

In PROSSTT, progress through differentiation is measured by pseudotime. Therefore, the mixture of cells from different developmental stages corresponds to a mixture of cells with different pseudotime. Given a set of sample points, PROSSTT samples from a Gaussian distribution around them with a given standard deviation. For each sampled pseudotime, it picks one of the possible branches randomly.

An `example <https://github.com/soedinglab/prosstt/blob/master/examples/sample_pseudotime_series.ipynb>` and the corresponding plots are available at the github repository.

Sample density
--------------

An additional way to sample cells is to impose a certain density on the lineage tree. Consider transdifferentiations like `the one described by Treutlein *et al.* in Nature<https://www.nature.com/articles/nature18323>`. Most cells will follow one path, but some could progress into a different cell fate than the intended. The difference in efficiency or in abundance of each branch can be represented by a different density.

Differences in differentiation speed can also be represented by different densities. Imagine a differentiation where cells progress quickly to the cell fate decision point and then linger there, after which they progress quickly to their respective cell fates. For PROSSTT, this would mean having a higher density around the branch point and a lower one towards the end points of the lineage tree, as demonstrated in `another example notebook<https://github.com/soedinglab/prosstt/blob/master/examples/density_sampling.ipynb>`.

Combinations and manual sampling
--------------------------------

Density sampling and pseudotime series can be combined; after determining the density of each branch, users can then sample from the tree in a pseudotime series experiment (for an example, see the corresponding `jupyter notebook<https://github.com/soedinglab/prosstt/blob/master/examples/combined_sampling.ipynb>`)

Finally, if users are not covered by the available sampling options and want to create something more elaborate, there is always the option to sample pairs of time points and branch assignments from the lineage tree and let PROSSTT simulate the cells for this input.

```
from prosstt import tree
from prosstt import simulation as sim

newick_string = "(A:50,B:50)C:50;"
G = 500
lineage = tree.Tree.from_newick(newick_string, genes=G)
sample_times = np.array([0, 25, 49, 51, 149, 149])
sample_branches = np.array(['C', 'C', 'C', 'B', 'A', 'A'])
sim._sample_data_at_times(lineage, sample_times, branches=sample_branches)
```