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




sample_pseudotime_series
sample_density
sample_whole_tree
_sample_data_at_times