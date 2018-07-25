Generate tree
=============

The aim of this step is to generate a ``Tree`` object. For illustration purposes, we will try to create the lineage tree for a differentiation with a single cell fate decision (single bifurcation) where one side of the differentiation lasts longer than the other. We want the simulated data to have a total of 1000 genes and be regulated by 15 expression programs::

                  100
           70   ---B------
        ---A---|         
                ---C--
                  60

A single bifurcation, where branch A is connected to branches B and C. The pseudotime lengths of the branches are A:70, B:100 and C:60, respectively. No part of the lineage tree is populated more densely than others - that means that if we picked a cell from the tree in a random manner, all positions are equally probable. 

This can happen in various ways:

Manual definition
-----------------

The ``Tree()`` constructor can be called directly with a topology that describes branch connectivity, a pseudotime length for each branch, the total number of branches and branch points, as well as values for the number of genes and expression programs (modules) (see :ref:`step2-label`)::

    from prosstt.tree import Tree
    t = Tree(topology=[["A", "B"], ["A", "C"]],
             time={"A":70, "B":100, "C":60},
             num_branches=3,
             branch_points=1,
             modules=15,
             G=1000,
             density=None,
             root="A")

If ``density=None`` is passed (or if no value for the parameter is given), PROSSTT will automatically calculate a uniform density for the lineage tree.

From a Newick-formatted string
------------------------------

A Newick-formatted string can be used instead; in this case PROSSTT will parse the topology from the Newick string. It is important to note that just like in the manual definition, the Newick tree describes branches and their connections::

    import newick
    from prosstt.tree import Tree

    newick_string = "(B:100,C:60)A:70;"
    newick_tree = newick.loads(newick_string)
    t = Tree.from_newick(newick_tree, modules=15, genes=1000, density=None)

The ``topology``, ``time``, ``num_branches``, ``branch_points``, and ``root`` parameters can be parsed out of the newick tree. The same considerations about ``density`` as before apply.

Automatic tree generation
-------------------------

For higher order topologies (above 3 leaves in the tree) multiple topological arrangements are possible (also see Supplemental Material of the paper). As was the case in the MERLoT_ benchmark, sometimes users want to generate a tree topology with the number of bifurcations as the only input parameter. For example, for a lineage tree with 4 bifurcations (9 branches) where each branch has the same length, 50 pseudotime units::

    from prosstt.tree import Tree
    import string

    branch_names = list(string.ascii_uppercase[:9])
    time = {branch: 50 for branch in branch_names}
    t = Tree.from_random_topology(branch_points=4, time=time, modules=15, genes=1000)

.. _MERLoT: https://www.biorxiv.org/content/early/2018/02/08/261768