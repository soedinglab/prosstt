Getting started
===============

Installation
------------

PROSSTT can be installed using the ``pip`` package manager, or any package manager compatible with ``pip``. ::

    git clone https://github.com/soedinglab/prosstt.git
    cd prosstt
    pip install prosstt


A simple example
----------------

This is a short program that would generate a simple single bifurcation. The biology behind this example is a differentiation that leads to a simple cell fate decision. Cells start in the undifferentiated type ``0`` and develop into either type ``1`` or type ``2``::
    
    import numpy as np
    from prosstt import simulation as sim
    from prosstt.tree import Tree

    # topology: three branches, one branching point
    #         -- 1 --
    # -- 0 --|
    #         -- 2 --
    top = [[0,1], [0,2]]

    # make tree using default parameters for branch
    # length and complexity
    t = Tree.from_topology(top)

    # simulate using default parameters and sampling
    # 2 cells for each pseudotime point
    X, pseudotimes, branches = sim.restricted_simulation(t)

Of the return values, ``X`` is an :math:`N \times G` matrix, where :math:`N` is the number of cells and :math:`G` the number of genes. ``pseudotimes`` is a vector of length :math:`N` that holds the pseudotime of each cell, quantifying progress through the differentiation. ``branches`` is a vector of length :math:`N` that contains the branch assignment for each cell.


Use cases
---------

For additional use cases look in <scripts> and <use documentation>.