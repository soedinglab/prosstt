Simulate average gene expression along tree
===========================================

Having obtained a ``Tree`` object (we will be calling it ``tree``), the next step is then to simulate average gene expression along each branch. In PROSSTT, gene expression is controlled by a small number of gene expression programs. These can be imagined as instructions for gene expression; they describe change in relative gene expression (from 0% to 100% or even more).

Contribution coefficients
-------------------------

Each expression program can influence multiple genes, and the expression of each gene is the weighted sum of all expression programs. These contributions are controlled by coefficients sampled (by default) from a Gamma distribution::

    coefficients = simulate_coefficients(tree, a=0.05)

If a simpler model is desired, the coefficients can also be sampled from a Beta distribution. In this case, each gene will be controlled by a maximum of two expression programs, with the contribution of all others set to zero::

    coefficients = simulate_coefficients(tree, a=2, b=2)

Expression programs
-------------------

Given the program-gene coefficients, we visit the branches of the tree breadth-first and simulate the expression programs. We can then calculate the relative average gene expression over pseudotime as the matrix product of the expression programs and the coefficient matrix. We adjust it so that the relative expression for each gene starts where it ended in the preceding branch::

    from prosstt import simulation as sim
    from prosstt import sim_utils as sut
    import numpy as np

    bfs = sut.breadth_first_branches(tree)

    for branch in bfs:
        programs[branch] = sim.sim_expr_branch(tree.time[branch], tree.modules)
        rel_means[branch] = np.dot(programs[branch], coefficients)
        rel_means[branch] = sut.adjust_to_parent(rel_means, branch, topology)


Each expression program in each branch is a random walk with a momentum term that takes place in log space. PROSSTT requires that expression programs have a low correlation with each other (default: 0.2). This makes sure that a small number of (non-redundant) programs can encode all the variability needed.

At this point, PROSSTT enforces two quality checks: first, it requires that no expression program has a value above a cutoff (default is 8; with a relative expression value of ``exp(8)``, a gene is upregulated to almost 3000 times its original expression level). Second, it requires that expression programs move differently in branches that are parallel in pseudotime (branches that share a branch point).

More specifically, PROSSTT calculates the Pearson correlation coefficient between the time course of each expression program in both branches. It then requires that a proportion of the expression programs are anticorrelated. This, by default, is set on 0. This check is useful when using a small number of modules, since otherwise it might not be easy to generate enough variability to differentiate the branches in the dimensionality reduction step::

    above_cutoff = (np.max(rel_means[branch]) > rel_exp_cutoff)
    parallels = sut.find_parallel(tree, programs, branch)
    diverges = sut.diverging_parallel(parallels, rel_means, tree.G, tol=inter_branch_tol)

This process results in B matrices with the size :math:`G \times T_b`, where :math:`G` is the number of expression programs and :math:`T_b` the length of branch B in pseudotime units. Each row contains the time evolution (:math:`T_b` steps) of the relative expression of gene :math:`g` for the current branch.

Getting average gene expression along tree
------------------------------------------

Translating the relative average expression to absolute expression is straightforward; the matrices of relative expression just need to be scaled by the base expression of each gene (the "100%" of expression they would have when unperturbed by the differentiation procedure that is being simulated)::

    import sim_utils as sut
    gene_scale = sut.simulate_base_gene_exp(tree, programs)
    Ms = {}
    for branch in tree.branches:
        Ms[branch] = np.exp(programs[branch]) * gene_scale
    tree.add_genes(Ms)

After that, the simulation of the lineage tree is ready and it is time to sample cells from it.