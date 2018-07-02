Simulate average gene expression along tree
===========================================

Having obtained a ``Tree`` object (we will be calling it ``tree``), the next step is then to simulate average gene expression along each branch. In PROSSTT, gene expression is controlled by a small number of gene expression programs. These can be imagined as instructions for gene expression; they describe change in relative gene expression (from 0% to 100% or even more).

Contribution coefficients
-------------------------

Each expression program can influence multiple genes, and the expression of each gene is the weighted sum of all expression programs. These contributions are controlled by coefficients sampled (by default) from a Gamma distribution:

```
coefficients = simulate_coefficients(tree, a=0.05)
```

If a simpler model is desired, the coefficients can also be sampled from a Beta distribution. In this case, each gene will be controlled by a maximum of two expression programs, with the contribution of all others set to zero.

```
coefficients = simulate_coefficients(tree, a=2, b=2)
```

