### PROSSTT

PROSSTT (PRObabilistic Simulations of ScRNA-seq Tree-like Topologies) is a package with code for the simulation of scRNAseq data for dynamic processes such as cell differentiation. PROSSTT is open source GPL-licensed software implemented in Python.

Single-cell RNAseq is revolutionizing cellular biology, and many algorithms are developed for the analysis of scRNAseq data. PROSSTT provides an easy way to test the performance of trajectory inference methods on realistic data with a known "gold standard". The algorithm can produce datasets with arbitrary topologies while simulating an arbitrary number of sampled cells and genes.

### Installation

PROSSTT can be installed using the `pip` package manager or any `pip`-compatible package manager:

	git clone https://github.com/soedinglab/prosstt.git
	cd prosstt
	pip install .

### Dependencies

PROSSTT was developed and tested in Python 3.5 and 3.6. While older Python 3 versions should work, there is no guarantee that they will. PROSSTT itself, the scripts, and notebooks included in the package, require the `numpy`, `scipy`, `pandas` and `matplotlib` libraries. We recommend using [scanpy](https://github.com/theislab/scanpy) to calculate diffusion maps in order to visualize the simulations, which requires [anndata](https://github.com/theislab/anndata) and Python 3.6 to work.

### How to use

We provide jupyter notebooks with a [baseline example](https://github.com/soedinglab/prosstt/blob/master/examples/minimal_example.ipynb), a more involved example that explains the choice of [variance parameters](https://github.com/soedinglab/prosstt/blob/master/examples/variance_sim.ipynb), and a notebook that showcases the different [sampling strategies](https://github.com/soedinglab/prosstt/blob/master/examples/density_sampling.ipynb).

For more information, please refer to the [documentation](http://www.404errorpages.com/).
