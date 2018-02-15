### PROSSTT

PROSSTT (PRObabilistic Simulations of ScRNA-seq Tree-like Topologies) is a package with code for the simulation of scRNAseq data for dynamic processes such as cell differentiation. PROSSTT is open source GPL-licensed software implemented in Python.

Single-cell RNAseq is revolutionizing cellular biology, and many algorithms are developed for the analysis of scRNAseq data. PROSSTT provides an easy way to test the performance of trajectory inference methods on realistic data with a known "gold standard". The algorithm can produce datasets with arbitrary topologies while simulating an arbitrary number of sampled cells and genes.

### Installation

PROSSTT can be installed using the `pip` package manager or any `pip`-compatible package manager:

	git clone https://github.com/soedinglab/prosstt.git
	cd prosstt
	pip install .

### Dependencies

PROSSTT was developed and tested in [Python 3.5](https://www.python.org/downloads/release/python-350/) and [3.6](https://www.python.org/downloads/release/python-360/). While older Python 3 versions should work, there is no guarantee that they will. PROSSTT requires:

* [`numpy`](www.numpy.org), for data structures
* [`scipy`](https://www.scipy.org/), for probabilistic distributions and special functions
* [`pandas`](https://pandas.pydata.org/), for I/O

We also recommend the following libraries:

* [`matplotlib`](https://matplotlib.org/), for plotting
* [`jupyter`](http://jupyter.readthedocs.io/en/latest/index.html) notebooks, for demonstration and development purposes
* [`scanpy`](https://github.com/theislab/scanpy), for the visualization of simulations via diffusion maps. This requires [anndata](https://github.com/theislab/anndata) and Python 3.6 to work.

### How to use

We provide jupyter notebooks with a [baseline example](https://github.com/soedinglab/prosstt/blob/master/examples/minimal_example.ipynb), a more involved example that explains the choice of [variance parameters](https://github.com/soedinglab/prosstt/blob/master/examples/variance_sim.ipynb), and a notebook that showcases the different [sampling strategies](https://github.com/soedinglab/prosstt/blob/master/examples/density_sampling.ipynb).

Alternatively, we include a python script that can be run on the command line (examples/generate_simN.py) to produce simulations like the ones used in the [MERLoT](https://www.biorxiv.org/content/early/2018/02/08/261768) paper.