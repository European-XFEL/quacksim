# QuackSim Documentation

QuackSim is a Julia module, which produces simulations of angular streaking experiments.
FEL characteristics can be reconstructed from angular streaking experiments using [quack](https://git.xfel.eu/machineLearning/quack),
but a critical element in there is the production of a basis, which is performed in QuackSim.

This page documents several API functions in QuackSim below. A good starting for usage is to peruse the notebooks in the `notebooks` directory.
The following provides a quick description of the notebook's contents:
* `notebooks/Example.ipynb` shows how one can obtain a basis simulation to be used in `quack`. It illustrates how this is done assuming an example configuration from a simulated dataset, but all configurations may be adapted as needed.
* `notebooks/Simulation Example.ipynb` shows how one may produce a simulated dataset, in which an known FEL condition assumed and the effect of the angular streaking is simulated. It may be a nice starting point to produce initial example datasets for tests. In it, the width of each FEL spike is assumed to be the same.
* `notebooks/Simulation Example with Width Scan.ipynb` shows another set of steps to produce simulated datasets, in which a single FEL spike is produced with varying widths.
* `notebooks/Visualization of the Propagation Effect.ipynb` illustrates the effect of different streaking fields visually with simple plots.

## Installation

This package may be installed as follows:

```
$ julia
julia> ]
pkg> add https://git.xfel.eu/machineLearning/quack-sim.git
```

Alternatively, the package may be downloaded and one may simply do the following from a Jupyter notebook, adapting the path to `QuackSim.jl` as appropriate:
```
include("../src/QuackSim.jl")
using .QuackSim
```

In this case, the necessary dependencies may need to be installed manually. The following minimum packages are needed:
```
$ julia
julia> ]
pkg> add Statistics EllipsisNotation Format CUDA KernelAbstractions DifferentialEquations HDF5 Tullio
pkg> add LaTeXStrings Plots Random Distributions
```

## Contents

```@contents
```

## Functions and structures

```@docs
SFASim
simulate
simulate_one
calculate
```

## Constants

```@docs
eV_per_au
c
fs_per_au
nm_per_au
```

## Index

```@index
```
