# QuackSim

Julia module used to simulate angular streaking experiments.
This is necessary for building the basis used in the `quack` repository to reconstruct FEL times and energies. In `quack`, one must provide
the H5 file produced here as input.

See the example notebooks in the `notebooks` directory to see how one could:
* use `Producing_simulation_basis.ipynb` to produce an H5 file containing basis for the reconstruction of the angular streaking observations (needed for QUACK);
* use `Simulation.ipynb` to simulate a FEL-gas interaction producing a photo-electron, and the interaction with the streaking laser field.

## Installation

This package is written in Julia, but it can be used from a Julia notebook, or a Python notebook.

### Usage in Python

One can simply do the following in this directory to install the package. It will use `JuliaCall` to interface between Python and Julia.

```
pip install .
```

See how to produce a basis set below.

### Producing a basis set for QUACK

This can be done by following the steps in the notebook `Producing_simulation_basis.ipynb`. Please follow the instructions and comments in the notebook.
The output file should be provided to Quack.

### Usage in Julia

If in Maxwell, the Julia installation may be loaded with `module load exfel julia` and this package is already available.

If installing Julia locally, it is recommended to use (juliaup)[https://github.com/JuliaLang/juliaup] to manage Julia versions.
Julia version greater than 1.10 is recommended.
Otherwise, one may install this package as follows:

```
julia
]
add https://github.com/European-XFEL/quacksim.git
```

The following modules in Julia are necessary to run the notebooks.

```
$ julia
julia> ]
pkg> add Statistics EllipsisNotation Format CUDA KernelAbstractions DifferentialEquations HDF5 Tullio
pkg> add LaTeXStrings Plots Random Distributions
```

The second set of packages shown are only needed for the notebooks themselves, but not for the `QuackSim` module.

### Command line experiment simulation

This can be done with the notebook `Simulation_Julia.ipynb` or using the command line, as follows. For the notebook,
open the notebook and follow the given instructions.

Given a time-parametrization of the electric field, one can simulate the effect of angular streaking using the script `simulate_experiment.jl`.

A basis set simulation may be done as follows for the Coherent QUACK algorithm:

```
cd quacksim/scripts
julia simulate_experiment.jl --input input.h5 \
                             --output output.h5 \
                             --laser-wavelength 1030.0 \
                             --gas Ne1s \
                             --Up 0.4
```

The parameters have the following meaning:
  * `--input`: an input H5 file produced with the SASE imitator.
  * `--output`: the output file with observations after angular streaking.
  * `--laser-wavelength` indicates the laser wavelength in nm.
  * `--gas` indicates the gas and orbital type. It must match one of the entries in QuackSim internal database. Often Ne 1s is used, which corresponds to `Ne1s`.
  * `--Up` indicates the ponderomotive potential, in eV.


