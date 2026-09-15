module QuackSim

using Adapt: adapt
using CUDA: CUDA, CuArray, cu
using Dates: Dates
using EllipsisNotation: (..)
using FFTW: fftfreq
using Format: format
using HDF5: h5open
using Interpolations: Interpolations, Flat, linear_interpolation
using OrdinaryDiffEqDefault: OrdinaryDiffEqDefault
using OrdinaryDiffEqLowOrderRK: Euler
using OrdinaryDiffEqRosenbrock: Rodas4P, Rodas5P
using OrdinaryDiffEqSDIRK: KenCarp4
using OrdinaryDiffEqTsit5: DiscreteCallback, ODEProblem, Tsit5, solve
using Tullio: @tullio

include("constants.jl")
include("types.jl")

include("xatom_data.jl")
include("utils.jl")

include("amplitude_derivatives.jl")
include("semiclassical_derivatives.jl")

include("semiclassical_simulation.jl")
include("amplitude_simulation.jl")

export simulate;
export simulate_one;
export calculate;
export orbitals;

export simulate_amplitude;

end
