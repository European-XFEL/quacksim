
module QuackSim

using Statistics
using EllipsisNotation
using Format
using CUDA, KernelAbstractions
using OrdinaryDiffEqDefault, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqTsit5
using Tullio

using HDF5

import Dates
import Interpolations
using Adapt
import DelimitedFiles

using FFTW

using ReadableRegex
using Interpolations

include("constants.jl")
include("types.jl")

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

"""
Initialize the module.
"""
function __init__()
    current_dir = @__DIR__
    fname = joinpath(current_dir, "..", "data", "samples.h5")
    println("Reading sample information from $(fname)")
    global orbitals = get_sample_data(fname)
end

end
