
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
    load_cross_section()
    #load_tda()
    load_xatom_tda()
end

end
