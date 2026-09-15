
"""Show message `msg` including the current task number."""
function log(msg::AbstractString)
    #println(format("[{1:3d}] {2}", Threads.threadid(), msg))
    #println(format("[{1}] {2}", objectid(current_task()), msg))
    time = Dates.format(Dates.now(), "e, dd u yyyy HH:MM:SS")
    println(format("[{1}] {2}", time, msg))
end


"""
Build the orbital information and interpolations from `xatom_data`.
"""
function get_sample_data()
    o = Dict{String, Dict{String, Any}}()
    energy = collect(xatom_energy) ./ eV_per_au
    for (orb, d) in xatom_data
        o[orb] = Dict{String, Any}(
            "species" => d.species,
            "orbital" => d.orbital,
            "n" => d.n,
            "l" => d.l,
            "m" => d.m,
            "Ip" => d.Ip,
            "tda_int" => Interpolations.linear_interpolation(energy, d.tda, extrapolation_bc=Interpolations.Line()),
            "cross_section" => Interpolations.linear_interpolation(energy, d.cross_section, extrapolation_bc=Interpolations.Line()),
        )
    end
    return o
end

"""
Dictionary of species known for the gas.
"""
const orbitals = get_sample_data()


"""
Simple Runge-Kutta 4 integrator.

Much slower than implementations in DifferentialEquations.jl, but useful for cross-checking.
"""
function rk4!(f!::Function, u::AbstractArray{<:Real}, p::SFASim, t_axis::Vector{<:Real})
    h = t_axis[2] - t_axis[1]
    k = [similar(u) for i = 1:4]
    for (i, t) in enumerate(t_axis)
        if i % 1000 == 0
            log("At iteration $i/$(length(t_axis))...")
        end
        f!(k[1], u, p, t)
        f!(k[2], u + h*k[1]/2, p, t + h/2)
        f!(k[3], u + h*k[2]/2, p, t + h/2)
        f!(k[4], u + h*k[3], p, t + h)
        u .+= h*(k[1] + 2*k[2] + 2*k[3] + k[4]) ./ 6
    end
    nothing
end

"""
Return |b|^2 after the integration is done.
"""
function calculate_B(u::AbstractArray{<:Real})::AbstractArray{<:Real}
    n = size(u)[1]
    @inbounds Re_du = view(u, n-1, ..)
    @inbounds Im_du = view(u, n, ..)
    return Re_du.^2 .+ Im_du.^2
end

"""
Show status during integration.
"""
function callback_printout!(integrator)
    tspan = integrator.sol.prob.tspan
    T = tspan[2] - tspan[1]
    log(format("Progress: {1:3.0f} % (@ t = {2:7.3f} fs, span: {3:7.3f} fs to {4:7.3f} fs)",
               (integrator.t - tspan[1])/T*100,
               integrator.t*fs_per_au,
               tspan[1]*fs_per_au, tspan[2]*fs_per_au))
    return nothing
end

