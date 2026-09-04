#!/usr/bin/env julia

using ArgParse
using HDF5
using LaTeXStrings
using Plots
using Format

using Random, Distributions
using EllipsisNotation

using LinearInterpolations
using Adapt
using CUDA

using QuackSim

function main()
    s = ArgParseSettings(description = "Simulate observations in angular streaking for a particular electric FEL field previously simulated with the SASE imitator.",
                         commands_are_required = false
                         )
    @add_arg_table s begin
        "--input"
            arg_type = String
            required = true
            action = :store_arg
            help = "Name of the input file, output of sase_imitator.py."
        "--output"
            arg_type = String
            required = true
            action = :store_arg
            help = "Name of the output HDF5 file name."
        "--ellipticity"
            arg_type = Float64
            required = false
            action = :store_arg
            default = 1.0
            help = "Optical laser ellipticity."
        "--min-W"
            arg_type = Float64
            required = true
            action = :store_arg
            help = "Minimum photo-electron kinetic energy, in eV."
        "--max-W"
            arg_type = Float64
            required = true
            action = :store_arg
            help = "Maximum photo-electron kinetic energy, in eV."
        "--nbins-W"
            arg_type = Int32
            required = true
            action = :store_arg
            help = "Number of bins for the photo-electron kinetic energy."
        "--laser-wavelength"
            arg_type = Float64
            required = true
            metavar = "WAVELENGTH"
            help = "Laser wavelength, in nm."
        "--use-svae"
            action = :store_true
            help = "Use SVAE input data?"
        "--gas"
            arg_type = String
            required = false
            action = :store_arg
            default = "Ne1s"
            metavar = "GAS"
            help = "Gas and orbital."
        "--Up"
            arg_type = Float64
            default = 0.4
            metavar = "ENERGY"
            help = "Ponderomotive potential, in eV."
        "--random-Up"
            action = :store_true
            help = "If set, produce random Up values from zero up to the value of --Up."
        "--beta"
            arg_type = Float64
            default = 2.0
            metavar = "BETA"
            help = "Beta for the angular distribution."
        "--use-gpu"
            action = :store_true
            help = "Use GPU?"
        "--max-events"
            arg_type = Int32
            default = -1
            action = :store_arg
            help = "Maximum number of events to process."
    end
    args = parse_args(ARGS, s)

    use_svae = false
    if args["use-svae"] != nothing
        if args["use-svae"]
            use_svae = true
        end
    end

    use_gpu = false
    if args["use-gpu"] != nothing
        if args["use-gpu"]
            use_gpu = true
        end
    end
    println(format("Use GPU? {1}", use_gpu))

    E = nothing
    t = nothing
    td_amp = nothing
    td_phase = nothing
    td_axis = nothing
    # input
    fid = h5open(args["input"], "r")
    # what matters
    if !use_svae
        # E-field
        E = read(fid["E"])
        # time in au
        t = read(fid["t"])
    else
        td_amp = read(fid["time_domain_amplitude"])
        td_phase = read(fid["time_domain_phase"])
        td_axis = read(fid["time_domain_axis"])
    end

    # for cross checks
    # spectra
    E_spec = read(fid["E_spec"])
    t_spec = read(fid["t_spec"])
    spec = read(fid["spec"])
    # axes
    E_axis = read(fid["spec_E_axis"])
    t_axis = read(fid["spec_t_axis"])

    E_photon = read(fid["E_photon"])
    width = read(fid["width"])
    envelope_width = read(fid["envelope_width"])
    bw = read(fid["bw"])

    close(fid)

    #W_axis = collect(LinRange(100.0/eV_per_au, 160.0/eV_per_au, 160))
    W_axis = collect(LinRange(args["min-W"]/eV_per_au, args["max-W"]/eV_per_au, args["nbins-W"]))
    θ_axis = collect(range(0, 2π, step=2π/16.0))[1:16]
    θ_axis_deg = rad2deg.(θ_axis)
    NW, Nθ = length(W_axis), length(θ_axis)

    output = args["output"]

    gas = args["gas"]
    Ip = orbitals[gas]["Ip"]/eV_per_au
    tda = orbitals[gas]["tda_int"]
    #σ = adapt(CuArray{Float32}, σ);
    ωx = E_photon/eV_per_au
    #amp_Al = args["photoelectron-spread"]/eV_per_au*c*sqrt(2)
    ω0 = ωx - Ip

    Nsim = 0
    if !use_svae
        Nsim = size(E)[2]
    else
        Nsim = size(td_amp)[2]
    end
    if args["max-events"] > 0
        Nsim = min(Nsim, args["max-events"])
        spec = spec[.., 1:Nsim]
        E_spec = E_spec[.., 1:Nsim]
        t_spec = t_spec[.., 1:Nsim]
        width = width[1:Nsim]
        bw = bw[1:Nsim]
        envelope_width = envelope_width[1:Nsim]
    end
    Random.seed!(123)

    ellipticity = args["ellipticity"]

    #photoelectron_sigma = args["photoelectron-spread"]/eV_per_au
    #amp_Al = (sqrt(2*ω0 + 2*photoelectron_sigma) - sqrt(2*ω0))*c*sqrt(2)
    #Up = (amp_Al/sqrt(2))^2/(4*c^2)
    #Up = args["Up"]/eV_per_au
    #amp_Al = sqrt(4Up)*c*sqrt(2)
    Up = ones((Nsim,)).*args["Up"]./eV_per_au
    if args["random-Up"] != nothing
        Up .= rand(Float64, (Nsim,)).*args["Up"]./eV_per_au
    end
    A_axis = sqrt.(4Up).*c.*sqrt(1+ellipticity^2)


    laser_wavelength = args["laser-wavelength"]
    λl = laser_wavelength/nm_per_au
    ωl = 2π*c/λl
    Tl = 2π/ωl
    tspan=(0.0, Tl)
    beta = args["beta"]
    println(format("beta = {1:.2f}", beta))

    if use_gpu
        W_axis = cu(W_axis)
        θ_axis = cu(θ_axis)
        tda = Dict(k=>adapt(CuArray{Float32}, v) for (k, v) in tda)
    end

    # create the object holding all information for the integration
    result = []
    for i in 1:Nsim
        this_A = [A_axis[i]]
        if use_gpu
            this_A = cu(this_A)
        end
        if !use_svae
            E_fel_i = Interpolate(t, E[:,i]; extrapolate=:replicate)
            E_fel = tx -> E_fel_i(tx)
            if use_gpu
                E_fel = adapt(CuArray{Float32}, tx->E_fel_i(tx))
            end
        else
            ϕ_int = Interpolate(td_axis, td_phase[:,i]; extrapolate=:replicate)
            A_int = Interpolate(td_axis, td_amp[:,i]; extrapolate=:replicate)
            ϕ_fel = tx -> ϕ_int(tx)
            A_fel = tx -> A_int(tx)
            if use_gpu
                ϕ_fel = adapt(CuArray{Float32}, tx -> ϕ_int(tx))
                A_fel = adapt(CuArray{Float32}, tx -> A_int(tx))
            end
            E_fel = tx -> A_fel(tx)*sin(ωx*tx + ϕ_fel(tx))
        end
        p = CustomFELSFASim(W_axis, θ_axis,                   # observation Grid (eV, rad)
                            ωx,                               # FEL mean energy
                            this_A,                           # streaking laser vector potential amplitude (V)
                            ωl,                               # streaking laser angular frequency
                            Ip,                               # binding energy
                            ellipticity,                      # ellipticity
                            E_fel,                            # FEL electric field
                            tda,
                            false,
                            beta,
                            :horizontal
            )
        # returns an array of size: (NW, Nθ, NA) = (NW, Nθ, 1)
        # sum over all pairs of FEL times and energies
        @time B = calculate(p,
                            tspan                 # time span
                            )
        # normalize
        B_cpu = Array(B)./sum(B);
        # take the only vector potential
        B_cpu = B_cpu[:,:,1]
        # store with an extra index, so we can concatenate over it later
        push!(result, B_cpu[[CartesianIndex()], ..])
    end
    result = cat(result..., dims=1);

    if use_gpu
        W_axis = Array(W_axis)
        θ_axis = Array(θ_axis)
    end

    h5open(output, "w") do fid
        # observation
        fid["obs", chunk=(1,size(result)[2:end]...), deflate=5] = result
        # photo-electron kinetic energy
        fid["W_axis"] = W_axis*eV_per_au;
        # photo-electron angle
        fid["theta_axis"] = θ_axis

        fid["true_Up"] = Up*eV_per_au
        fid["amp_Al"] = A_axis
        fid["Ip"] = Ip*eV_per_au
        fid["Tl"] = Tl*fs_per_au
        fid["wl"] = ωl/fs_per_au

        fid["ellipticity"] = ellipticity

        # spectra
        fid["E_spec", chunk=(1, size(E_spec)[2:end]...), deflate=5] = E_spec
        fid["t_spec", chunk=(1, size(t_spec)[2:end]...), deflate=5] = t_spec
        fid["spec", chunk=(1, size(spec)[2:end]...), deflate=5] = spec
        fid["spec_E_axis", chunk=(1, size(E_axis)[2:end]...), deflate=5] = E_axis
        fid["spec_t_axis", chunk=(1, size(t_axis)[2:end]...), deflate=5] = t_axis

        fid["width"] = width
        fid["envelope_width"] = envelope_width
        fid["bw"] = bw

        fid["gas"] = gas
        fid["beta"] = beta
    end
end

main()
