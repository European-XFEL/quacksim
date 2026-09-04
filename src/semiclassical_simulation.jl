
"""
First calculate the integral of the action S from t=0 to t=∞ for ϕ and save it in ϕT.
Then calculate the integral of the scattering amplitude b up until detection.
Return the probability B of observing the corresponding photo-electron in the eTOFs.

The calculations are done using the DifferentialEquations.jl package using the Runge-Kutta solver
with adaptive time step Tsit5. This has been observed to produce best results.

The axes of the returned `AbstractArray{<:Real}` mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and
4. FEL energy (bins as in `p.ω_axis`, in a.u.).

Atomic units can be converted to eV, by multiplying values by `eV_per_au`, and converted to fs,
by multiplying by `fs_per_au`.

# Arguments
- `κ::Real`: Time (in a.u.) when the FEL pulse happened, relative to the start of the streaking laser.
- `τ::Real`: Pulse width (in a.u.) of the FEL pulse.
- `p::SFASim`: Simulation configuration, which establishes the Grid setup.
- `tspan::Tuple{<:Real, <:Real}`: Span of time (in a.u.) where to perform the simulation.
"""
function calculate(κ::Real, τ::Real, p::SFASim, tspan::Tuple{<:Real, <:Real}; integration::Symbol=:auto, analytic_action::Bool=true)::AbstractArray{<:Real}
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
    Nω = length(p.ω_axis)
    NA = length(p.Al_axis)
    Tl = 2π/p.ωl
    # println("Settings:")
    # println("  - FEL mean energy, hbar ωx/eV = $(ωx.*eV_per_au)")
    # println("  - Time span/fs = $(tspan.*fs_per_au)")
    # println("  - Time step, dt/fs = $(dt.*fs_per_au)")
    # println("  - Photoelectron kinetic energy, W/eV = $(p.W_axis.*eV_per_au)")
    # println("  - Ponderomotive potential, Up/eV = $(p.Al_axis.^2/(4*c.^2).*eV_per_au)")
    # println("  - Streaking time, κ/fs = $(κ*fs_per_au)")
    # println("  - FEL energy, hbar ω/eV = $(p.ω_axis.*eV_per_au)")
    # println("  - FEL time width, τ/fs = $(p.τ.*fs_per_au)")
    # println("  - Observation's bins: NW, Nθ = $(NW), $(Nθ)")
    # println("  - Spectrogram's bins: Nω = $(Nω)")
    # println("  - Ponderomotive potential bins: NA = $(NA)")

    if integration == :manual
        dt = 2π/p.ωx/20
        Ts, Te = 0.0, Tl
        all_t_axis = collect(range(Ts, Te, step=dt))
        t_axis = collect(range(tspan[1], tspan[2], step=dt))
        ϕT = similar(p.W_axis, NW, Nθ, NA)
        fill!(ϕT, 0.0)

        log("Calculating ϕ(t=∞)...")
        rk4!(dϕ!, ϕT, p, all_t_axis)
        func! = (du, u, p, t) -> db!(du, u, p, t, ϕT, κ, τ)

        uT = similar(p.W_axis, 3, NW, Nθ, NA, Nω)
        fill!(uT, 0.0)

        log("Calculating b(t=∞)...")
        rk4!(func!, uT, p, t_axis)
    else
        # trigger message every time we go over 1 percent of the laser period
        T = tspan[2] - tspan[1]
        condition_printout(u, t, integrator) = (Int(round((integrator.t - tspan[1])/T*10)) - Int(round((integrator.tprev - tspan[1])/T*10)) >= 1)
        cb = DiscreteCallback(condition_printout, callback_printout!; save_positions=(false,false));

        dt = 2π/p.ωx/100
        Ts, Te = 0, Tl
        alg = nothing
        kwargs = Dict{Symbol, Any}(:dt => dt,
                                   :abstol => 1e-6,
                                   :reltol => 1e-4
                                   )
        if integration == :Euler
            alg = Euler()
        elseif integration == :Tsit5
            alg = Tsit5()
        elseif integration == :Rodas4P
            alg = Rodas4P()
        elseif integration == :Rodas5P
            alg = Rodas5P()
        elseif integration == :KenCarp4
            alg = KenCarp4()
        end

        ϕT = similar(p.W_axis, NW, Nθ, NA)
        if analytic_action
            log("Using analytical result for ϕ(t=Tl)...")
            ϕT .= analyticϕ(p, Te)
        else
            log("Calculating ϕ(t=Tl)...")
            fill!(ϕT, 0.0)
            prob = ODEProblem(dϕ!,                         # integrand
                              ϕT,                          # initial condition
                              (Ts, Te),                    # from which time to which time?
                              p
            )
            # solve ODE
            ϕT .= solve(prob, alg; save_everystep=false, saveat=Te, kwargs...).u[end]
        end

        log("Calculating b(t=Tl)...")
        func = analytic_action ? ((du, u, p, t) -> db_analyticϕ!(du, u, p, t, ϕT, κ, τ)) : ((du, u, p, t) -> db!(du, u, p, t, ϕT, κ, τ))
        uT = analytic_action ? similar(p.W_axis, 2, NW, Nθ, NA, Nω) : similar(p.W_axis, 3, NW, Nθ, NA, Nω)
        fill!(uT, 0.0)

        (Ts, Te) = tspan
        prob = ODEProblem(func,                         # integrand
                          uT,                           # initial condition
                          (Ts, Te),                     # from which time to which time?
                          p
         )
        # solve ODE
        uT = solve(prob, alg; save_everystep=false, saveat=Te, kwargs...).u[end]
        ϕT = nothing
    end

    # axes:
    # 1. observed eTOF kinetic energy (W)
    # 2. observed eTOF angle (θ)
    # 3. streaking vector potential amplitude (amp_Al)
    # 4. FEL frequency (ω)

    # square b solution
    log("Calculating observation...")
    # calculate it and place it in the first index
    # this avoids using more memory
    return calculate_B(uT)
end

"""
First calculate the integral of the action S from t=0 to t=∞ for ϕ and save it in ϕT.
Then calculate the integral of the scattering amplitude b up until detection.
Return the probability B of observing the corresponding photo-electron in the eTOFs.

The calculations are done using the DifferentialEquations.jl package using the Runge-Kutta solver
with adaptive time step Tsit5. This has been observed to produce best results.

The axes of the returned `AbstractArray{<:Real}` mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and
4. FEL energy (bins as in `p.ω_axis`, in a.u.).

Atomic units can be converted to eV, by multiplying values by `eV_per_au`, and converted to fs,
by multiplying by `fs_per_au`.

# Arguments
- `κ::Real`: Time (in a.u.) when the FEL pulse happened, relative to the start of the streaking laser.
- `τ::Real`: Pulse width (in a.u.) of the FEL pulse.
- `p::SFASim`: Simulation configuration, which establishes the Grid setup.
- `tspan::Tuple{<:Real, <:Real}`: Span of time (in a.u.) where to perform the simulation.
"""
function calculate(κ::Real, τ::Real, p::SFASimMomentum, tspan::Tuple{<:Real, <:Real}; integration::Symbol=:auto)::AbstractArray{<:Real}
    Nx, Ny, Nz = length(p.px_axis), length(p.py_axis), length(p.pz_axis)
    Nω = length(p.ω_axis)
    NA = length(p.Al_axis)
    Tl = 2π/p.ωl
    # println("Settings:")
    # println("  - FEL mean energy, hbar ωx/eV = $(ωx.*eV_per_au)")
    # println("  - Time span/fs = $(tspan.*fs_per_au)")
    # println("  - Time step, dt/fs = $(dt.*fs_per_au)")
    # println("  - Photoelectron kinetic energy, W/eV = $(p.W_axis.*eV_per_au)")
    # println("  - Ponderomotive potential, Up/eV = $(p.Al_axis.^2/(4*c.^2).*eV_per_au)")
    # println("  - Streaking time, κ/fs = $(κ*fs_per_au)")
    # println("  - FEL energy, hbar ω/eV = $(p.ω_axis.*eV_per_au)")
    # println("  - FEL time width, τ/fs = $(p.τ.*fs_per_au)")
    # println("  - Observation's bins: NW, Nθ = $(NW), $(Nθ)")
    # println("  - Spectrogram's bins: Nω = $(Nω)")
    # println("  - Ponderomotive potential bins: NA = $(NA)")

    if integration == :manual
        dt = 2π/p.ωx/20
        Ts, Te = 0.0, Tl
        all_t_axis = collect(range(Ts, Te, step=dt))
        t_axis = collect(range(tspan[1], tspan[2], step=dt))
        ϕT = similar(p.px_axis, Nx, Ny, Nz, NA)
        fill!(ϕT, 0.0)

        log("Calculating ϕ(t=∞)...")
        rk4!(dϕ!, ϕT, p, all_t_axis)
        func! = (du, u, p, t) -> db!(du, u, p, t, ϕT, κ, τ)

        uT = similar(p.px_axis, 3, Nx, Ny, Nz, NA, Nω)
        fill!(uT, 0.0)

        log("Calculating b(t=∞)...")
        rk4!(func!, uT, p, t_axis)
    else
        # trigger message every time we go over 1 percent of the laser period
        T = tspan[2] - tspan[1]
        condition_printout(u, t, integrator) = (Int(round((integrator.t - tspan[1])/T*10)) - Int(round((integrator.tprev - tspan[1])/T*10)) >= 1)
        cb = DiscreteCallback(condition_printout, callback_printout!; save_positions=(false,false));

        dt = 2π/p.ωx/100
        Ts, Te = 0, Tl
        alg = nothing
        kwargs = Dict{Symbol, Any}(:dt => dt,
                                   :abstol => 1e-6,
                                   :reltol => 1e-4
                                   )
        if integration == :Euler
            alg = Euler()
        elseif integration == :Tsit5
            alg = Tsit5()
        elseif integration == :Rodas4P
            alg = Rodas4P()
        elseif integration == :Rodas5P
            alg = Rodas5P()
        elseif integration == :KenCarp4
            alg = KenCarp4()
        end

        ϕT = similar(p.px_axis, Nx, Ny, Nz, NA)
        log("Using analytical result for ϕ(t=Tl)...")
        ϕT .= analyticϕ(p, Te)

        log("Calculating b(t=Tl)...")
        func = ((du, u, p, t) -> db_analyticϕ!(du, u, p, t, ϕT, κ, τ))
        uT = similar(p.px_axis, 2, Nx, Ny, Nz, NA, Nω)
        fill!(uT, 0.0)

        (Ts, Te) = tspan
        prob = ODEProblem(func,                         # integrand
                          uT,                           # initial condition
                          (Ts, Te),                     # from which time to which time?
                          p
         )
        # solve ODE
        uT = solve(prob, alg; save_everystep=false, saveat=Te, kwargs...).u[end]
        ϕT = nothing
    end

    # axes:
    # 1. observed eTOF px
    # 2. observed eTOF py
    # 3. observed eTOF pz
    # 4. streaking vector potential amplitude (amp_Al)
    # 5. FEL frequency (ω)

    # square b solution
    log("Calculating observation...")
    # calculate it and place it in the first index
    # this avoids using more memory
    return calculate_B(uT)
end

"""
Return the probability B of observing the corresponding photo-electron in the eTOFs,
given an initial FEL produced at time κ (in a.u.) relative to the start of the streaking laser period,
a pulse width τ (in a.u.), and the streaking setup in `p`.

Calculations done in `calculate`. This is only a wrapper.

The axes of the returned `AbstractArray{<:Real}` mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and
4. FEL energy (bins as in `p.ω_axis`, in a.u.).

Atomic units can be converted to eV, by multiplying values by `eV_per_au`, and converted to fs,
by multiplying by `fs_per_au`.

# Arguments
- `κ::Real`: Time (in a.u.) when the FEL pulse happened, relative to the start of the streaking laser.
- `τ::Real`: Pulse width (in a.u.) of the FEL pulse.
- `p::SFASim`: Simulation configuration, which establishes the Grid setup.
- `tspan::Tuple{<:Real, <:Real}`: Span of time (in a.u.) where to perform the simulation.
"""
function simulate_one(κ::Real, τ::Real, p::SFASim; integration::Symbol=:auto)::AbstractArray{<:Real}
    log(format("Calculating for streaking time κ = {1:6.3f} fs, width = {2:6.3f} fs", κ*fs_per_au, τ*fs_per_au))
    tmax = κ + 3τ
    tmin = κ - 3τ
    B = calculate(κ, τ, p, (tmin, tmax); integration=integration);
    B = B./maximum(B, dims=(1,2));
    # p = heatmap(Array(p.W_axis).*eV_per_au,
    #         rad2deg.(Array(p.θ_axis)),
    #         Array(B)[:,:,end,Int(length(p.ω_axis)/2)]',
    #         left_margin=3Plots.mm,
    #         size=(1000, 600),
    #         xlabel="Photo-electron kinetic energy [eV]",
    #         ylabel=L"$\theta$ [deg]")
    # display(p)
    return B;
end

"""

Return the probability B of observing the corresponding photo-electron in the eTOFs,
given an initial FEL produced at time κ (in a.u.) relative to the start of the streaking laser period,
a pulse width τ (in a.u.), and the streaking setup in `p`.

Calculations done in `calculate`. This is only a wrapper.

The axes of the returned `AbstractArray{<:Real}` mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and
4. FEL energy (bins as in `p.ω_axis`, in a.u.).

Atomic units can be converted to eV, by multiplying values by `eV_per_au`, converted to fs,
by multiplying by `fs_per_au`, and converted to nm, by multiplying by `nm_per_au`.

# Arguments
- `W_axis`: Array of photo-electron kinetic energies (in eV) observed in the eTOF.
- `theta_axis`: Array of photo-electron detection angles (in rd) observed in the eTOF.
- `laser_wavelength`: Streaking laser wavelength (in nm).
- `time_resolution`: Array of FEL time widths to simulate (in fs).
- `photoelectron_spread`: Maximum span of observed streaking amplitude (in eV) used to estimate maximum ponderomotive potential.
- `omega_spread`: FEL energy axis is built as the center energy +/- this value (in eV).
- `output_filename`: Name of the H5 file in which to store the simulations.
- `gas`: Gas and orbital setup.
- `ellipticity`: Ellipticity of the streaking laser (1.0 means circular polarization).
- `nbins_omega`: Number of bins in the FEL energy axis.
- `nbins_time`: Number of bins in the FEL time axis.
- `nbins_Up`: Number of bins in the ponderomotive potential.
- `n_parallel`: Number of asynchronous workers to use.
- `use_gpu`: If true, use GPU processing.
"""
function simulate(;
                  W_axis::AbstractVector{<:Real},
                  theta_axis::AbstractVector{<:Real},
                  laser_wavelength::Real,
                  time_resolution::AbstractVector{<:Real},
                  photoelectron_spread::Real,
                  omega_spread::Real,
                  output_filename::AbstractString,
                  gas::AbstractString="Ne1s",
                  ellipticity::Real=1.0,
                  nbins_omega::Integer=32,
                  nbins_time::Integer=64,
                  nbins_Up::Integer=40,
                  n_parallel::Integer=2,
                  use_gpu::Bool=true,
                  integration::Symbol=:auto,
                  beta::Real=2.0,
                  polarization=:horizontal,
                  )

    #eval(macroexpand(Distributed, quote @everywhere using .QuackSim end))
    θ_axis = theta_axis
    # experimental inputs
    W_axis = W_axis./eV_per_au

    # PES Ne 1s binding energy
    Ip = orbitals[gas]["Ip"]/eV_per_au

    # Laser (4.75 um)
    λl = laser_wavelength/nm_per_au

    # rough size of +/- photo-electron shift
    photoelectron_sigma = photoelectron_spread/eV_per_au

    τ_axis = time_resolution./fs_per_au #0.5./(energy_resolution./eV_per_au)
    Nτ = length(τ_axis)

    # end of experimental parameters
    # start of output configuration (also configurable)

    # laser settings
    Tl = λl/c
    ωl = 2π/Tl

    # FEL set photon energy
    ωx = W_axis[Int(length(W_axis)//2)] + Ip

    # FEL photon energy axis
    ω_axis = collect(LinRange(ωx - omega_spread/eV_per_au, ωx + omega_spread/eV_per_au, nbins_omega))

    # streaking field
    κ_axis = collect(LinRange(0.0, Tl, nbins_time))

    # amplitude of the A-field for streaking: let's bump electrons
    # A = -c int_-inf^t E(t') dt'
    #amp_Al = photoelectron_sigma*c/sqrt(ωx - Ip)*sqrt(2)
    # this would be correct:
    amp_Al = (sqrt(2*ω0 + 2*photoelectron_sigma) - sqrt(2*ω0))*c*sqrt(2)

    Al_axis = collect(LinRange(0.0, 1.5*amp_Al, nbins_Up))
    # end of output configuration

    # equivalently
    Up_axis = Al_axis.^2/(4*c.^2).*eV_per_au;
    kick_axis = sqrt.(4Up_axis./eV_per_au).*c.*eV_per_au;

    NW, Nθ = length(W_axis), length(θ_axis)
    Nω, Nκ = length(ω_axis), length(κ_axis)
    NA = length(Al_axis)

    # create simulation configuration
    itp = orbitals[gas]["tda_int"]
    if use_gpu
        cuitp = Dict(k=>adapt(CuArray{Float32}, v) for (k, v) in itp)
        p = SFASim(cu(W_axis),
                   cu(θ_axis),
                   ωx,
                   cu(ω_axis),
                   cu(Al_axis),
                   ωl,
                   Ip,
                   ellipticity,
                   cuitp,
                   false,
                   beta,
                   polarization
            );
        CUDA.synchronize()
    else
        p = SFASim(W_axis,
                   θ_axis,
                   ωx,
                   ω_axis,
                   Al_axis,
                   ωl,
                   Ip,
                   ellipticity,
                   itp,
                   false,
                   beta,
                   polarization
            );
    end

    # calculate for each time point in parallel
    # chunk them in n_parallel to be sure it fits in memory
    #chunk_size = Int(Nκ/n_parallel)
    #κ_chunks = collect(Iterators.partition(collect(1:Nκ), chunk_size))
    #B = Mem.pin(similar(W_axis, NW, Nθ, NA, Nω, Nκ))
    #Threads.@threads for iκ_batch in κ_chunks
    #    for i = iκ_batch
    #        copyto!(B[..,i], simulate_one(κ_axis[i], p))
    #    end
    #end

    @sync begin
        B = asyncmap(κ -> cat([Array(simulate_one(κ, τ, p; integration=integration)) for τ in τ_axis]...; dims=5), κ_axis;
                     ntasks=n_parallel);
        # concatenate results over time axis
        B = cat(B...; dims=6);
    end

    # move A axis last and time width before last
    B = permutedims(B, (1, 2, 4, 6, 5, 3));

    # be sure it is the correct size
    if size(B) != (NW, Nθ, Nω, Nκ, Nτ, NA)
        println("Unexpected simulation result size. Size: $(size(B)), expected $((NW, Nθ, Nω, Nκ, Nτ, NA)).");
    end

    # this is ready to be used in reconstruction
    B_reshape = reshape(permutedims(B, (6, 5, 4, 3, 2, 1)), (NA, Nτ, Nκ*Nω, NW*Nθ));

    println("Calculation finished. Pre-calculating B^T B.");

    # this is pre-calculated for use in the QP method
    # do it in the GPU for speed
    if use_gpu
        Br = cu(B_reshape);
    else
        Br = B_reshape;
    end
    # does not fit in GPU memory
    #@tullio BtB[ia, ib, j, l] := Br[ia, ib, j, k] * Br[ia, ib, l, k]
    BtB = cu(zeros(NA, Nτ, Nκ*Nω, Nκ*Nω))
    for ia = 1:NA
        #@tullio BtB[$ia, ib, j, l] = Br[$ia, ib, j, k] * Br[$ia, ib, l, k]
        for ib = 1:Nτ
            BtB_ = view(BtB, ia, ib, :, :)
            Br_ = view(Br, ia, ib, :, :)
            @tullio BtB_[j, l] = Br_[j, k] * Br_[l, k]
        end
    end
    BtB = Array(BtB);
    Br = nothing;

    # permute before saving in HDF5, because Julia and Python have different conventions
    p_B = permutedims(B, (6, 5, 4, 3, 2, 1))
    #B_reshape = permutedims(B_reshape, (3, 2, 1));
    #BtB_reshape = permutedims(BtB, (3, 2, 1));

    # save it
    println(format("Saving result in the output file '{1}'", output_filename))
    h5open(output_filename, "w") do fid
        fid["B", chunk=(1, 1, Nκ, Nω, Nθ, NW)] = p_B;
        fid["B_reshape", chunk=(1, 1, Nκ*Nω, NW*Nθ), deflate=5] = B_reshape;
        fid["BtB_reshape", chunk=(1, 1, Nκ*Nω, Nκ*Nω), deflate=5] = BtB;

        # photo-electron kinetic energy
        fid["W_axis", chunk=size(W_axis), deflate=5] = W_axis*eV_per_au;
        # photo-electron angle
        fid["theta_axis", chunk=size(θ_axis), deflate=5] = θ_axis;
        # FEL energy
        fid["omega_axis", chunk=size(ω_axis), deflate=5] = ω_axis*eV_per_au;
        # FEL time width
        fid["tau", chunk=size(τ_axis), deflate=5] = τ_axis.*fs_per_au;
        # streaking laser time of interaction
        fid["kappa_axis", chunk=size(κ_axis), deflate=5] = κ_axis*fs_per_au;
        # streaking amplitude
        fid["Al_axis", chunk=size(Al_axis), deflate=5] = Al_axis;
        # equivalently: ponderomotive potential
        fid["Up_axis", chunk=size(Up_axis), deflate=5] = Up_axis;
        # equivalently: kick
        fid["kick_axis", chunk=size(kick_axis), deflate=5] = kick_axis;
        # the angular frequency of the laser
        fid["wl"] = ωl/fs_per_au;
        # the wavelength of the laser
        fid["ll"] = λl*nm_per_au;
        fid["Tl"] = Tl*fs_per_au;
        # Ip
        fid["Ip"] = Ip*eV_per_au;
        # gas
        fid["gas"] = gas;
        fid["beta"] = beta;
    end

    nothing
end

