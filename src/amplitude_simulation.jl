
"""
Similar to `calculate` above, but perform calculation for the given FEL electric field.

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
- `p::CustomFELSFASim`: Simulation configuration, which establishes the Grid setup.
"""
function calculate(p::CustomFELSFASim, tspan::Tuple{<:Real, <:Real}; integration::Symbol=:auto, analytic_action::Bool=true)::AbstractArray{<:Real}
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
    NA = length(p.Al_axis)
    Tl = 2π/p.ωl
    # println("Settings:")
    # println("  - FEL mean energy, hbar ωx/eV = $(ωx.*eV_per_au)")
    # println("  - Time span/fs = $(tspan.*fs_per_au)")
    # println("  - Time step, dt/fs = $(dt.*fs_per_au)")
    # println("  - Photoelectron kinetic energy, W/eV = $(p.W_axis.*eV_per_au)")
    # println("  - Ponderomotive potential, Up/eV = $(p.Al_axis.^2/(4*c.^2).*eV_per_au)")
    # println("  - FEL energy, hbar ω/eV = $(p.ω_axis.*eV_per_au)")
    # println("  - FEL time width, τ/fs = $(p.τ.*fs_per_au)")
    # println("  - Observation's bins: NW, Nθ = $(NW), $(Nθ)")
    # println("  - Spectrogram's bins: Nω = $(Nω)")
    # println("  - Ponderomotive potential bins: NA = $(NA)")

    if integration == :manual
        ϕT = similar(p.W_axis, NW, Nθ, NA)
        fill!(ϕT, 0.0)

        Ts, Te = 0.0, Tl
        all_t_axis = collect(range(Ts, Te, step=dt))
        t_axis = collect(range(tspan[1], tspan[2], step=dt))

        log("Calculating ϕ(t=∞)...")
        rk4!(dϕ!, ϕT, p, all_t_axis)
        func! = (du, u, p, t) -> db!(du, u, p, t, ϕT)

        uT = similar(p.W_axis, 3, NW, Nθ, NA)
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
        elseif integration == :Rodas5P
            alg = Rodas5P()
        elseif integration == :Rodas4P
            alg = Rodas4P()
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
        func = analytic_action ? ((du, u, p, t) -> db_analyticϕ!(du, u, p, t, ϕT)) : ((du, u, p, t) -> db!(du, u, p, t, ϕT))
        uT = analytic_action ? similar(p.W_axis, 2, NW, Nθ, NA) : similar(p.W_axis, 3, NW, Nθ, NA)
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

# from here on: for the amplitude based calculation

"""
Return (by reference in the first argument) du(t)/dt, where u(t) is a 6-dimensional array.

Its dimensions have the following meaning:
1. u[1] is Re{b(t)}, u[2] is Im{b(t)};
2. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
3. observed eTOF angle (bins as in `p.θ_axis`, in radians);
4. streaking vector potential amplitude (bins as in `p.amp_Al`); and
5. FEL energy (bins as in `p.ω_axis`, in a.u.).

This allows for solving differential equations simultaneously for several bins of W, θ, Al, and ω.

The equation for b depends on ϕ through the integral going from the current time to the time of detection.
This is accomplished by calculating first ϕT, which gives ϕ(t=∞) and calculating the integral starting
at the current time as ϕT - ϕ(t), where ϕ(t) is the integral from zero to the current time.

# Arguments
- `du::AbstractArray{<:Real}`: The return value.
- `u::AbstractArray{<:Real}`: The current value of u(t).
- `p::SFAAmpSim`: The simulation configuration.
- `t::Real`: Current time.
- `ϕT::AbstractArray{<:Real}`: The integrated ϕ from 0 to infinity. Its dimensions correspond to the photo-electron kinetic energy, the photo-electron angle and the streaking laser amplitude, resectively. It does not change during integration.
"""
function dbamp!(du::AbstractArray{<:Real}, u::AbstractArray{<:Real}, p::SFAAmpSim, t::Real, ϕT::AbstractArray{<:Real}, ω::AbstractArray{<:Real})
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
    Nω = size(ω)[1] #length(p.ω_axis)
    NA = length(p.Al_axis)

    # create views for each output variable selecting
    # one dimension
    @inbounds Re_du = view(du, 1, ..)
    @inbounds Im_du = view(du, 2, ..)

    # get the observed grid
    W = reshape(p.W_axis,:,1,1)
    θ = reshape(p.θ_axis,1,:,1)
    Al = reshape(p.Al_axis,1,1,:)
    eps = p.ellipticity
    # get the gas binding energy
    Ip = p.Ip
    # get the laser photon energy
    ωl = p.ωl

    # calculate d
    qx = sqrt.(2 .* W).*cos.(θ) .- Al./sqrt(1+eps^2)./c.*cos.(ωl*t)
    qy = sqrt.(2 .* W).*sin.(θ) .- Al.*eps./sqrt(1+eps^2)./c.*sin.(ωl*t)
    q2 = (qx.^2 + qy.^2)
    # calculate the electric dipole for the gas atom
    #dx = qx ./ (q2.^3)
    # use √cross section times cos(qθ) instead
    #dx = sqrt.(p.σ.(q2/2 .+ Ip)./sqrt.(q2)).*(qx)
    qθ = π/2 # always in the z axis
    qϕ = atan.(qy, qx)
    dx = similar(q2)
    if !p.flat_tda
        dx .= sum(tda.(q2/2 .+ Ip) for tda ∈ values(p.tda))
    else
        dx .= 1.0
    end
    if p.polarization == :vertical
        qϕ = π/2 .- qϕ
    end
    dx .*= sqrt.(1/2*(1.0 .+ p.beta*(3*cos.(qϕ).^2 .- 1.0)./2.0))

    # b = int_0^t dt' E . d exp(-i int_t'^t dt'' H)
    # b(inf) = int_0^inf dt' E.d exp(-i int_t'^inf dt'' H)
    # db/dt = E.d exp(-i (int_0^inf dt'' H - int_0^t H dt''))
    # dphi/dt = H
    # db/dt = E.d exp(-i phi(inf) + i phi(t))
    # calculate the FEL electric wave packet as E . d
    v = reshape(-ω.*t, 1, 1, 1, Nω)
    ewp = reshape(dx, NW, Nθ, NA, 1)

    # axes:
    # 1. observed eTOF kinetic energy (W)
    # 2. observed eTOF angle (θ)
    # 3. streaking vector potential amplitude (amp_Al)
    # 4. FEL frequency (ω)

    # get the current ϕ
    ϕ_current = reshape(ϕT - analyticϕ(p, t), NW, Nθ, NA, 1)

    # calculate db/dt in its real and imaginary parts
    Re_du .= ewp.*cos.(v .- ϕ_current)
    Im_du .= ewp.*sin.(v .- ϕ_current)
    nothing
end

"""
Calculate the amplitude of a photon observation.

First calculate the integral of the action S from t=0 to t=∞ for ϕ and save it in ϕT.
Then calculate the integral of the scattering amplitude b up until detection.
Return b.

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
- `p::SFAAmpSim`: Simulation configuration, which establishes the Grid setup.
"""
function calculate_amplitude(p::SFAAmpSim, ω::AbstractArray{<:Real}; integration::Symbol=:auto)::AbstractArray{<:Complex}
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
    Nω = size(ω)[1]
    NA = length(p.Al_axis)
    Tl = 2π/p.ωl

    # trigger message every time we go over 1 percent of the laser period
    T = Tl
    condition_printout(u, t, integrator) = (Int(round(integrator.t/T*10)) - Int(round(integrator.tprev/T*10)) >= 1)
    cb = DiscreteCallback(condition_printout, callback_printout!; save_positions=(false,false));

    ωmax = maximum(ω)
    dt = 2π/ωmax/100
    Ts, Te = 0, Tl
    alg = nothing
    kwargs = Dict{Symbol, Any}(:dt => dt,
                               :abstol => 1e-7,
                               :reltol => 1e-5,
                               :callback => cb
                               )
    if integration == :Euler
        alg = Euler()
    elseif integration == :Tsit5
        alg = Tsit5()
    elseif integration == :Rodas5P
        alg = Rodas5P()
    elseif integration == :Rodas4P
        alg = Rodas4P()
    elseif integration == :KenCarp4
        alg = KenCarp4()
    end

    ϕT = similar(p.W_axis, NW, Nθ, NA)
    log("Using analytical result for ϕ(t=Tl)...")
    ϕT .= analyticϕ(p, Te)

    log("Calculating b(t=Tl)...")
    func = ((du, u, p, t) -> dbamp!(du, u, p, t, ϕT, ω))
    uT = similar(p.W_axis, 2, NW, Nθ, NA, Nω)
    fill!(uT, 0.0)

    prob = ODEProblem(func,                         # integrand
                      uT,                           # initial condition
                      (Ts, Te),                     # from which time to which time?
                      p
     )
    # solve ODE
    uT = solve(prob, alg; save_everystep=false, saveat=Te, kwargs...).u[end]
    ϕT = nothing

    # axes:
    # 1. observed eTOF kinetic energy (W)
    # 2. observed eTOF angle (θ)
    # 3. streaking vector potential amplitude (amp_Al)
    # 4. FEL frequency (ω)

    return uT[1,..] + uT[2,..].*1im;
end

"""

Return the amplitude b of observing the corresponding photo-electron in the eTOFs,
given an initial FEL produced with energy ħω.
and the streaking setup in `p`.

Calculations done in `calculate_amplitude`. This is only a wrapper.

The axes of the returned `AbstractArray{<:Real}` mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and
4. FEL energy (bins as in `p.ω_axis`, in a.u.).

Atomic units can be converted to eV, by multiplying values by `eV_per_au`, converted to fs,
by multiplying by `fs_per_au`, and converted to nm, by multiplying by `nm_per_au`.

# Arguments
-  fel_energy: Set FEL energy, in eV.
- `W_axis`: Array of photo-electron kinetic energies (in eV) observed in the eTOF.
- `theta_axis`: Array of photo-electron detection angles (in rd) observed in the eTOF.
- `Ip`: Gas binding energy.
- `laser_wavelength`: Streaking laser wavelength (in nm).
- `maxUp`: Maximum ponderomotive potential, in eV.
- `omega_spread`: FEL energy axis is built as the center energy +/- this value (in eV).
- `output_filename`: Name of the H5 file in which to store the simulations.
- `ellipticity`: Ellipticity of the streaking laser (1.0 means circular polarization).
- `nbins_omega`: Number of bins in the FEL energy axis.
- `nbins_Up`: Number of bins in the ponderomotive potential.
- `n_parallel`: Number of asynchronous workers to use.
- `use_gpu`: If true, use GPU processing.
"""
function simulate_amplitude(;
                  fel_energy::Real,
                  W_axis::AbstractVector{<:Real},
                  theta_axis::AbstractVector{<:Real},
                  laser_wavelength::Real,
                  maxUp::Real,
                  output_filename::AbstractString,
                  gas::AbstractString="Ne1s",
                  ellipticity::Real=1.0,
                  nbins_time::Integer=64,
                  nbins_Up::Integer=40,
                  n_parallel::Integer=1,
                  use_gpu::Bool=true,
                  oversampling::Int64=2,
                  integration::Symbol=:auto,
                  beta::Real=2.0,
                  polarization=:horizontal
                  )

    #eval(macroexpand(Distributed, quote @everywhere using .QuackSim end))
    θ_axis = theta_axis
    # experimental inputs
    W_axis = W_axis./eV_per_au

    # PES Ne 1s binding energy
    Ip = orbitals[gas]["Ip"]/eV_per_au

    # Laser (4.75 um)
    λl = laser_wavelength/nm_per_au

    # end of experimental parameters
    # start of output configuration (also configurable)

    # laser settings
    Tl = λl/c
    ωl = 2π/Tl

    # FEL set photon energy
    ωx = fel_energy/eV_per_au

    # FEL photon energy axis
    # ifft of energy axis is the time axis and maximum time is the energy sampling rate/2
    # E(t) = F^-1 [E(ω)]
    # t axis: (0, Tl) with nbins_time
    # ω axis: (-0.5ωs, 0.5ωs) for ωs = 1/δt
    # maybe divide by 2π, since time is in units of 2π, as there is no 2π in the FT of E(t)
    δt = Tl/nbins_time
    ωs = 2π/δt
    ω_axis = ωx .+ fftfreq(oversampling*nbins_time, ωs)
    time_axis = collect(LinRange(0.0, oversampling*Tl, oversampling*nbins_time))

    # amplitude of the A-field for streaking: let's bump electrons
    # A = -c int_-inf^t E(t') dt'
    #amp_Al = photoelectron_sigma*c/sqrt(ωx - Ip)*sqrt(2)
    # this would be correct:
    #amp_Al = (sqrt(2*ω0 + 2*photoelectron_sigma) - sqrt(2*ω0))*c*sqrt(2)

    #Al_axis = collect(LinRange(0.0, 1.5*amp_Al, nbins_Up))

    # Up axis
    Up_axis = collect(LinRange(0.0, maxUp, nbins_Up));
    # in terms of the vector potential amplitude
    # Up = Al_axis.^2/(4*c.^2).*eV_per_au;
    Al_axis = sqrt.(4Up_axis./eV_per_au).*c.*sqrt(1+ellipticity^2)

    # end of output configuration

    #kick_axis = sqrt.(4*(ωx - Ip).*eV_per_au.*Up_axis);
    # this would be correct:
    kick_axis = sqrt.(4Up_axis./eV_per_au).*c.*eV_per_au;

    NW, Nθ = length(W_axis), length(θ_axis)
    Nω = length(ω_axis)
    NA = length(Al_axis)

    # create simulation configuration
    itp = orbitals[gas]["tda_int"]
    if use_gpu
        cuitp = Dict(k=>adapt(CuArray{Float32}, v) for (k, v) in itp)
        p = SFAAmpSim(cu(W_axis),
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
        p = SFAAmpSim(W_axis,
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

    chunk_size = Int(Nω)
    k_chunks = collect(Iterators.partition(1:Nω, chunk_size))
    @sync begin
        B = asyncmap(k -> Array(calculate_amplitude(p, p.ω_axis[k]; integration=integration)), k_chunks;
                     ntasks=n_parallel);
        # concatenate results over energy axis
        B = cat(B...; dims=4);
    end
    #B = calculate_amplitude(p, p.ω_axis; integration=integration);
    #B = Array(B)

    # move A axis last
    B = permutedims(B, (1, 2, 4, 3));

    # be sure it is the correct size
    if size(B) != (NW, Nθ, Nω, NA)
        println("Unexpected simulation result size. Size: $(size(B)), expected $((NW, Nθ, Nω, NA)).");
    end

    # permute before saving in HDF5, because Julia and Python have different conventions
    p_B = permutedims(B, (4, 3, 2, 1))

    # save it
    println(format("Saving result in the output file '{1}'", output_filename))
    h5open(output_filename, "w") do fid
        fid["b", chunk=(NA, Nω, Nθ, 1)] = p_B;

        # photo-electron kinetic energy
        fid["W_axis", chunk=size(W_axis), deflate=5] = W_axis*eV_per_au;
        # photo-electron angle
        fid["theta_axis", chunk=size(θ_axis), deflate=5] = θ_axis;
        # FEL energy
        fid["omega_axis", chunk=size(ω_axis), deflate=5] = ω_axis*eV_per_au;
        # time axis
        fid["time_axis", chunk=size(time_axis), deflate=5] = time_axis*fs_per_au;
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
        fid["gas"] = gas
        # oversampling
        fid["oversampling"] = oversampling
        fid["beta"] = beta;
    end

    nothing
end


