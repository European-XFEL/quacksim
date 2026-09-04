
convert_Et(x) = x

"""
Return (by reference in the first argument) dϕ(t)/dt, where ϕ(t) is a 3-dimensional array.
Its dimensions mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and

This is the exponent in the b-equation (except for the factor of -i/ħ).
The magnitude of the photo-electron momentum is sqrt(2mW), where W is the photo-electron
kinetic energy and m is the electron mass (in atomic units, 1).

# Arguments
- `du::AbstractArray{<:Real}`: The return value.
- `u::AbstractArray{<:Real}`: The current value of ϕ(t).
- `p::SFASim`: The simulation configuration.
- `t::Real`: Current time.
"""
function dϕ!(du::AbstractArray{<:Real}, u::AbstractArray{<:Real}, p::AbstractEnergySimSettings, t::Real)
    W = reshape(photoelectron_energy_axis(p),:,1,1)
    θ = reshape(photoelectron_angle_axis(p),1,:,1)
    Al = reshape(streaking_field_axis(p),1,1,:)
    eps = ellipticity(p)
    Ip = binding_energy(p)
    ωl = streaking_field_energy(p)
    qx = sqrt.(2 .* W).*cos.(θ) .- Al./sqrt(1+eps^2)./c.*cos.(ωl*t)
    qy = sqrt.(2 .* W).*sin.(θ) .- Al.*eps./sqrt(1+eps^2)./c.*sin.(ωl*t)
    du .= (qx.^2 + qy.^2 )./2 .+ Ip
    nothing
end

"""
Return ∫ϕ(t)dt, where ϕ(t) is a 3-dimensional array.
Its dimensions mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and

This is the time-integrated ϕ exponent b-equation (except for the factor of -i/ħ).
The magnitude of the photo-electron momentum is sqrt(2mW), where W is the photo-electron
kinetic energy and m is the electron mass (in atomic units, 1).

# Arguments
- `p::SFASim`: The simulation configuration.
- `t::Real`: Current time.
"""
function analyticϕ(p::AbstractEnergySimSettings, t::Real)::AbstractArray{<:Real}
    W = reshape(photoelectron_energy_axis(p),:,1,1)
    θ = reshape(photoelectron_angle_axis(p),1,:,1)
    Al = reshape(streaking_field_axis(p),1,1,:)
    eps = ellipticity(p)
    Ip = binding_energy(p)
    ωl = streaking_field_energy(p)
    px = sqrt.(2 .* W).*cos.(θ)
    Ax = Al./sqrt(1+eps^2)./c
    py = sqrt.(2 .* W).*sin.(θ)
    Ay = Al.*eps./sqrt(1+eps^2)./c
    u = ωl*t
    # q = (px, py) - (Ax cos(u), Ay sin(u))
    # We need int q^2/2 + Ip dt
    # int q^2/2 + Ip dt = 1/2 int p^2 dt + 1/2 int A^2 dt - 2 1/2 int p A dt + Ip t
    # int p^2 dt ==
    p2 = (px.^2 + py.^2).*t
    # int A^2 dt = int Ax^2 cos^2(u) dt + int Ay^2 sin^2(u) dt
    # int f(u) dt = 1/ωl int f(u) du
    # int cos^2(u) du = 1/2 u + 1/4 sin 2u
    # int sin^2(u) du = 1/2 u - 1/4 sin 2u
    A2 = (Ax.^2).*(0.5.*u .+ 0.25.*sin.(2*u))./ωl .+ (Ay.^2).*(0.5.*u .- 0.25.*sin.(2*u))./ωl
    # int p A dt = int Ax px cos(u) dt + int Ay py sin(u) dt
    # = Ax px (sin(u)/ωl) + Ay py (-cos(u)/ωl)
    pA = Ax.*px.*(sin.(u)./ωl) .+ Ay.*py.*(-cos.(u)/ωl)
    return (p2 .+ A2 .- 2*pA)./2 .+ Ip*t
end

"""
Return ∫ϕ(t)dt, where ϕ(t) is a 3-dimensional array.
Its dimensions mean:
1. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
2. observed eTOF angle (bins as in `p.θ_axis`, in radians);
3. streaking vector potential amplitude (bins as in `p.amp_Al`); and

This is the time-integrated ϕ exponent b-equation (except for the factor of -i/ħ).
The magnitude of the photo-electron momentum is sqrt(2mW), where W is the photo-electron
kinetic energy and m is the electron mass (in atomic units, 1).

# Arguments
- `p::SFASim`: The simulation configuration.
- `t::Real`: Current time.
"""
function analyticϕ(p::AbstractMomentumSimSettings, t::Real)::AbstractArray{<:Real}
    px = reshape(photoelectron_px(p),:,1,1,1)
    py = reshape(photoelectron_py(p),1,:,1,1)
    pz = reshape(photoelectron_pz(p),1,1,:,1)
    Al = reshape(streaking_field_axis(p),1,1,1,:)
    eps = ellipticity(p)
    Ip = binding_energy(p)
    ωl = streaking_field_energy(p)
    Ax = Al./sqrt(1+eps^2)./c
    Ay = Al.*eps./sqrt(1+eps^2)./c
    u = ωl*t
    # q = (px, py) - (Ax cos(u), Ay sin(u))
    # We need int q^2/2 + Ip dt
    # int q^2/2 + Ip dt = 1/2 int p^2 dt + 1/2 int A^2 dt - 2 1/2 int p A dt + Ip t
    # int p^2 dt ==
    p2 = (px.^2 .+ py.^2 .+ pz.^2).*t
    # int A^2 dt = int Ax^2 cos^2(u) dt + int Ay^2 sin^2(u) dt
    # int f(u) dt = 1/ωl int f(u) du
    # int cos^2(u) du = 1/2 u + 1/4 sin 2u
    # int sin^2(u) du = 1/2 u - 1/4 sin 2u
    A2 = (Ax.^2).*(0.5.*u .+ 0.25.*sin.(2*u))./ωl .+ (Ay.^2).*(0.5.*u .- 0.25.*sin.(2*u))./ωl
    # int p A dt = int Ax px cos(u) dt + int Ay py sin(u) dt
    # = Ax px (sin(u)/ωl) + Ay py (-cos(u)/ωl)
    pA = Ax.*px.*(sin.(u)./ωl) .+ Ay.*py.*(-cos.(u)/ωl)
    return (p2 .+ A2 .- 2*pA)./2 .+ Ip*t
end


"""
Specialized verson of db!, where the FEL is custom.

Return (by reference in the first argument) du(t)/dt, where u(t) is a 4-dimensional array.

Its dimensions have the following meaning:
1. u[1,...] is ϕ(t), u[2] is Re{b(t)}, u[3] is Im{b(t)};
2. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
3. observed eTOF angle (bins as in `p.θ_axis`, in radians);
4. streaking vector potential amplitude (bins as in `p.amp_Al`); and

This allows for solving differential equations simultaneously for several bins of W, θ, and Al.

The equation for b depends on ϕ through the integral going from the current time to the time of detection.
This is accomplished by calculating first ϕT, which gives ϕ(t=∞) and calculating the integral starting
at the current time as ϕT - ϕ(t), where ϕ(t) is the integral from zero to the current time.

# Arguments
- `du::AbstractArray{<:Real}`: The return value.
- `u::AbstractArray{<:Real}`: The current value of u(t).
- `p::CustomFELSFASim`: The simulation configuration.
- `t::Real`: Current time.
- `ϕT::AbstractArray{<:Real}`: The integrated ϕ from 0 to infinity. Its dimensions correspond to the photo-electron kinetic energy, the photo-electron angle and the streaking laser amplitude, resectively. It does not change during integration.
"""
function db!(du::AbstractArray{<:Real}, u::AbstractArray{<:Real}, p::CustomFELSFASim, t::Real, ϕT::AbstractArray{<:Real})
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
    NA = length(p.Al_axis)

    # create views for each output variable selecting
    # one dimension
    @inbounds dϕ_t = view(du, 1, ..)
    @inbounds Re_du = view(du, 2, ..)
    @inbounds Im_du = view(du, 3, ..)

    # get the observed grid
    W = reshape(p.W_axis,:,1,1)
    θ = reshape(p.θ_axis,1,:,1)
    Al = reshape(p.Al_axis,1,1,:)
    eps = p.ellipticity
    # get the gas binding energy
    Ip = p.Ip
    # get the laser photon energy
    ωl = p.ωl

    # the following uses dϕ_t to calculate q^2 first to avoid allocating
    # space for qx and qy

    # calculate dϕ/dt
    # (avoid calling dϕ to avoid extra memory allocations, but this is the same)
    qx = sqrt.(2 .* W).*cos.(θ) .- Al./sqrt(1+eps^2)./c.*cos.(ωl*t)
    qy = sqrt.(2 .* W).*sin.(θ) .- Al.*eps./sqrt(1+eps^2)./c.*sin.(ωl*t)
    q2 = (qx.^2 + qy.^2)
    dϕ_t .= q2./2 .+ Ip
    # calculate the electric dipole for the gas atom
    #dx = qx ./ (q2.^3)
    # use √cross section times cos(qθ) instead
    #dx = sqrt.(p.σ.(q2/2 .+ Ip)./sqrt.(q2)).*(qx)
    #qϕ = atan.(qy, qx)
    #dx = orbitals[p.gas]["d"].(sqrt.(q2), π/2, qϕ)
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
    #Ex_t = Ex(p, t)
    #ω = reshape(p.ω_axis, :, 1)
    #Ex_t = exp.(-(t .- κ).^2 ./ (2*(τ.^2))) .* sin.(ω .* (t .- κ))
    Ex_t = [convert_Et(p.E(t))]
    ewp = reshape(dx, NW, Nθ, NA) .* reshape(Ex_t, 1, 1, 1)

    # axes:
    # 1. observed eTOF kinetic energy (W)
    # 2. observed eTOF angle (θ)
    # 3. streaking vector potential amplitude (amp_Al)

    # get the current ϕ
    @inbounds ϕ = view(u, 1, ..)

    # calculate db/dt in its real and imaginary parts
    Re_du .= ewp.*cos.(-(ϕT .- ϕ))
    Im_du .= ewp.*sin.(-(ϕT .- ϕ))
    nothing
end

"""
Specialized verson of db!, where the FEL is custom and ϕ is calculated analytically.

Return (by reference in the first argument) du(t)/dt, where u(t) is a 4-dimensional array.

Its dimensions have the following meaning:
1. u[1] is Re{b(t)}, u[2] is Im{b(t)};
2. observed eTOF kinetic energy (bins as in `p.W_axis`, in a.u.);
3. observed eTOF angle (bins as in `p.θ_axis`, in radians);
4. streaking vector potential amplitude (bins as in `p.amp_Al`); and

This allows for solving differential equations simultaneously for several bins of W, θ, and Al.

The equation for b depends on ϕ through the integral going from the current time to the time of detection.
This is accomplished by calculating first ϕT, which gives ϕ(t=∞) and calculating the integral starting
at the current time as ϕT - ϕ(t), where ϕ(t) is the integral from zero to the current time.

# Arguments
- `du::AbstractArray{<:Real}`: The return value.
- `u::AbstractArray{<:Real}`: The current value of u(t).
- `p::CustomFELSFASim`: The simulation configuration.
- `t::Real`: Current time.
- `ϕT::AbstractArray{<:Real}`: The integrated ϕ from 0 to infinity. Its dimensions correspond to the photo-electron kinetic energy, the photo-electron angle and the streaking laser amplitude, resectively. It does not change during integration.
"""
function db_analyticϕ!(du::AbstractArray{<:Real}, u::AbstractArray{<:Real}, p::CustomFELSFASim, t::Real, ϕT::AbstractArray{<:Real})
    NW, Nθ = length(p.W_axis), length(p.θ_axis)
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
    #qϕ = atan.(qy, qx)
    #dx = orbitals[p.gas]["d"].(sqrt.(q2), π/2, qϕ)
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
    #Ex_t = Ex(p, t)
    #ω = reshape(p.ω_axis, :, 1)
    #Ex_t = exp.(-(t .- κ).^2 ./ (2*(τ.^2))) .* sin.(ω .* (t .- κ))
    Ex_t = [convert_Et(p.E(t))]
    ewp = reshape(dx, NW, Nθ, NA) .* reshape(Ex_t, 1, 1, 1)

    # axes:
    # 1. observed eTOF kinetic energy (W)
    # 2. observed eTOF angle (θ)
    # 3. streaking vector potential amplitude (amp_Al)

    # get the current ϕ
    ϕ_current = reshape(ϕT - analyticϕ(p, t), NW, Nθ, NA, 1)

    # calculate db/dt in its real and imaginary parts
    Re_du .= ewp.*cos.(-ϕ_current)
    Im_du .= ewp.*sin.(-ϕ_current)
    nothing
end

