
export SFASim;
export SFASimMomentum;

export CustomFELSFASim;
export SFAAmpSim;
export create_custom_fel_from_svea;

"""Abstract type representing simulation conditions."""
abstract type AbstractSimSettings end

"""Abstract type representing simulation conditions based on energy."""
abstract type AbstractEnergySimSettings end

"""Abstract type representing simulation conditions based on momentum."""
abstract type AbstractMomentumSimSettings end

"""
Structure containing all information needed for the simulation, including the
observation grid and all initial condition axes.

# Elements:
- `W_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron kinetic energy grid.
- `θ_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron angle grid.
- `ωx::Real`: Average FEL photon energy.
- `ω_axis::AbstractVector{<:Real}`: FEL energy grid.
- `Al_axis::AbstractVector{<:Real}`: Streaking laser amplitude grid.
- `ωl::Real`: Streaking laser angular frequency.
- `Ip::Real`: Gas binding energy.
- `ellipticity::Real`: Ellipticity of the streaking laser.
"""
struct SFASim <: AbstractEnergySimSettings
    # grid of observed energy and azimuthal angle
    W_axis::AbstractVector{<:Real}
    θ_axis::AbstractVector{<:Real}
    # FEL mean angular frequency
    ωx::Real
    # FEL grid
    ω_axis::AbstractVector{<:Real}
    # Streaking field
    Al_axis::AbstractVector{<:Real}
    ωl::Real
    # Ionization potential
    Ip::Real
    # Ellipticity
    ellipticity::Real
    # cross section
    #σ
    # tda
    tda
    # whether to ignore the value of tda
    flat_tda::Bool
    # beta
    beta::Real
    # polarization
    polarization
end

photoelectron_energy_axis(p::SFASim)::AbstractVector{<:Real} = p.W_axis
photoelectron_angle_axis(p::SFASim)::AbstractVector{<:Real} = p.θ_axis
streaking_field_axis(p::SFASim)::AbstractVector{<:Real} = p.Al_axis
ellipticity(p::SFASim)::Real = p.ellipticity
binding_energy(p::SFASim)::Real = p.Ip
streaking_field_energy(p::SFASim)::Real = p.ωl

"""
Structure containing all information needed for the simulation, including the
observation grid and all initial condition axes.

# Elements:
- `px_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron x momentum.
- `py_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron y momentum.
- `ωx::Real`: Average FEL photon energy.
- `ω_axis::AbstractVector{<:Real}`: FEL energy grid.
- `Al_axis::AbstractVector{<:Real}`: Streaking laser amplitude grid.
- `ωl::Real`: Streaking laser angular frequency.
- `Ip::Real`: Gas binding energy.
- `ellipticity::Real`: Ellipticity of the streaking laser.
"""
struct SFASimMomentum <: AbstractMomentumSimSettings
    # grid of observed energy and azimuthal angle
    px_axis::AbstractVector{<:Real}
    py_axis::AbstractVector{<:Real}
    pz_axis::AbstractVector{<:Real}
    # FEL mean angular frequency
    ωx::Real
    # FEL grid
    ω_axis::AbstractVector{<:Real}
    # Streaking field
    Al_axis::AbstractVector{<:Real}
    ωl::Real
    # Ionization potential
    Ip::Real
    # Ellipticity
    ellipticity::Real
    # tda
    tda
    # whether to ignore the value of tda
    flat_tda::Bool
    # beta
    beta::Real
    # polarization
    polarization
end

photoelectron_px(p::SFASimMomentum)::AbstractVector{<:Real} = p.px_axis
photoelectron_py(p::SFASimMomentum)::AbstractVector{<:Real} = p.py_axis
photoelectron_pz(p::SFASimMomentum)::AbstractVector{<:Real} = p.pz_axis
streaking_field_axis(p::SFASimMomentum)::AbstractVector{<:Real} = p.Al_axis
ellipticity(p::SFASimMomentum)::Real = p.ellipticity
binding_energy(p::SFASimMomentum)::Real = p.Ip
streaking_field_energy(p::SFASimMomentum)::Real = p.ωl

"""
Structure containing all information needed for the custom E-field simulation, including the
observation grid and all initial condition axes.

# Elements:
- `W_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron kinetic energy grid.
- `θ_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron angle grid.
- `Al_axis::AbstractVector{<:Real}`: Streaking laser amplitude grid.
- `ωl::Real`: Streaking laser angular frequency.
- `Ip::Real`: Gas binding energy.
- `ellipticity::Real`: Ellipticity of the streaking laser.
- `E::Function: FEL electric field. Should receive only a time argument.
"""
struct CustomFELSFASim <:AbstractEnergySimSettings
    # grid of observed energy and azimuthal angle
    W_axis::AbstractVector{<:Real}
    θ_axis::AbstractVector{<:Real}
    # FEL mean angular frequency
    ωx::Real
    # Streaking field
    Al_axis::AbstractVector{<:Real}
    ωl::Real
    # Ionization potential
    Ip::Real
    # Ellipticity
    ellipticity::Real
    # custom electric field and time axis
    E #::Function
    # cross section
    #σ
    # tda
    tda
    # whether to ignore the value of tda
    flat_tda::Bool
    # beta
    beta::Real
    # polarization
    polarization
end

function create_custom_fel_from_svea(t::AbstractVector{<:Real},
                                     A_td::AbstractVector{<:Real},
                                     ϕ_td::AbstractVector{<:Real},
                                     W_axis::AbstractVector{<:Real},
                                     θ_axis::AbstractVector{<:Real},
                                     ωx::Real,
                                     Al_axis::AbstractVector{<:Real},
                                     ωl::Real,
                                     Ip::Real,
                                     ellipticity::Real,
                                     tda,
                                     flat_tda::Bool,
                                     beta::Real,
                                     polarization)::CustomFELSFASim
    ϕ_int = linear_interpolation(t, ϕ_td; extrapolation_bc=Flat())
    A_int = linear_interpolation(t, A_td; extrapolation_bc=Flat())
    ϕ_fel = tx -> ϕ_int(tx)
    A_fel = tx -> A_int(tx)
    E = tx -> A_fel(tx)*sin(ωx*tx + ϕ_fel(tx))
    CustomFELSFASim(W_axis,
                    θ_axis,
                    ωx,
                    Al_axis,
                    ωl,
                    Ip,
                    ellipticity,
                    E,
                    tda,
                    flat_tda,
                    beta,
                    polarization);
end

photoelectron_energy_axis(p::CustomFELSFASim)::AbstractVector{<:Real} = p.W_axis
photoelectron_angle_axis(p::CustomFELSFASim)::AbstractVector{<:Real} = p.θ_axis
streaking_field_axis(p::CustomFELSFASim)::AbstractVector{<:Real} = p.Al_axis
ellipticity(p::CustomFELSFASim)::Real = p.ellipticity
binding_energy(p::CustomFELSFASim)::Real = p.Ip
streaking_field_energy(p::CustomFELSFASim)::Real = p.ωl


"""
Structure containing all information needed for the simulation based on amplitudes and the frequency/phase FEL,
including the observation grid and all initial condition axes.

# Elements:
- `W_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron kinetic energy grid.
- `θ_axis::AbstractArray{<:Real}`: Observed eTOF photo-electron angle grid.
- `ωx::Real`: Average FEL photon energy.
- `ω_axis::AbstractVector{<:Real}`: FEL energy grid.
- `α_axis::AbstractVector{<:Real}`: FEL phase grid.
- `Al_axis::AbstractVector{<:Real}`: Streaking laser amplitude grid.
- `ωl::Real`: Streaking laser angular frequency.
- `Ip::Real`: Gas binding energy.
- `ellipticity::Real`: Ellipticity of the streaking laser.
"""
struct SFAAmpSim <: AbstractEnergySimSettings
    # grid of observed energy and azimuthal angle
    W_axis::AbstractVector{<:Real}
    θ_axis::AbstractVector{<:Real}
    # FEL mean angular frequency
    ωx::Real
    # FEL grid
    ω_axis::AbstractVector{<:Real}
    # Streaking field
    Al_axis::AbstractVector{<:Real}
    ωl::Real
    # Ionization potential
    Ip::Real
    # Ellipticity
    ellipticity::Real
    # cross section
    #σ
    # tda
    tda
    # whether to ignore the value of tda
    flat_tda::Bool
    # beta
    beta::Real
    # polarization
    polarization
end

photoelectron_energy_axis(p::SFAAmpSim)::AbstractVector{<:Real} = p.W_axis
photoelectron_angle_axis(p::SFAAmpSim)::AbstractVector{<:Real} = p.θ_axis
streaking_field_axis(p::SFAAmpSim)::AbstractVector{<:Real} = p.Al_axis
ellipticity(p::SFAAmpSim)::Real = p.ellipticity
binding_energy(p::SFAAmpSim)::Real = p.Ip
streaking_field_energy(p::SFAAmpSim)::Real = p.ωl


