#!/usr/bin/env julia

using ArgParse
using HDF5
using LaTeXStrings
using Plots
using Format

include("../src/QuackSim.jl")
using .QuackSim

function main()
    s = ArgParseSettings(description = "Simulate observations in angular streaking for a particular detector setting.",
                         commands_are_required = false
                         )
    @add_arg_table s begin
        "--read-axes-from"
            arg_type = String
            required = false
            action = :store_arg
            help = "Name of the file with angle and energy axes as observation."
        "--output"
            arg_type = String
            required = true
            action = :store_arg
            help = "Name of the output HDF5 file name."
        "--laser-wavelength"
            arg_type = Float64
            required = true
            metavar = "WAVELENGTH"
            help = "Laser wavelength, in nm."
        "--max-time-resolution"
            arg_type = Float64
            required = false
            action = :store_arg
            default = -1.0
            metavar = "TIME"
            help = "Maximum time resolution, in fs."
        "--gas"
            arg_type = String
            required = false
            action = :store_arg
            default = "Ne1s"
            metavar = "GAS"
            help = "Gas and orbital."
        "--photoelectron-spread"
            arg_type = Float64
            default = 7.5
            metavar = "ENERGY"
            help = "Spread of observed photoelectron energy, in eV. This is used to estimate the range of ponderomotive potentials to be probed."
        "--omega-spread"
            arg_type = Float64
            default = 15.0
            metavar = "ENERGY"
            help = "Spread around the center W axis in which to estimate the FEL energy, in eV."
        "--nbins-omega"
            arg_type = Int64
            default = 32
            metavar = "BINS"
            help = "Number of bins of FEL energy."
        "--nbins-time"
            arg_type = Int64
            default = 64
            metavar = "BINS"
            help = "Number of bins of FEL time."
        "--nbins-Up"
            arg_type = Int64
            default = 40
            metavar = "BINS"
            help = "Number of bins of ponderomotive potential."
    end
    args = parse_args(ARGS, s)

    # extra input
    if args["read-axes-from"] != nothing
        fid = h5open(args["read-axes-from"], "r")
        # simulated data
        if haskey(fid, "W_axis")
            W_axis = read(fid["W_axis"])
            θ_axis = read(fid["theta_axis"])
        else
            # EuXFEL PES
            W_axis = read(fid["energy_axis"])
            θ_axis = collect(LinRange(0, 2π, 16))
        end
        close(fid)
    else
        fid = h5open("../../quack/calibrated-data-p2828/231_v2.h5", "r")
        used_detectors = [0, 1, 2, 3, 6, 8, 10, 15]
        W_axis = Array(fid["energy_axis"])
        detector_angle = [180. , 156.6, 134.4, 123.3, 112.2,  67.8,  56.7,  45.6,  34.5,
                 11.1,   0. , 348.9, 292.2, 281.1, 247.8, 203.4]
        θ_axis_deg = detector_angle[used_detectors.+1]
        θ_axis = deg2rad.(θ_axis_deg)
        close(fid)
    end

    laser_wavelength = args["laser-wavelength"]
    Tl = laser_wavelength/c/nm_per_au*fs_per_au
    max_tau = args["max-time-resolution"]
    if max_tau < 0.0
        max_tau = 0.5*Tl
    end
    time_resolution = collect(LinRange(0.2f0, max_tau, 10))
    #time_resolution = [0.2f0]

    output = args["output"]

    #W_axis = collect(LinRange(args["start-w"], args["end-w"], args["nbins-w"]))
    #used_detectors = args["used-detectors"]
    #detector_angle = args["detector-angle"]
    #used_detectors = used_detectors .+ 1
    #θ_axis_deg = detector_angle[used_detectors]
    #θ_axis = deg2rad.(θ_axis_deg)

    gas = args["gas"]
    photoelectron_spread = args["photoelectron-spread"]
    omega_spread = args["omega-spread"]

    nbins_omega = args["nbins-omega"]
    nbins_time = args["nbins-time"]
    nbins_Up = args["nbins-Up"]
    simulate(W_axis=W_axis,
             theta_axis=θ_axis,
             laser_wavelength=laser_wavelength,
             time_resolution=time_resolution,
             photoelectron_spread=photoelectron_spread,
             output_filename=output,
             gas=gas,
             omega_spread=omega_spread,
             nbins_omega=nbins_omega,
             nbins_time=nbins_time,
             nbins_Up=nbins_Up,
             )
end

main()
