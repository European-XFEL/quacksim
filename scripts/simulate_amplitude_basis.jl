#!/usr/bin/env julia

using ArgParse
using HDF5
using LaTeXStrings
using Plots
using Format
import H5Zzstd

using QuackSim

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
        "--fel-energy"
            arg_type = Float64
            required = true
            action = :store_arg
            default = 1000.0
            help = "FEL set mea energy, in eV."
        "--ellipticity"
            arg_type = Float64
            required = false
            action = :store_arg
            default = 1.0
            help = "Ellipticity."
        "--min-W"
            arg_type = Float64
            required = false
            action = :store_arg
            help = "Minimum photo-electron kinetic energy, in eV."
        "--max-W"
            arg_type = Float64
            required = false
            action = :store_arg
            help = "Maximum photo-electron kinetic energy, in eV."
        "--nbins-W"
            arg_type = Int64
            required = false
            action = :store_arg
            help = "Number of bins in the photo-electron kinetic energy axis."
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
        "--gas"
            arg_type = String
            required = false
            action = :store_arg
            default = "Ne1s"
            metavar = "GAS"
            help = "Gas and orbital."
        "--max-Up"
            arg_type = Float64
            metavar = "PONDEROMOTIVE_POTENTIAL"
            help = "Maximum ponderomotive potential, in eV."
        "--nbins-time"
            arg_type = Int64
            default = 64
            metavar = "BINS"
            help = "Number of bins in time."
        "--nbins-Up"
            arg_type = Int64
            default = 40
            metavar = "BINS"
            help = "Number of bins of ponderomotive potential."
        "--oversampling"
            arg_type = Int64
            default = 1
            metavar = "INTEGER"
            help = "How many times to oversample energy axis."
        "--beta"
            arg_type = Float64
            default = 2.0
            metavar = "BETA"
            help = "Beta value for the angular distribution."
        "--vertical"
            action = :store_true
            help = "Vertical polarization?"
        "--use-gpu"
            action = :store_true
            help = "Use GPU?"
    end
    args = parse_args(ARGS, s)

    # extra input
    if args["min-W"] != nothing && args["max-W"] != nothing && args["nbins-W"] != nothing
        minW = args["min-W"]
        maxW = args["max-W"]
        nbinsW = args["nbins-W"]
        W_axis = collect(LinRange(minW, maxW, nbinsW))
        θ_axis = collect(range(0, 2π, step=2π/16.0))[1:16]
    elseif args["read-axes-from"] != nothing
        fid = h5open(args["read-axes-from"], "r")
        # simulated data
        if haskey(fid, "W_axis")
            W_axis = read(fid["W_axis"])
            θ_axis = read(fid["theta_axis"])
        else
            # EuXFEL PES
            W_axis = read(fid["energy_axis"])
            #θ_axis = collect(LinRange(0, 2π, 16))
            θ_axis = collect(range(0, 2π, step=2π/16.0))[1:16]
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

    use_gpu = false
    if args["use-gpu"] != nothing
        use_gpu = true
    end
    output = args["output"]
    epsilon = args["ellipticity"]

    gas = args["gas"]
    laser_wavelength = args["laser-wavelength"]
    maxUp = args["max-Up"]

    nbins_time = args["nbins-time"]
    nbins_Up = args["nbins-Up"]
    oversampling = args["oversampling"]
    beta = args["beta"]
    polarization=:horizontal
    if args["vertical"] != nothing
        if args["vertical"]
            polarization = :vertical
        end
    end
    fel_energy = args["fel-energy"]
    simulate_amplitude(
                       fel_energy=fel_energy,
                       W_axis=W_axis,
                       theta_axis=θ_axis,
                       laser_wavelength=laser_wavelength,
                       maxUp=maxUp,
                       output_filename=output,
                       gas=gas,
                       nbins_time=nbins_time,
                       nbins_Up=nbins_Up,
                       oversampling=oversampling,
                       use_gpu=use_gpu,
                       beta=beta,
                       polarization=polarization,
                       ellipticity=epsilon,
                       )
end

main()

