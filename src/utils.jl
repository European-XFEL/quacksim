
"""Show message `msg` including the current task number."""
function log(msg::AbstractString)
    #println(format("[{1:3d}] {2}", Threads.threadid(), msg))
    #println(format("[{1}] {2}", objectid(current_task()), msg))
    time = Dates.format(Dates.now(), "e, dd u yyyy HH:MM:SS")
    println(format("[{1}] {2}", time, msg))
end


function read_tda(filename::String, n::Integer, l::Integer)
    #                       P.E. =  1050.00     n,l -        E,l             R_int
    #                                           1 0 -   192.99 1      -3.13987E-02
    #                                           2 0 -  1006.91 1       6.74669E-03
    #                                           2 1 -  1029.99 0       1.31832E-03
    #                                           2 1 -  1029.99 2       3.18585E-03
    reg_title = look_for(
                           BEGIN * char_in("#") * one_or_more(WHITESPACE) *
                           "P.E." * one_or_more(WHITESPACE) * "=" * one_or_more(WHITESPACE) *
                           ReadableRegex.capture(one_or_more(DEC_DIGIT_NUMBER) * maybe(".") * one_or_more(DEC_DIGIT_NUMBER); as="energy")
                           )
    reg_end = look_for(
                           BEGIN * char_in("#") * maybe(one_or_more(WHITESPACE)) * END
                           )
    reg_data = look_for(
                            BEGIN * char_in("#") * maybe(one_or_more(WHITESPACE)) *
                            ReadableRegex.capture(DEC_DIGIT_NUMBER; as="n") *
                            one_or_more(WHITESPACE) *
                            ReadableRegex.capture(DEC_DIGIT_NUMBER; as="l") *
                            one_or_more(WHITESPACE) * "-" * one_or_more(WHITESPACE) *
                            ReadableRegex.capture(one_or_more(DEC_DIGIT_NUMBER) * maybe(".") * one_or_more(DEC_DIGIT_NUMBER); as="E") *
                            one_or_more(WHITESPACE) *
                            ReadableRegex.capture(DEC_DIGIT_NUMBER; as="tda_l") *
                            one_or_more(WHITESPACE) *
                            ReadableRegex.capture(maybe(either("+", "-")) * one_or_more(DEC_DIGIT_NUMBER) *
                                    maybe(".") * one_or_more(DEC_DIGIT_NUMBER) *
                                    maybe(either("E", "e")) * maybe(either("+", "-")) *
                                    one_or_more(DEC_DIGIT_NUMBER); as="tda")
                        )
    e = Dict()
    tda = Dict()
    tda_int = Dict()
    start_table = false
    start = false
    energy = 0.0
    open(filename) do f
        while !eof(f)
            line = readline(f)
            # ignore initial comments
            if !start_table
                if occursin("# Photoabsorption cross section (in Mb):", line)
                    start_table = true
                end
                # nothing to do so far
                continue
            end
            # if the main data started, process it
            if start_table
                # it has several blocks marked by reg_title
                # check if we found the title so far
                if !start
                    # no title so, check for it
                    m = match(reg_title, line)
                    if m != nothing
                        # title found and we have not seen it so far
                        # get photon energy and mark as the start of section
                        energy = parse(Float64, m["energy"])
                        start = true
                    end
                    # if the title was not found, do nothing
                    # either way, go to the next line, we are done here
                else # a start of block was found previously
                    # check if it ends now
                    m = match(reg_end, line)
                    if m != nothing
                        # it did end, so mark the end of block, and go to the next line
                        start = false
                    else
                        # it did not end, so let's try some actual data!
                        m_data = match(reg_data, line)
                        if m_data != nothing
                            # is it for the orbital we care?
                            found_n, found_l = parse(Int32, m_data["n"]), parse(Int32, m_data["l"])
                            tda_l = parse(Int32, m_data["tda_l"])
                            tda_value = parse(Float64, m_data["tda"])
                            if found_n == n && found_l == l
                                # yes, so save it
                                if tda_l ∉ keys(e)
                                    e[tda_l] = Vector{Float64}()
                                end
                                if tda_l ∉ keys(tda)
                                    tda[tda_l] = Vector{Float64}()
                                end
                                push!(e[tda_l], energy)
                                push!(tda[tda_l], tda_value)
                            end
                        end
                    end
                end
            end
        end
    end
    # interpolate
    tda_int = Dict(p=>Interpolations.linear_interpolation(e[p]/eV_per_au, tda[p], extrapolation_bc=Interpolations.Line()) for p ∈ keys(tda))
    return tda_int
end


"""
Dictionary of species known for the gas.
"""
orbitals = Dict(
            "" => Dict("species"=> "unknown", "n"=>1, "l"=>0, "m"=>0, "Ip"=>870.2, "filename"=>"", "tda"=>"", "orbital"=>"1s"),
            "Ne1s" => Dict("species"=> "Ne", "n"=>1, "l"=>0, "m"=>0, "Ip"=>870.2, "filename"=>"pcs_Ne.txt", "tda"=>"pcs_Ne_v.txt", "orbital"=>"1s"),
            "N1s" => Dict("species"=> "N", "n"=>1, "l"=>0, "m"=>0, "Ip"=>409.9, "filename"=>"pcs_N.txt", "tda"=>"pcs_N_v.txt", "orbital"=>"1s"),
            "Xe3d" => Dict("species"=> "Xe", "n"=>3, "l"=>2, "m"=>0, "Ip"=>689.0, "filename"=>"pcs_Xe.txt", "tda"=>"pcs_Xe3d_v.txt", "orbital"=>"3d"),
            "Kr1s" => Dict("species"=> "Kr", "n"=>1, "l"=>0, "m"=>0, "Ip"=>14.326e3, "filename"=>"pcs_Kr1s.txt", "tda"=>"pcs_Kr1s_v.txt", "orbital"=>"1s"),
            # Xe5p Ip is actually 12.1 eV, but we are assuming multi-photon excitation for a 5 eV laser
            "Xe5p_multi" => Dict("species"=> "Xe", "n"=>5, "l"=>1, "m"=>0, "Ip"=>2.0, "filename"=>"pcs_Xe.txt", "tda"=>"pcs_Xe5p_v.txt", "orbital"=>"5p"),
            # Kr4p Ip is actually 14.1 eV, but we are assuming multi-photon excitation for a 5 eV laser
            "Kr4p_multi" => Dict("species"=> "Kr", "n"=>4, "l"=>1, "m"=>0, "Ip"=>3.2, "filename"=>"pcs_Kr1s.txt", "tda"=>"pcs_Kr1s_v.txt", "orbital"=>"1s"),
            "Xe2p" => Dict("species"=> "Xe", "n"=>2, "l"=>1, "m"=>0, "Ip"=>4776.21, "filename"=>"pcs_Xe2p.txt", "tda"=>"pcs_Xe2p_v.txt", "orbital"=>"2p"),
            "Ar2s" => Dict("species"=> "Ar", "n"=>2, "l"=>0, "m"=>0, "Ip"=>326.0, "filename"=>"pcs_Ar.txt", "tda"=>"pcs_Ar_v.txt", "orbital"=>"2s"),
            "Xe3d52" => Dict("species"=> "Xe", "n"=>3, "l"=>2, "m"=>0, "Ip"=>676.0, "filename"=>"pcs_Xe.txt", "tda"=>"pcs_Xe3d_v.txt", "orbital"=>"3d"),
            "He1s" => Dict("species"=> "He", "n"=>1, "l"=>0, "m"=>0, "Ip"=>24.6, "filename"=>"pcs_He.txt", "tda"=>"pcs_He_v.txt", "orbital"=>"1s"),
           );
"""
Load modules cross section from files in data.
"""
function load_cross_section()
    for name in keys(orbitals)
        if name == ""
            # default for no gas specified
            orbitals[name]["σ"] =  (E -> (E./(E.^3)).^2);
        else
            orbital = orbitals[name]
            current_dir = @__DIR__
            fname = joinpath(current_dir, "..", "data", orbital["filename"])
            log(format("Loading file '{1}' with data for {2}...", fname, name))
            pcs, header = DelimitedFiles.readdlm(fname; header=true, comments=true)
            pe = pcs[:,1]./eV_per_au;
            σ = pcs[:,findfirst(==(orbital["orbital"]), header)[2]];
            # do linear interpolation, so that we can integrate it later
            fit_σ = Interpolations.linear_interpolation(pe, σ, extrapolation_bc=Interpolations.Line());
            orbitals[name]["σ"] = fit_σ;
        end
    end
end

"""
Load transition dipole amplitudes from files in data.
"""
function load_tda()
    for name in keys(orbitals)
        if name == ""
            # default for no gas specified
            orbitals[name]["d"] = (E -> (E./(E.^3)));
        else
            orbital = orbitals[name]
            current_dir = @__DIR__
            fname = joinpath(current_dir, "..", "data", orbital["tda"])
            log(format("Loading file '{1}' with data for {2}...", fname, name))
            h5open(fname, "r") do fid
                p = read(fid["pr"])
                θ = read(fid["ptheta"])
                ϕ = read(fid["pphi"])
                d = read(fid["fHi"])
                # do linear interpolation
                fit_d = Interpolations.linear_interpolation((p, θ, ϕ), d, extrapolation_bc=Interpolations.Line());
                orbitals[name]["d"] = fit_d;
            end
        end
    end
end

"""
Read TDA from XATOM
"""
function load_xatom_tda()
    for name in keys(orbitals)
        if name == ""
            # default for no gas specified
            orbitals[name]["σ"] =  (E -> (E./(E.^3)).^2);
        else
            orbital = orbitals[name]
            current_dir = @__DIR__
            fname = joinpath(current_dir, "..", "data", orbital["tda"])
            log(format("Loading file '{1}' with data for {2}...", fname, name))
            orbitals[name]["tda_int"] = read_tda(fname, orbital["n"], orbital["l"])
        end
    end
end


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

