# Read the XATOM cross section and transition dipole amplitude files in data/,
# sample them on a common photon energy grid and dump the arrays as a Julia
# source file.
#
# Usage: julia --project scripts/xatom_to_julia.jl

using DelimitedFiles: DelimitedFiles
using Interpolations: Interpolations
using ReadableRegex: ReadableRegex, BEGIN, DEC_DIGIT_NUMBER, END, WHITESPACE,
    char_in, either, look_for, maybe, one_or_more

datadir = joinpath(@__DIR__, "..", "data")
output = joinpath(@__DIR__, "..", "src", "xatom_data.jl")

# Photon energies (eV)
energy = 0.0:1.0:3000.0

orbitals = Dict(
    "He1s" => (species="He", n=1, l=0, m=0, Ip=24.6, filename="pcs_He.txt", tda="pcs_He_v.txt", orbital="1s"),
    "Ne1s" => (species="Ne", n=1, l=0, m=0, Ip=870.2, filename="pcs_Ne.txt", tda="pcs_Ne_v.txt", orbital="1s"),
    "N1s" => (species="N", n=1, l=0, m=0, Ip=409.9, filename="pcs_N.txt", tda="pcs_N_v.txt", orbital="1s"),
    "Xe3d" => (species="Xe", n=3, l=2, m=0, Ip=689.0, filename="pcs_Xe.txt", tda="pcs_Xe3d_v.txt", orbital="3d"),
    "C1s" => (species="C", n=1, l=0, m=0, Ip=291.0, filename="pcs_C.txt", tda="pcs_C_v.txt", orbital="1s"),
    "O1s" => (species="O", n=1, l=0, m=0, Ip=537.0, filename="pcs_O.txt", tda="pcs_O_v.txt", orbital="1s"),
    "S1s" => (species="S", n=1, l=0, m=0, Ip=2448.0, filename="pcs_S.txt", tda="pcs_S_v.txt", orbital="1s"),
    "F1s" => (species="F", n=1, l=0, m=0, Ip=688.0, filename="pcs_F.txt", tda="pcs_F_v.txt", orbital="1s"),
    # others, not used so often
    "Kr1s" => (species="Kr", n=1, l=0, m=0, Ip=14.326e3, filename="pcs_Kr1s.txt", tda="pcs_Kr1s_v.txt", orbital="1s"),
    # Xe5p Ip is actually 12.1 eV, but we are assuming multi-photon excitation for a 5 eV laser
    "Xe5p_multi" => (species="Xe", n=5, l=1, m=0, Ip=2.0, filename="pcs_Xe.txt", tda="pcs_Xe5p_v.txt", orbital="5p"),
    # Kr4p Ip is actually 14.1 eV, but we are assuming multi-photon excitation for a 5 eV laser
    "Kr4p_multi" => (species="Kr", n=4, l=1, m=0, Ip=3.2, filename="pcs_Kr1s.txt", tda="pcs_Kr1s_v.txt", orbital="1s"),
    "Xe2p" => (species="Xe", n=2, l=1, m=0, Ip=4776.21, filename="pcs_Xe2p.txt", tda="pcs_Xe2p_v.txt", orbital="2p"),
    "Ar2s" => (species="Ar", n=2, l=0, m=0, Ip=326.0, filename="pcs_Ar.txt", tda="pcs_Ar_v.txt", orbital="2s"),
)

"""
Read the transition dipole amplitudes for orbital `n`, `l` from an XATOM file,
returning the photon energies and values per angular momentum of the continuum state.
"""
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
    e = Dict{Int, Vector{Float64}}()
    tda = Dict{Int, Vector{Float64}}()
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
                continue
            end
            # the table has several blocks, each starting with reg_title and ending with reg_end
            if !start
                m = match(reg_title, line)
                if !isnothing(m)
                    energy = parse(Float64, m["energy"])
                    start = true
                end
            elseif !isnothing(match(reg_end, line))
                start = false
            else
                m_data = match(reg_data, line)
                if !isnothing(m_data)
                    found_n, found_l = parse(Int, m_data["n"]), parse(Int, m_data["l"])
                    if found_n == n && found_l == l
                        tda_l = parse(Int, m_data["tda_l"])
                        push!(get!(e, tda_l, Float64[]), energy)
                        push!(get!(tda, tda_l, Float64[]), parse(Float64, m_data["tda"]))
                    end
                end
            end
        end
    end
    return Dict(p => (; pe=e[p], value=tda[p]) for p in keys(tda))
end

data = Dict{String, Any}()
for (name, orbital) in orbitals
    println("Reading $name...")

    pcs, header = DelimitedFiles.readdlm(joinpath(datadir, orbital.filename); header=true, comments=true)
    σ = pcs[:, findfirst(==(orbital.orbital), header)[2]]
    σ_int = Interpolations.linear_interpolation(pcs[:, 1], σ, extrapolation_bc=Interpolations.Line())

    # The amplitudes are zero outside the calculated range, and are summed
    # over the angular momenta of the continuum state.
    tda = zeros(length(energy))
    for t in values(read_tda(joinpath(datadir, orbital.tda), orbital.n, orbital.l))
        tda .+= Interpolations.linear_interpolation(t.pe, t.value, extrapolation_bc=0).(energy)
    end

    data[name] = (; orbital.species, orbital.orbital, orbital.n, orbital.l, orbital.m, orbital.Ip,
                  tda, cross_section=σ_int.(energy))
end

function write_array(io, x::AbstractVector, indent; per_line=6)
    println(io, "[")
    for chunk in Iterators.partition(x, per_line)
        println(io, indent, "    ", join(repr.(chunk), ", "), ",")
    end
    print(io, indent, "]")
end

function write_dict(io, d::AbstractDict, indent)
    println(io, "Dict(")
    for key in sort(collect(keys(d)))
        print(io, indent, "    ", repr(key), " => ")
        write_value(io, d[key], indent * "    ")
        println(io, ",")
    end
    print(io, indent, ")")
end

function write_tuple(io, t::NamedTuple, indent)
    println(io, "(")
    for (key, value) in pairs(t)
        print(io, indent, "    ", key, " = ")
        write_value(io, value, indent * "    ")
        println(io, ",")
    end
    print(io, indent, ")")
end

function write_value(io, x, indent)
    if x isa AbstractVector
        write_array(io, x, indent)
    elseif x isa AbstractDict
        write_dict(io, x, indent)
    elseif x isa NamedTuple
        write_tuple(io, x, indent)
    else
        print(io, repr(x))
    end
end

open(output, "w") do io
    println(io, "# Generated by scripts/xatom_to_julia.jl from the XATOM files in data/. Do not edit.")
    println(io, "# tda: transition dipole amplitude summed over the continuum angular momenta,")
    println(io, "# cross_section: photoionization cross section (Mb), both sampled at `xatom_energy`.")
    println(io)
    println(io, "\"Photon energies (eV) at which `xatom_data` is sampled.\"")
    println(io, "const xatom_energy = ", repr(energy))
    println(io)
    print(io, "const xatom_data = ")
    write_value(io, data, "")
    println(io)
end

println("Wrote $(length(data)) orbitals to $output")
