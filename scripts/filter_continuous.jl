#!/usr/bin/env julia

using DSP
using IniFile            # ] add IniFile
using Glob              # ] add Glob
using JSON3             # ] add JSON3
using NPZ               # ] add NPZ
using Printf, Dates
using FilePathsBase: Path, joinpath, parent, basename, splitext
using Statistics
using GLMakie           # ] add GLMakie

# ---------------- Filters ----------------

function make_bandpass(fs; hp=0.5, lp=150.0, order=8)
    @assert lp < fs/2 "LP debe ser < fs/2 (Nyquist)."
    # Usar frecuencias normalizadas (0-1)
    hp_norm = hp / (fs/2)
    lp_norm = lp / (fs/2)
    digitalfilter(Bandpass(hp_norm, lp_norm), Butterworth(order))
end

function make_notch(fs; f0=50.0, Q=35, order=2)
    bw = f0 / Q
    f1 = max(0.1, f0 - bw/2)
    f2 = min(fs/2 - 1e-3, f0 + bw/2)
    # Usar frecuencias normalizadas (0-1)
    f1_norm = f1 / (fs/2)
    f2_norm = f2 / (fs/2)
    digitalfilter(Bandstop(f1_norm, f2_norm), Butterworth(order))
end

function filtfilt_channels(filt, X::AbstractMatrix)
    C, T = size(X)
    Y = similar(X)
    @inbounds for c in 1:C
        Y[c, :] = filtfilt(filt, @view X[c, :])
    end
    return Y
end

# ---------------- BIDS helpers ----------------

"""
    load_eeg(path::AbstractString) -> data::Matrix{Float32}, fs::Float64, ch_names::Vector{String}

Lee continuo desde archivos BrainVision (.eeg/.vhdr/.vmrk) (canales × muestras).
Basado en la implementación de plot_raw_traces.jl
"""
function load_eeg(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    
    if ext == ".eeg"
        # Para archivos BrainVision - leer el .vhdr correspondiente
        vhdr_path = replace(path, ".eeg" => ".vhdr")
        @assert isfile(vhdr_path) "Archivo .vhdr no encontrado: $vhdr_path"
        
        # Leer metadatos desde .vhdr usando IniFile
        ini = read(Inifile(), String(vhdr_path))
        
        # Acceder a las secciones del INI
        common = ini.sections["Common Infos"]
        chsec = haskey(ini.sections, "Channel Infos") ? ini.sections["Channel Infos"] : Dict{String,String}()
        
        # Obtener información básica
        nchan = parse(Int, common["NumberOfChannels"])
        orient = get(common, "DataOrientation", "MULTIPLEXED")
        
        # Obtener formato binario
        binary_section = haskey(ini.sections, "Binary Infos") ? ini.sections["Binary Infos"] : Dict{String,String}()
        binfmt = get(binary_section, "BinaryFormat", "IEEE_FLOAT_32")
        
        # Calcular frecuencia de muestreo
        sintstr = replace(get(common, "SamplingInterval", "2000"), "us"=>"")
        sint_us = parse(Float64, replace(sintstr, "µs"=>""))
        fs = 1e6 / sint_us
        
        # Extraer nombres de canales
        labels = String[]
        for k in sort(collect(keys(chsec)))  # Ch1, Ch2, ...
            v = chsec[k]                     # "Fz,,..." -> primer campo
            push!(labels, split(v, ",")[1])
        end
        length(labels) == nchan || (labels = [ "Ch$(i)" for i=1:nchan ])
        
        # Tipo de dato
        T = binfmt == "INT_16" ? Int16 : Float32
        
        # Leer archivo binario
        io = open(path, "r")
        file_size = Base.filesize(io)
        nels = div(file_size, sizeof(T))
        nsamp = div(nels, nchan)
        seekstart(io)
        raw = read!(io, Vector{T}(undef, nsamp*nchan))
        close(io)
        
        # Reordenar datos según orientación
        if orient == "MULTIPLEXED"
            X = reshape(raw, (nchan, nsamp))    # columnas = tiempo
        else
            X = permutedims(reshape(raw, (nsamp, nchan)))  # VECTORIZED
        end
        
        # Convertir a Float32
        X32 = Float32.(X)
        
        return X32, float(fs), labels
    else
        error("Formato no soportado: $ext. Use .eeg")
    end
end

function bids_entities(fname::String)
    # sub-XX_ses-YY_task-ZZZ_run-1_eeg.edf -> Dict("sub"=>"XX", "ses"=>"YY", "task"=>"ZZZ", "run"=>"1")
    ent = Dict{String,String}()
    base = splitext(basename(fname))[1]
    for part in split(base, '_')
        if occursin('-', part)
            k, v = split(part, '-')
            ent[k] = v
        end
    end
    return ent
end

function out_paths(root_deriv::String, ents::Dict, desc::String)
    sub = get(ents, "sub", "NA")
    ses = get(ents, "ses", "1")
    rel = joinpath("derivatives", "preproc", "sub-$sub", "ses-$ses", "eeg")
    outdir = joinpath(root_deriv, rel)
    mkpath(outdir)
    stem = joinpath(outdir,
        @sprintf("sub-%s_ses-%s%s_desc-%s_eeg",
                 sub, ses,
                 get(ents,"task","")=="" ? "" : "_task-"*ents["task"],
                 desc))
    return stem * ".npz", stem * ".json"
end

# ---------------- Pipeline ----------------

"""
Aplica HP+notch en continuo (fase cero). No aplica LP final.
"""
function filter_continuous!(infile::String; hp=0.5, lp_for_ica=150.0, order_bp=8,
                            notch_f0=50.0, notch_Q=35, notch_order=2,
                            root_deriv::String=".")

    X, fs, ch_names = load_eeg(infile)

    @assert fs > 2lp_for_ica "fs=$(fs) es insuficiente para LP=$(lp_for_ica) Hz."
    bp = make_bandpass(fs; hp=hp, lp=lp_for_ica, order=order_bp)
    X_bp = filtfilt_channels(bp, X)

    nz = make_notch(fs; f0=notch_f0, Q=notch_Q, order=notch_order)
    X_f = filtfilt_channels(nz, X_bp)

    ents = bids_entities(infile)
    # Usar la ruta del proyecto NeuroSmart-EEG como root_deriv
    project_root = dirname(dirname(abspath(@__FILE__)))
    stem_npz, stem_json = out_paths(project_root, ents, @sprintf("hp%0.1fnotch%0.0f", hp, notch_f0))

    # Guarda matriz como NPZ (canales × muestras) - solo datos EEG
    NPZ.npzwrite(stem_npz, X_f)
    meta = Dict(
        "BIDSVersion" => "1.9.0",
        "GeneratedBy" => [Dict("Name"=>"neurosmart-eeg", "Version"=>"0.1", "Description"=>"HP + Notch pre-epoch")],
        "SamplingFrequency" => fs,
        "EEGChannelCount" => size(X_f,1),
        "ChannelNames" => ch_names,
        "EEGReference" => "as-recorded (update if re-referenced)",
        "Filters" => Dict(
            "Highpass" => Dict("Type"=>"Butterworth", "Order"=>order_bp, "CutoffHz"=>hp, "ZeroPhase"=>true),
            "Lowpass_temp" => Dict("Type"=>"Butterworth", "Order"=>order_bp, "CutoffHz"=>lp_for_ica, "ZeroPhase"=>true, "Purpose"=>"stabilize ICA / pre-epoch"),
            "Notch" => Dict("Type"=>"Butterworth bandstop", "Order"=>notch_order, "FrequencyHz"=>notch_f0, "Q"=>notch_Q, "ZeroPhase"=>true)
        ),
        "DateTime" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    )
    open(stem_json, "w") do io
        JSON3.write(io, meta; allow_inf=true)
    end

    return stem_npz, stem_json, size(X_f)
end

# ---------------- Funciones de visualización ----------------

"""
    plot_filtered_traces(X::AbstractMatrix, fs::Float64, ch_names::Vector{String}, 
                         title::String; t_start::Real=0.0, win_s::Real=10.0, 
                         step_uV::Union{Nothing,Real}=nothing, grid_s::Real=1.0, 
                         pngpath::Union{Nothing,AbstractString}=nothing)

Genera un plot apilado de las trazas EEG filtradas, similar a plot_raw_traces.jl
"""
function plot_filtered_traces(X::AbstractMatrix, fs::Float64, ch_names::Vector{String}, title::String;
                             t_start::Real=0.0, win_s::Real=10.0, 
                             step_uV::Union{Nothing,Real}=nothing, grid_s::Real=1.0, 
                             pngpath::Union{Nothing,AbstractString}=nothing)
    
    nchan, ns = size(X)
    i0 = clamp(round(Int, t_start*fs)+1, 1, ns)
    i1 = clamp(i0 + round(Int, win_s*fs) - 1, 1, ns)
    seg = @view X[:, i0:i1]
    t = ((i0-1):(i1-1)) ./ fs

    # Escala vertical automática: mediana de std × factor
    sds = mapslices(x->std(x; corrected=false), seg; dims=2)[:]
    Δ = isnothing(step_uV) ? 6*median(sds) : Float64(step_uV)
    offsets = collect(0:Δ:Δ*(nchan-1))

    fig = Figure(resolution=(1400, 800), fontsize=12)
    ax = Axis(fig[1,1], title=title, xlabel="Tiempo (s)", ylabel="Canal (offset)", yreversed=true)

    for ch in 1:nchan
        y = seg[ch, :] .+ offsets[ch]
        lines!(ax, t, y)
        text!(ax, t[1], offsets[ch], text=ch_names[ch], align=(:right,:center), offset=(-8,0))
    end

    # Rejilla vertical cada grid_s
    for τ in t[1]:grid_s:t[end]
        vlines!(ax, τ, color=:gray, linewidth=0.5)
    end

    hidespines!(ax, :l)
    ax.yticksvisible = false
    ax.yticklabelsvisible = false

    if pngpath !== nothing
        mkpath(dirname(String(pngpath)))
        save(String(pngpath), fig)
    end
    
    return fig
end

"""
    generate_filtered_preview(npz_path::String, json_path::String, 
                             outdir::String; win_s::Real=10.0)

Genera una visualización de los datos filtrados y la guarda como PNG
"""
function generate_filtered_preview(npz_path::String, json_path::String, outdir::String; win_s::Real=10.0)
    # Leer datos filtrados
    X_f = NPZ.npzread(npz_path)
    
    # Leer metadatos
    meta = JSON3.read(read(json_path), Dict)
    fs = meta["SamplingFrequency"]
    ch_names = meta["ChannelNames"]
    
    # Extraer información del nombre de archivo
    base_name = splitext(basename(npz_path))[1]
    png_name = "$(base_name)_desc-filteredpreview_t0-$(Int(win_s))s.png"
    png_path = joinpath(outdir, png_name)
    
    # Crear título descriptivo
    title = "FILTRADO (HP 0.5Hz + Notch 50Hz) — $base_name — fs=$(round(fs, digits=2)) Hz, $(length(ch_names)) chans"
    
    # Generar plot
    plot_filtered_traces(X_f, fs, ch_names, title; t_start=0, win_s=win_s, pngpath=png_path)
    
    return png_path
end

# ---------------- CLI mínima ----------------

function find_bids_eeg_files(bidsroot::AbstractString)
    """Encuentra archivos EEG en estructura BIDS, similar a plot_raw_traces.jl"""
    paths = String[]
    for (root, _, files) in walkdir(bidsroot)
        for f in files
            # Buscar archivos .eeg siguiendo convenciones BIDS
            if endswith(f, "_eeg.eeg") || 
               (endswith(f, ".eeg") && contains(f, "sub-") && contains(f, "ses-"))
                push!(paths, joinpath(root, f))
            end
        end
    end
    return sort(paths)
end

function main(; generate_plots::Bool=true, win_s::Real=10.0)
    root = abspath(get(ENV, "BIDS_ROOT", "bids"))
    files = find_bids_eeg_files(root)
    isempty(files) && error("No se encontraron archivos EEG en $root")

    # Crear directorio para plots filtrados
    project_root = dirname(dirname(abspath(@__FILE__)))
    plots_dir = joinpath(project_root, "derivatives", "qc", "filtered_preview")
    if generate_plots
        mkpath(plots_dir)
        @info "Generando plots en: $plots_dir"
    end

    @info "Encontrados $(length(files)) archivos EEG"
    for f in sort(files)
        @info "Filtrando" f
        npz, js, sz = filter_continuous!(f; hp=0.5, lp_for_ica=150.0, order_bp=8,
                                         notch_f0=50.0, notch_Q=35, notch_order=2)
        @info "Guardado" npz js sz
        
        # Generar preview si está habilitado
        if generate_plots
            try
                png_path = generate_filtered_preview(npz, js, plots_dir; win_s=win_s)
                @info "Plot generado: $(basename(png_path))"
            catch e
                @warn "Error generando plot para $(basename(npz)): $e"
            end
        end
    end
    
    if generate_plots
        @info "Plots de datos filtrados guardados en: $plots_dir"
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Permitir argumentos de línea de comandos
    generate_plots = true
    win_s = 10.0
    
    # Parsear argumentos simples
    for arg in ARGS
        if arg == "--no-plots"
            generate_plots = false
        elseif startswith(arg, "--window=")
            win_s = parse(Float64, split(arg, "=")[2])
        end
    end
    
    main(generate_plots=generate_plots, win_s=win_s)
end