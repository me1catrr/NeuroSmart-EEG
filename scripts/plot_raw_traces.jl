#!/usr/bin/env julia
# scripts/plot_raw_traces.jl
using GLMakie, IniFile, FilePathsBase, DataFrames, CSV, Dates, Statistics
import JSON3
try; import EDF; has_edf=true; catch; has_edf=false; end

# ---------- Utilidades ----------
struct EEGRaw
    data::Matrix{Float32}   # (chan × samples)
    fs::Float64
    chlabels::Vector{String}
    info::Dict{String,Any}
end

function read_brainvision(vhdr_path::AbstractString)::EEGRaw
    ini = read(Inifile(), String(vhdr_path))
    
    # Acceder correctamente a las secciones del INI
    common = ini.sections["Common Infos"]
    chsec = haskey(ini.sections, "Channel Infos") ? ini.sections["Channel Infos"] : Dict{String,String}()

    # El archivo de datos debe tener el mismo nombre base que el .vhdr pero con extensión .eeg
    base_name = replace(basename(String(vhdr_path)), ".vhdr" => "")
    datfile = normpath(joinpath(dirname(String(vhdr_path)), base_name * ".eeg"))
    nchan   = parse(Int, common["NumberOfChannels"])
    orient  = get(common, "DataOrientation", "MULTIPLEXED")
    
    # Obtener formato binario de la sección Binary Infos
    binary_section = haskey(ini.sections, "Binary Infos") ? ini.sections["Binary Infos"] : Dict{String,String}()
    binfmt  = get(binary_section, "BinaryFormat", "IEEE_FLOAT_32")
    
    sintstr = replace(get(common, "SamplingInterval", "2000"), "us"=>"")
    sint_us = parse(Float64, replace(sintstr, "µs"=>""))
    fs      = 1e6 / sint_us

    # nombres de canal en orden
    labels = String[]
    for k in sort(collect(keys(chsec)))  # Ch1, Ch2, ...
        v = chsec[k]                     # "Fz,,..." -> primer campo
        push!(labels, split(v, ",")[1])
    end
    length(labels) == nchan || (labels = [ "Ch$(i)" for i=1:nchan ])

    # tipo de dato
    T = binfmt == "INT_16" ? Int16 : Float32

    # lee binario y reordena
    io = open(datfile, "r")
    file_size = Base.filesize(io)
    nels = div(file_size, sizeof(T))
    nsamp = div(nels, nchan)
    seekstart(io)
    raw = read!(io, Vector{T}(undef, nsamp*nchan))
    close(io)

    if orient == "MULTIPLEXED"
        X = reshape(raw, (nchan, nsamp))    # columnas = tiempo
    else
        X = permutedims(reshape(raw, (nsamp, nchan)))  # VECTORIZED
    end

    # escala a Float32 (si era Int16, lo dejas en unidades de ADC;
    # para BrainVision típico IEEE_FLOAT_32 ya está en µV)
    X32 = Float32.(X)
    info = Dict("path"=>String(vhdr_path), "format"=>"BrainVision")
    return EEGRaw(X32, fs, labels, info)
end

function read_edf(edf_path::AbstractString)::EEGRaw
    @assert has_edf "Instala EDF.jl para leer EDF (Pkg.add(\"EDF\"))."
    e = EDF.read(String(edf_path))
    fs = e.signals[1].sampling_frequency
    chlabels = [s.label for s in e.signals]
    data = reduce(hcat, e.signals) |> x->permutedims(hcat([Float32.(s.data) for s in e.signals]...))
    # EDF.read no siempre devuelve formato homogéneo según versión; forzamos chan×samples
    X = reduce(hcat, [Float32.(e.signals[i].data) for i in 1:length(e.signals)])
    X = permutedims(X)  # chan × samples
    return EEGRaw(X, fs, chlabels, Dict("path"=>String(edf_path), "format"=>"EDF"))
end

function find_bids_eeg(bidsroot::AbstractString)
    paths = String[]
    for (root, _, files) in walkdir(bidsroot)
        for f in files
            # Buscar archivos EEG siguiendo convenciones BIDS
            if endswith(f, "_eeg.vhdr") || endswith(f, "_eeg.edf") || 
               (endswith(f, ".vhdr") && contains(f, "sub-") && contains(f, "ses-")) ||
               (endswith(f, ".edf") && contains(f, "sub-") && contains(f, "ses-"))
                push!(paths, joinpath(root, f))
            end
        end
    end
    sort(paths)
end

# Función para leer metadatos BIDS
function read_bids_metadata(bidsroot::AbstractString)
    # Leer participants.tsv
    participants_file = joinpath(bidsroot, "participants.tsv")
    participants_df = DataFrame()
    if isfile(participants_file)
        participants_df = CSV.read(participants_file, DataFrame)
    end
    
    # Leer dataset_description.json
    dataset_desc_file = joinpath(bidsroot, "dataset_description.json")
    dataset_info = Dict()
    if isfile(dataset_desc_file)
        dataset_info = JSON3.read(read(dataset_desc_file), Dict)
    end
    
    return participants_df, dataset_info
end

# Función para extraer información del nombre de archivo BIDS
function parse_bids_filename(filename::AbstractString)
    parts = split(basename(filename), "_")
    info = Dict{String,String}()
    
    for part in parts
        if startswith(part, "sub-")
            info["subject"] = part
        elseif startswith(part, "ses-")
            info["session"] = part
        elseif startswith(part, "task-")
            info["task"] = part
        elseif startswith(part, "run-")
            info["run"] = part
        end
    end
    
    # Extraer extensión y tipo de archivo
    if endswith(filename, ".vhdr")
        info["format"] = "BrainVision"
        info["extension"] = ".vhdr"
    elseif endswith(filename, ".edf")
        info["format"] = "EDF"
        info["extension"] = ".edf"
    end
    
    return info
end

# ---------- Función para crear layout de electrodos ----------
function create_electrode_layout(chlabels::Vector{String}, outpath::AbstractString)
    # Coordenadas aproximadas de electrodos 10-20 estándar
    electrode_coords = Dict(
        "Fp1" => (-0.35, 0.75), "Fp2" => (0.35, 0.75),
        "Fz" => (0, 0.6), "F3" => (-0.25, 0.5), "F4" => (0.25, 0.5),
        "F7" => (-0.5, 0.4), "F8" => (0.5, 0.4),
        "FC1" => (-0.15, 0.3), "FC2" => (0.15, 0.3),
        "FC5" => (-0.35, 0.25), "FC6" => (0.35, 0.25),
        "FT9" => (-0.6, 0.3), "FT10" => (0.6, 0.3),
        "Cz" => (0, 0.1), "C3" => (-0.25, 0.1), "C4" => (0.25, 0.1),
        "T7" => (-0.5, 0), "T8" => (0.5, 0),
        "CP1" => (-0.15, -0.1), "CP2" => (0.15, -0.1),
        "CP5" => (-0.35, -0.15), "CP6" => (0.35, -0.15),
        "TP9" => (-0.6, -0.1), "TP10" => (0.6, -0.1),
        "Pz" => (0, -0.3), "P3" => (-0.25, -0.3), "P4" => (0.25, -0.3),
        "P7" => (-0.5, -0.4), "P8" => (0.5, -0.4),
        "O1" => (-0.25, -0.6), "O2" => (0.25, -0.6), "Oz" => (0, -0.6)
    )
    
    fig = Figure(resolution=(800, 800), fontsize=12)
    ax = Axis(fig[1,1], title="Layout de Electrodos EEG - Sistema 10-20", 
              xlabel="", ylabel="", aspect=DataAspect())
    
    # Dibujar cabeza
    θ = range(0, 2π, length=100)
    head_x = cos.(θ) * 0.8
    head_y = sin.(θ) * 0.8
    lines!(ax, head_x, head_y, linewidth=3, color=:black)
    
    # Dibujar nariz
    nose_x = [0, -0.05, 0.05]
    nose_y = [0.8, 0.75, 0.75]
    lines!(ax, nose_x, nose_y, linewidth=2, color=:black)
    
    # Dibujar orejas
    ear_left_x = [-0.8, -0.85, -0.8]
    ear_left_y = [0.1, 0, -0.1]
    ear_right_x = [0.8, 0.85, 0.8]
    ear_right_y = [0.1, 0, -0.1]
    lines!(ax, ear_left_x, ear_left_y, linewidth=2, color=:black)
    lines!(ax, ear_right_x, ear_right_y, linewidth=2, color=:black)
    
    # Dibujar electrodos
    used_electrodes = String[]
    for label in chlabels
        if haskey(electrode_coords, label)
            x, y = electrode_coords[label]
            scatter!(ax, [x], [y], markersize=20, color=:red, strokewidth=2, strokecolor=:darkred)
            text!(ax, x, y-0.08, text=label, align=(:center, :center), fontsize=10, color=:black)
            push!(used_electrodes, label)
        end
    end
    
    # Líneas de referencia
    vlines!(ax, 0, color=:gray, linewidth=1, linestyle=:dash)
    hlines!(ax, 0, color=:gray, linewidth=1, linestyle=:dash)
    
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xticksvisible = false
    ax.yticksvisible = false
    
    # Información adicional
    info_text = "Total electrodos: $(length(used_electrodes))\nElectrodos encontrados: $(join(used_electrodes, ", "))"
    text!(ax, 0, -0.9, text=info_text, align=(:center, :center), fontsize=10, color=:gray)
    
    save(outpath, fig)
    @info "Layout de electrodos guardado en: $(outpath)"
    return fig
end

# ---------- Plot apilado tipo Analyzer ----------
function stacked_traces(figtitle::AbstractString, X::AbstractMatrix{<:Real}, fs::Real, chlabels::Vector{String};
                        t_start::Real=0.0, win_s::Real=10.0, step_uV::Union{Nothing,Real}=nothing,
                        grid_s::Real=1.0, pngpath::Union{Nothing,AbstractString}=nothing)

    nchan, ns = size(X)
    i0 = clamp(round(Int, t_start*fs)+1, 1, ns)
    i1 = clamp(i0 + round(Int, win_s*fs) - 1, 1, ns)
    seg = @view X[:, i0:i1]
    t = ((i0-1):(i1-1)) ./ fs

    # escala vertical automática: mediana de std × factor
    sds = mapslices(x->std(x; corrected=false), seg; dims=2)[:]
    Δ = isnothing(step_uV) ? 6*median(sds) : Float64(step_uV)  # 6*sd ~ cómodo para 500 Hz/µV
    offsets = collect(0:Δ:Δ*(nchan-1))

    fig = Figure(resolution=(1400, 800), fontsize=12)
    ax  = Axis(fig[1,1], title=figtitle, xlabel="Tiempo (s)", ylabel="Canal (offset)", yreversed=true)

    for ch in 1:nchan
        y = seg[ch, :] .+ offsets[ch]
        lines!(ax, t, y)
        text!(ax,  t[1], offsets[ch], text=chlabels[ch], align=(:right,:center), offset=(-8,0))
    end

    # rejilla vertical cada grid_s
    for τ in t[1]:grid_s:t[end]
        vlines!(ax, τ, color=:gray, linewidth=0.5)
    end

    hidespines!(ax, :l)
    ax.yticksvisible = false
    ax.yticklabelsvisible = false
    fig[1,1] = ax

    if pngpath !== nothing
        mkpath(dirname(String(pngpath)))
        save(String(pngpath), fig)
    end
    return fig
end

# ---------- Runner ----------
function main(bidsroot::AbstractString; win_s=10.0, grid_s=1.0, step_uV=nothing, limit::Union{Nothing,Int}=nothing)
    files = find_bids_eeg(bidsroot)
    isempty(files) && error("No se encontraron archivos EEG BIDS en $(bidsroot)")

    # Leer metadatos BIDS
    participants_df, dataset_info = read_bids_metadata(bidsroot)
    
    @info "Encontrados $(length(files)) registros EEG" files
    @info "Dataset: $(get(dataset_info, "Name", "No especificado"))"
    @info "Versión BIDS: $(get(dataset_info, "BIDSVersion", "No especificada"))"
    
    n = isnothing(limit) ? length(files) : min(limit, length(files))
    # La carpeta derivatives debe estar en la raíz del proyecto
    project_root = dirname(dirname(bidsroot))  # NeuroSmart-EEG/
    outroot = joinpath(project_root, "derivatives", "qc", "raw_preview")
    mkpath(outroot)

    # Generar layout de electrodos una sola vez
    layout_generated = false
    electrode_layout_path = joinpath(outroot, "electrode_layout.png")

    for (k, f) in enumerate(files[1:n])
        # Parsear información BIDS del archivo
        bids_info = parse_bids_filename(f)
        
        # Leer datos EEG
        rec = endswith(f, ".edf") ? read_edf(f) : read_brainvision(f)
        
        # Generar nombre de archivo siguiendo convenciones BIDS
        base = replace(basename(f), r"\.(vhdr|edf)$" => "")
        png_name = "$(base)_desc-rawpreview_t0-$(Int(win_s))s.png"
        png = joinpath(outroot, png_name)
        
        # Crear título más informativo
        subject_id = get(bids_info, "subject", "unknown")
        session_id = get(bids_info, "session", "unknown")
        task_id = get(bids_info, "task", "unknown")
        
        # Buscar información del participante si está disponible
        participant_info = ""
        if !isempty(participants_df) && hasproperty(participants_df, :participant_id)
            subject_row = findfirst(x -> x == subject_id, participants_df.participant_id)
            if !isnothing(subject_row)
                row = participants_df[subject_row, :]
                if hasproperty(participants_df, :sex)
                    participant_info = " — Sexo: $(row.sex)"
                end
                if hasproperty(participants_df, :age)
                    participant_info *= " — Edad: $(row.age)"
                end
                if hasproperty(participants_df, :group)
                    participant_info *= " — Grupo: $(row.group)"
                end
            end
        end
        
        title = "$(rec.info["format"]) — $subject_id $session_id $task_id$participant_info — fs=$(round(rec.fs, digits=2)) Hz, $(length(rec.chlabels)) chans"
        
        # Generar plot
        stacked_traces(title, rec.data, rec.fs, rec.chlabels; 
                      t_start=0, win_s=win_s, grid_s=grid_s, step_uV=step_uV, pngpath=png)
        @info "Guardado ($k/$n)" png
        
        # Generar layout de electrodos una sola vez
        if !layout_generated
            create_electrode_layout(rec.chlabels, electrode_layout_path)
            layout_generated = true
        end
    end
    
    # Generar reporte resumen
    generate_qc_report(bidsroot, files[1:n], outroot, dataset_info)
    
    @info "Listo. PNGs guardados en: $(outroot)"
end

# Función para generar reporte de calidad
function generate_qc_report(bidsroot::AbstractString, files::Vector{String}, outroot::AbstractString, dataset_info::Dict)
    report_file = joinpath(outroot, "qc_report.txt")
    
    open(report_file, "w") do io
        println(io, "=== REPORTE DE CALIDAD DE DATOS EEG ===")
        println(io, "Dataset: $(get(dataset_info, "Name", "No especificado"))")
        println(io, "Fecha: $(now())")
        println(io, "Total archivos procesados: $(length(files))")
        println(io, "")
        println(io, "Archivos procesados:")
        
        for (i, f) in enumerate(files)
            bids_info = parse_bids_filename(f)
            println(io, "$i. $(basename(f))")
            println(io, "   - Sujeto: $(get(bids_info, "subject", "unknown"))")
            println(io, "   - Sesión: $(get(bids_info, "session", "unknown"))")
            println(io, "   - Tarea: $(get(bids_info, "task", "unknown"))")
            println(io, "   - Formato: $(get(bids_info, "format", "unknown"))")
        end
    end
    
    @info "Reporte de calidad guardado en: $(report_file)"
end

# ejecutar desde CLI: julia --project scripts/plot_raw_traces.jl bids/
if abspath(PROGRAM_FILE) == @__FILE__
    bidsroot = length(ARGS) >= 1 ? ARGS[1] : "bids"
    main(bidsroot)
end