#!/usr/bin/env julia
# scripts/plot_filtered_traces.jl - Versión simplificada basada en plot_raw_traces.jl
using GLMakie, NPZ, JSON3, FilePathsBase, Statistics, Dates

# ---------- Función para crear plot apilado ----------
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

    fig = Figure(size=(1400, 800), fontsize=12)
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
        @info "✓ Plot guardado: $(basename(String(pngpath)))"
    end
    return fig
end

# ---------- Buscar archivos filtrados ----------
function find_filtered_eeg_files(deriv_root::AbstractString)
    """Encuentra archivos EEG filtrados (.npz) en la estructura derivatives"""
    paths = String[]
    for (root, _, files) in walkdir(deriv_root)
        for f in files
            if endswith(f, "_desc-hp0.5notch50_eeg.npz")
                push!(paths, joinpath(root, f))
            end
        end
    end
    return sort(paths)
end

# ---------- Función principal ----------
function main(deriv_root::AbstractString; win_s=10.0, grid_s=1.0, step_uV=nothing, limit::Union{Nothing,Int}=nothing)
    files = find_filtered_eeg_files(deriv_root)
    isempty(files) && error("No se encontraron archivos EEG filtrados (.npz) en $(deriv_root)")
    
    @info "Encontrados $(length(files)) archivos EEG filtrados"
    
    # La carpeta derivatives debe estar en la raíz del proyecto
    project_root = dirname(dirname(deriv_root))  # NeuroSmart-EEG/
    outroot = joinpath(project_root, "derivatives", "qc", "filtered_preview")
    mkpath(outroot)
    
    n = isnothing(limit) ? length(files) : min(limit, length(files))
    
    for (k, f_npz) in enumerate(files[1:n])
        f_json = replace(f_npz, ".npz" => ".json")
        if !isfile(f_json)
            @warn "Saltando $(basename(f_npz)): archivo JSON de metadatos no encontrado."
            continue
        end
        
        @info "Procesando ($k/$n)" basename(f_npz)
        
        try
            # Leer datos filtrados
            X_f = NPZ.npzread(f_npz)
            
            # Leer metadatos
            meta = JSON3.read(read(f_json), Dict)
            fs = meta["SamplingFrequency"]
            ch_names = String.(meta["ChannelNames"])  # Convertir a Vector{String}
            
            # Generar nombre de archivo siguiendo convenciones BIDS
            base = replace(basename(f_npz), "_desc-hp0.5notch50_eeg.npz" => "")
            png_name = "$(base)_desc-filtered_t0-$(Int(win_s))s.png"
            png = joinpath(outroot, png_name)
            
            # Crear título descriptivo
            title = "FILTRADO (HP 0.5Hz + Notch 50Hz) — $base — fs=$(round(fs, digits=2)) Hz, $(length(ch_names)) chans"
            
            # Generar plot
            stacked_traces(title, X_f, fs, ch_names; 
                          t_start=0, win_s=win_s, grid_s=grid_s, step_uV=step_uV, pngpath=png)
        catch e
            @warn "Error procesando $(basename(f_npz)): $e"
        end
    end
    
    @info "Completado. PNGs guardados en: $(outroot)"
end

# ejecutar desde CLI: julia --project scripts/plot_filtered_traces.jl derivatives/preproc/
if abspath(PROGRAM_FILE) == @__FILE__
    deriv_root = "derivatives/preproc"
    win_s = 10.0
    limit = nothing
    
    # Parsear argumentos simples
    for arg in ARGS
        if startswith(arg, "--limit=")
            limit = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--window=")
            win_s = parse(Float64, split(arg, "=")[2])
        elseif !startswith(arg, "--")
            deriv_root = arg
        end
    end
    
    main(deriv_root; win_s=win_s, limit=limit)
end