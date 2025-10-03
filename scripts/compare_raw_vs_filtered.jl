#!/usr/bin/env julia
# scripts/compare_raw_vs_filtered.jl - Comparación directa de datos crudos vs filtrados
using GLMakie, NPZ, JSON3, FilePathsBase, Statistics, Dates, DSP, Printf, IniFile
using FFTW

# ---------- Función para leer datos crudos (BrainVision) ----------
function load_raw_eeg(eeg_file::AbstractString)
    """Carga datos EEG crudos desde archivos BrainVision"""
    vhdr_file = replace(eeg_file, ".eeg" => ".vhdr")
    
    # Leer metadatos del .vhdr usando la misma lógica que plot_raw_traces.jl
    ini = read(Inifile(), String(vhdr_file))
    common = ini.sections["Common Infos"]
    chsec = haskey(ini.sections, "Channel Infos") ? ini.sections["Channel Infos"] : Dict{String,String}()
    
    nchan = parse(Int, common["NumberOfChannels"])
    
    # Obtener formato binario
    binary_section = haskey(ini.sections, "Binary Infos") ? ini.sections["Binary Infos"] : Dict{String,String}()
    binfmt = get(binary_section, "BinaryFormat", "IEEE_FLOAT_32")
    
    # Frecuencia de muestreo
    sintstr = replace(get(common, "SamplingInterval", "2000"), "us"=>"")
    sint_us = parse(Float64, replace(sintstr, "µs"=>""))
    fs = 1e6 / sint_us
    
    # Nombres de canal en orden
    ch_names = String[]
    for k in sort(collect(keys(chsec)))  # Ch1, Ch2, ...
        v = chsec[k]                     # "Fz,,..." -> primer campo
        push!(ch_names, split(v, ",")[1])
    end
    length(ch_names) == nchan || (ch_names = [ "Ch$(i)" for i=1:nchan ])
    
    # Tipo de dato
    T = binfmt == "INT_16" ? Int16 : Float32
    
    # Leer archivo binario
    io = open(eeg_file, "r")
    file_size = Base.filesize(io)
    n_samples = file_size ÷ (nchan * sizeof(T))
    
    data = Array{T}(undef, nchan, n_samples)
    read!(io, data)
    close(io)
    
    # Convertir a Float32 para consistencia
    raw_data = Float32.(data)
    
    return raw_data, fs, ch_names
end

# ---------- Función para calcular PSD ----------
function compute_psd(x::AbstractVector{<:Real}, fs::Real; nperseg::Int=2048)
    """Calcula la densidad espectral de potencia usando FFT simple"""
    
    n = length(x)
    window = hanning(min(nperseg, n))
    
    # Aplicar ventana y calcular FFT
    x_windowed = x[1:min(nperseg, n)] .* window[1:min(nperseg, n)]
    fft_result = fft(x_windowed)
    psd = abs2.(fft_result)
    
    # Solo frecuencias positivas
    n_freq = length(psd) ÷ 2 + 1
    psd = psd[1:n_freq]
    freqs = (0:(n_freq-1)) .* fs / length(x_windowed)
    
    return freqs, psd
end

# ---------- Función para encontrar archivos pareados ----------
function find_paired_files(bids_root::AbstractString, deriv_root::AbstractString)
    """Encuentra archivos crudos y filtrados que corresponden"""
    paired_files = NamedTuple[]
    
    # Buscar archivos filtrados
    filt_files = String[]
    for (root, _, files) in walkdir(deriv_root)
        for f in files
            if endswith(f, "_desc-hp0.5notch50_eeg.npz")
                push!(filt_files, joinpath(root, f))
            end
        end
    end
    
    # Para cada archivo filtrado, buscar el archivo crudo correspondiente
    for filt_file in sort(filt_files)
        # Extraer información del nombre del archivo filtrado
        base_name = replace(basename(filt_file), "_desc-hp0.5notch50_eeg.npz" => "")
        
        # Buscar archivo crudo correspondiente
        raw_file = nothing
        for (root, _, files) in walkdir(bids_root)
            for f in files
                if endswith(f, ".eeg") && contains(f, base_name)
                    raw_file = joinpath(root, f)
                    break
                end
            end
            raw_file !== nothing && break
        end
        
        if raw_file !== nothing
            push!(paired_files, (raw=raw_file, filtered=filt_file))
        end
    end
    
    return paired_files
end

# ---------- Función principal de comparación ----------
function compare_raw_vs_filtered(raw_data::AbstractMatrix{<:Real}, 
                                filtered_data::AbstractMatrix{<:Real},
                                fs::Real, ch_names::Vector{String};
                                ch_idx::Int=1, win_s::Real=10.0, 
                                base_name::AbstractString="")
    """Compara datos crudos vs filtrados y genera métricas"""
    
    nchan, ns_raw = size(raw_data)
    _, ns_filt = size(filtered_data)
    
    # Usar el mínimo de muestras
    ns = min(ns_raw, ns_filt)
    
    # Extraer segmento de tiempo
    i0 = 1
    i1 = min(round(Int, win_s * fs), ns)
    
    raw_seg = raw_data[ch_idx, i0:i1]
    filt_seg = filtered_data[ch_idx, i0:i1]
    t = ((i0-1):(i1-1)) ./ fs
    
    # Calcular PSD para ambos
    freqs_raw, psd_raw = compute_psd(raw_seg, fs)
    freqs_filt, psd_filt = compute_psd(filt_seg, fs)
    
    # Métricas de comparación
    noise_raw = std(raw_seg)
    noise_filt = std(filt_seg)
    noise_reduction = (noise_raw - noise_filt) / noise_raw * 100
    
    # Encontrar frecuencias dominantes
    max_freq_idx_raw = argmax(psd_raw)
    max_freq_idx_filt = argmax(psd_filt)
    dominant_freq_raw = freqs_raw[max_freq_idx_raw]
    dominant_freq_filt = freqs_filt[max_freq_idx_filt]
    
    # Calcular energía en bandas de frecuencia
    alpha_band = (8.0 .<= freqs_raw .<= 13.0)
    beta_band = (13.0 .<= freqs_raw .<= 30.0)
    gamma_band = (30.0 .<= freqs_raw .<= 50.0)
    
    alpha_power_raw = sum(psd_raw[alpha_band])
    beta_power_raw = sum(psd_raw[beta_band])
    gamma_power_raw = sum(psd_raw[gamma_band])
    
    alpha_power_filt = sum(psd_filt[alpha_band])
    beta_power_filt = sum(psd_filt[beta_band])
    gamma_power_filt = sum(psd_filt[gamma_band])
    
    return (t=t, raw=raw_seg, filtered=filt_seg,
            freqs_raw=freqs_raw, psd_raw=psd_raw,
            freqs_filt=freqs_filt, psd_filt=psd_filt,
            noise_raw=noise_raw, noise_filt=noise_filt, noise_reduction=noise_reduction,
            dominant_freq_raw=dominant_freq_raw, dominant_freq_filt=dominant_freq_filt,
            alpha_power_raw=alpha_power_raw, beta_power_raw=beta_power_raw, gamma_power_raw=gamma_power_raw,
            alpha_power_filt=alpha_power_filt, beta_power_filt=beta_power_filt, gamma_power_filt=gamma_power_filt)
end

# ---------- Función para crear plot de comparación ----------
function plot_comparison(comparison_data, ch_name::String, base_name::String;
                         pngpath::Union{Nothing,AbstractString}=nothing)
    """Crea un plot comparativo de datos crudos vs filtrados"""
    
    fig = Figure(size=(1800, 1200), fontsize=12)
    
    # --- Señales temporales ---
    ax1 = Axis(fig[1,1], title="Señales Temporales - $ch_name",
               xlabel="Tiempo (s)", ylabel="Amplitud (μV)")
    
    lines!(ax1, comparison_data.t, comparison_data.raw, 
           label="Crudo", color=:blue, linewidth=1, alpha=0.8)
    lines!(ax1, comparison_data.t, comparison_data.filtered, 
           label="Filtrado", color=:red, linewidth=1.5)
    
    axislegend(ax1, position=(1, 1))
    
    # --- PSD Comparativo ---
    ax2 = Axis(fig[1,2], title="Densidad Espectral de Potencia",
               xlabel="Frecuencia (Hz)", ylabel="PSD (μV²/Hz)")
    
    lines!(ax2, comparison_data.freqs_raw, comparison_data.psd_raw, 
           label="Crudo", color=:blue, linewidth=2, alpha=0.8)
    lines!(ax2, comparison_data.freqs_filt, comparison_data.psd_filt, 
           label="Filtrado", color=:red, linewidth=2)
    
    # Líneas de referencia para filtros
    vlines!(ax2, [0.5, 50], color=:gray, linestyle=:dash, linewidth=1, alpha=0.7)
    text!(ax2, 0.5, maximum([maximum(comparison_data.psd_raw), maximum(comparison_data.psd_filt)])*0.9, 
          text="HP 0.5Hz", align=(:center, :bottom), color=:gray)
    text!(ax2, 50, maximum([maximum(comparison_data.psd_raw), maximum(comparison_data.psd_filt)])*0.9, 
          text="Notch 50Hz", align=(:center, :bottom), color=:gray)
    
    xlims!(ax2, 0, 60)
    axislegend(ax2, position=(1, 1))
    
    # --- Métricas de comparación ---
    ax3 = Axis(fig[2,1:2], title="Métricas de Comparación")
    
    metrics = [
        ("Ruido Crudo (σ)", comparison_data.noise_raw),
        ("Ruido Filtrado (σ)", comparison_data.noise_filt),
        ("Reducción de Ruido (%)", comparison_data.noise_reduction),
        ("Frec. Dom. Crudo (Hz)", comparison_data.dominant_freq_raw),
        ("Frec. Dom. Filtrado (Hz)", comparison_data.dominant_freq_filt),
        ("Alpha Crudo", comparison_data.alpha_power_raw),
        ("Alpha Filtrado", comparison_data.alpha_power_filt),
        ("Beta Crudo", comparison_data.beta_power_raw),
        ("Beta Filtrado", comparison_data.beta_power_filt),
        ("Gamma Crudo", comparison_data.gamma_power_raw),
        ("Gamma Filtrado", comparison_data.gamma_power_filt)
    ]
    
    labels = [m[1] for m in metrics]
    values = [m[2] for m in metrics]
    
    # Normalizar valores para visualización
    norm_values = values ./ maximum(abs.(values))
    colors = repeat([:blue, :red], 6)[1:length(metrics)]
    
    barplot!(ax3, 1:length(metrics), norm_values, color=colors, alpha=0.7)
    ax3.xticks = (1:length(metrics), labels)
    ax3.yticks = (Float64[], String[])
    
    # Agregar valores reales como texto
    for (i, (label, value)) in enumerate(metrics)
        text!(ax3, i, norm_values[i] + 0.1,
              text=@sprintf("%.2f", value),
              align=(:center, :bottom), fontsize=10)
    end
    
    # Título general
    suptitle = "Comparación Crudo vs Filtrado: $base_name - Canal $ch_name"
    Label(fig[0, :], suptitle, fontsize=14, font=:bold)
    
    # Guardar si se especifica ruta
    if pngpath !== nothing
        mkpath(dirname(String(pngpath)))
        save(String(pngpath), fig)
        @info "✓ Plot comparativo guardado: $(basename(String(pngpath)))"
    end
    
    return fig
end

# ---------- Función principal ----------
function main(bids_root::AbstractString="bids", deriv_root::AbstractString="derivatives/preproc";
              limit::Union{Nothing,Int}=nothing, ch_idx::Int=1, win_s::Real=10.0)
    """Función principal para comparar datos crudos vs filtrados"""
    
    project_root = dirname(dirname(abspath(@__FILE__)))
    plots_dir = joinpath(project_root, "derivatives", "qc", "filtering_comparison")
    mkpath(plots_dir)
    
    @info "Buscando archivos pareados..."
    paired_files = find_paired_files(bids_root, deriv_root)
    isempty(paired_files) && error("No se encontraron archivos pareados")
    
    @info "Encontrados $(length(paired_files)) archivos pareados"
    
    files_to_process = isnothing(limit) ? paired_files : paired_files[1:min(limit, length(paired_files))]
    
    for (i, (raw_file, filt_file)) in enumerate(files_to_process)
        @info "Procesando ($i/$(length(files_to_process))): $(basename(raw_file))"
        
        try
            # Cargar datos crudos
            raw_data, fs_raw, ch_names_raw = load_raw_eeg(raw_file)
            
            # Cargar datos filtrados
            filt_data = NPZ.npzread(filt_file)
            meta_filt = JSON3.read(read(replace(filt_file, ".npz" => ".json")), Dict)
            ch_names_filt = String.(meta_filt["ChannelNames"])
            fs_filt = meta_filt["SamplingFrequency"]
            
            # Verificar compatibilidad
            if fs_raw != fs_filt
                @warn "Frecuencias de muestreo no coinciden: $fs_raw vs $fs_filt"
            end
            
            if length(ch_names_raw) != length(ch_names_filt)
                @warn "Número de canales no coincide: $(length(ch_names_raw)) vs $(length(ch_names_filt))"
            end
            
            # Comparar datos
            comparison = compare_raw_vs_filtered(raw_data, filt_data, fs_raw, ch_names_raw;
                                               ch_idx=ch_idx, win_s=win_s)
            
            # Generar plot
            ch_name = ch_names_raw[ch_idx]
            base_name = replace(basename(raw_file), "_eeg.eeg" => "")
            png_name = "$(base_name)_ch$(ch_idx)_$(ch_name)_comparison.png"
            png_path = joinpath(plots_dir, png_name)
            
            plot_comparison(comparison, ch_name, base_name; pngpath=png_path)
            
            @info "✓ Completado: $base_name - Reducción de ruido: $(round(comparison.noise_reduction, digits=1))%"
            
        catch e
            @warn "Error procesando $(basename(raw_file)): $e"
        end
    end
    
    @info "Completado. Plots comparativos guardados en: $plots_dir"
end

# ---------- CLI entry point ----------
if abspath(PROGRAM_FILE) == @__FILE__
    bids_root = "bids"
    deriv_root = "derivatives/preproc"
    limit = 3  # Por defecto procesar solo 3 archivos
    ch_idx = 1
    win_s = 10.0
    
    # Parsear argumentos
    for arg in ARGS
        if startswith(arg, "--limit=")
            global limit = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--channel=")
            global ch_idx = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--window=")
            global win_s = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--bids=")
            global bids_root = split(arg, "=")[2]
        elseif startswith(arg, "--deriv=")
            global deriv_root = split(arg, "=")[2]
        elseif !startswith(arg, "--")
            # Si no es un flag, asumir que es el directorio BIDS
            global bids_root = arg
        end
    end
    
    main(bids_root, deriv_root; limit=limit, ch_idx=ch_idx, win_s=win_s)
end
