#!/usr/bin/env julia
# scripts/compare_filtering_effects.jl - Análisis comparativo del filtrado EEG
using GLMakie, NPZ, JSON3, FilePathsBase, Statistics, Dates, DSP, Printf, IniFile
using FFTW  # Para análisis espectral

# ---------- Función para calcular PSD ----------
function compute_psd(x::AbstractVector{<:Real}, fs::Real; nperseg::Int=2048)
    """Calcula la densidad espectral de potencia usando FFT simple"""
    
    # Usar FFT simple para PSD
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

# ---------- Función para analizar señal filtrada ----------
function analyze_filtered_signal(filtered_data::AbstractMatrix{<:Real},
                                fs::Real, ch_names::Vector{String};
                                ch_idx::Int=1, win_s::Real=10.0)
    """Analiza una señal filtrada y muestra sus características"""
    
    nchan, ns = size(filtered_data)
    
    # Extraer segmento de tiempo
    i0 = 1
    i1 = min(round(Int, win_s * fs), ns)
    
    filt_seg = filtered_data[ch_idx, i0:i1]
    t = ((i0-1):(i1-1)) ./ fs
    
    # Calcular PSD
    freqs_filt, psd_filt = compute_psd(filt_seg, fs)
    
    # Métricas de la señal filtrada
    noise_filt = std(filt_seg)
    mean_amplitude = mean(abs.(filt_seg))
    
    # Encontrar picos espectrales
    max_freq_idx = argmax(psd_filt)
    dominant_freq = freqs_filt[max_freq_idx]
    
    # Calcular energía en bandas de frecuencia
    alpha_band = (8.0 .<= freqs_filt .<= 13.0)
    beta_band = (13.0 .<= freqs_filt .<= 30.0)
    gamma_band = (30.0 .<= freqs_filt .<= 50.0)
    
    alpha_power = sum(psd_filt[alpha_band])
    beta_power = sum(psd_filt[beta_band])
    gamma_power = sum(psd_filt[gamma_band])
    
    return (t=t, filt=filt_seg,
            freqs_filt=freqs_filt, psd_filt=psd_filt,
            noise_filt=noise_filt, mean_amplitude=mean_amplitude,
            dominant_freq=dominant_freq,
            alpha_power=alpha_power, beta_power=beta_power, gamma_power=gamma_power)
end

# ---------- Función para crear plot de análisis ----------
function plot_filtering_analysis(analysis_data, ch_name::String, base_name::String;
                                 pngpath::Union{Nothing,AbstractString}=nothing)
    """Crea un plot de análisis de la señal filtrada"""
    
    fig = Figure(size=(1600, 1000), fontsize=12)
    
    # --- Señal temporal ---
    ax1 = Axis(fig[1,1], title="Señal Filtrada (HP 0.5Hz + Notch 50Hz) - $ch_name", 
               xlabel="Tiempo (s)", ylabel="Amplitud (μV)")
    
    lines!(ax1, analysis_data.t, analysis_data.filt, 
           label="Filtrado", color=:red, linewidth=1.5)
    
    axislegend(ax1, position=(1, 1))
    
    # --- Densidad espectral ---
    ax2 = Axis(fig[1,2], title="Densidad Espectral de Potencia", 
               xlabel="Frecuencia (Hz)", ylabel="PSD (μV²/Hz)")
    
    lines!(ax2, analysis_data.freqs_filt, analysis_data.psd_filt, 
           label="Filtrado", color=:red, linewidth=2)
    
    # Líneas de referencia para filtros
    vlines!(ax2, [0.5, 50], color=:gray, linestyle=:dash, linewidth=1, alpha=0.7)
    text!(ax2, 0.5, maximum(analysis_data.psd_filt)*0.9, text="HP 0.5Hz", 
          align=(:center, :bottom), color=:gray)
    text!(ax2, 50, maximum(analysis_data.psd_filt)*0.9, text="Notch 50Hz", 
          align=(:center, :bottom), color=:gray)
    
    # Bandas de frecuencia EEG
    band_colors = [:orange, :purple, :green]
    band_ranges = [(8, 13), (13, 30), (30, 50)]
    band_names = ["Alpha", "Beta", "Gamma"]
    
    for (i, ((f1, f2), color, name)) in enumerate(zip(band_ranges, band_colors, band_names))
        band_mask = (analysis_data.freqs_filt .>= f1) .& (analysis_data.freqs_filt .<= f2)
        if any(band_mask)
            # Usar band! en lugar de bandfill!
            band!(ax2, analysis_data.freqs_filt[band_mask], 
                  zeros(length(analysis_data.freqs_filt[band_mask])), 
                  analysis_data.psd_filt[band_mask], 
                  color=color, alpha=0.3)
            text!(ax2, (f1+f2)/2, maximum(analysis_data.psd_filt)*0.7, text=name, 
                  align=(:center, :center), color=color, font=:bold)
        end
    end
    
    xlims!(ax2, 0, 60)  # Enfocar en frecuencias relevantes
    axislegend(ax2, position=(1, 1))
    
    # --- Métricas de la señal ---
    ax3 = Axis(fig[2,1:2], title="Características de la Señal Filtrada")
    
    metrics = [
        ("Ruido (σ)", analysis_data.noise_filt),
        ("Amplitud Media", analysis_data.mean_amplitude),
        ("Frec. Dominante (Hz)", analysis_data.dominant_freq),
        ("Potencia Alpha", analysis_data.alpha_power),
        ("Potencia Beta", analysis_data.beta_power),
        ("Potencia Gamma", analysis_data.gamma_power)
    ]
    
    labels = [m[1] for m in metrics]
    values = [m[2] for m in metrics]
    
    # Normalizar valores para visualización
    norm_values = values ./ maximum(abs.(values))
    colors = [:blue, :green, :red, :orange, :purple, :brown]
    
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
    suptitle = "Análisis de Señal EEG Filtrada: $base_name - Canal $ch_name"
    Label(fig[0, :], suptitle, fontsize=14, font=:bold)
    
    # Guardar si se especifica ruta
    if pngpath !== nothing
        mkpath(dirname(String(pngpath)))
        save(String(pngpath), fig)
        @info "✓ Plot de análisis guardado: $(basename(String(pngpath)))"
    end
    
    return fig
end

# ---------- Función para encontrar archivos filtrados ----------
function find_filtered_files(deriv_root::AbstractString)
    """Encuentra archivos filtrados para análisis"""
    filt_files = String[]
    for (root, _, files) in walkdir(deriv_root)
        for f in files
            if endswith(f, "_desc-hp0.5notch50_eeg.npz")
                push!(filt_files, joinpath(root, f))
            end
        end
    end
    return sort(filt_files)
end

# ---------- Función para cargar datos raw ----------
function load_raw_eeg(eeg_file::AbstractString)
    """Carga datos EEG raw usando la lógica de plot_raw_traces.jl"""
    vhdr_file = replace(eeg_file, ".eeg" => ".vhdr")
    
    if !isfile(vhdr_file)
        error("Archivo .vhdr no encontrado para $eeg_file")
    end
    
    # Leer configuración usando Inifile
    ini = read(Inifile(), String(vhdr_file))
    
    # Acceder a las secciones
    common = ini.sections["Common Infos"]
    chsec = haskey(ini.sections, "Channel Infos") ? ini.sections["Channel Infos"] : Dict{String,String}()
    binary_section = haskey(ini.sections, "Binary Infos") ? ini.sections["Binary Infos"] : Dict{String,String}()
    
    # Configuración
    nchan = parse(Int, common["NumberOfChannels"])
    orient = get(common, "DataOrientation", "MULTIPLEXED")
    binfmt = get(binary_section, "BinaryFormat", "IEEE_FLOAT_32")
    
    # Frecuencia de muestreo
    sintstr = replace(get(common, "SamplingInterval", "2000"), "us"=>"")
    sint_us = parse(Float64, replace(sintstr, "µs"=>""))
    fs = 1e6 / sint_us
    
    # Nombres de canales
    ch_names = String[]
    for k in sort(collect(keys(chsec)))
        v = chsec[k]
        push!(ch_names, split(v, ",")[1])
    end
    length(ch_names) == nchan || (ch_names = ["Ch$(i)" for i=1:nchan])
    
    # Tipo de dato
    T = binfmt == "INT_16" ? Int16 : Float32
    
    # Leer datos binarios
    io = open(eeg_file, "r")
    file_size = Base.filesize(io)
    nels = div(file_size, sizeof(T))
    nsamp = div(nels, nchan)
    seekstart(io)
    raw = read!(io, Vector{T}(undef, nsamp*nchan))
    close(io)
    
    # Reordenar según orientación
    if orient == "MULTIPLEXED"
        X = reshape(raw, (nchan, nsamp))    # columnas = tiempo
    else
        X = permutedims(reshape(raw, (nsamp, nchan)))  # VECTORIZED
    end
    
    # Convertir a Float64
    X = Float64.(X)
    
    return X, fs, ch_names
end

# ---------- Función principal ----------
function main(deriv_root::AbstractString="derivatives/preproc"; 
              limit::Union{Nothing,Int}=nothing, ch_idx::Int=1, win_s::Real=10.0)
    """Función principal para analizar características del filtrado"""
    
    project_root = dirname(dirname(abspath(@__FILE__)))
    plots_dir = joinpath(project_root, "derivatives", "qc", "filtering_analysis")
    mkpath(plots_dir)
    
    @info "Buscando archivos filtrados..."
    filt_files = find_filtered_files(deriv_root)
    isempty(filt_files) && error("No se encontraron archivos filtrados")
    
    @info "Encontrados $(length(filt_files)) archivos filtrados"
    
    files_to_process = isnothing(limit) ? filt_files : filt_files[1:min(limit, length(filt_files))]
    
    for (i, filt_file) in enumerate(files_to_process)
        @info "Procesando ($i/$(length(files_to_process))): $(basename(filt_file))"
        
        try
            # Cargar datos filtrados
            X_filt = NPZ.npzread(filt_file)
            meta_filt = JSON3.read(read(replace(filt_file, ".npz" => ".json")), Dict)
            ch_names_filt = String.(meta_filt["ChannelNames"])
            fs = meta_filt["SamplingFrequency"]
            
            # Analizar señal
            analysis = analyze_filtered_signal(X_filt, fs, ch_names_filt; 
                                             ch_idx=ch_idx, win_s=win_s)
            
            # Generar plot
            ch_name = ch_names_filt[ch_idx]
            base_name = replace(basename(filt_file), "_desc-hp0.5notch50_eeg.npz" => "")
            png_name = "$(base_name)_ch$(ch_idx)_$(ch_name)_analysis.png"
            png_path = joinpath(plots_dir, png_name)
            
            plot_filtering_analysis(analysis, ch_name, base_name; pngpath=png_path)
            
            @info "✓ Completado: $base_name - Ruido: $(round(analysis.noise_filt, digits=2)) μV, Frec. dominante: $(round(analysis.dominant_freq, digits=1)) Hz"
            
        catch e
            @warn "Error procesando $(basename(filt_file)): $e"
        end
    end
    
    @info "Completado. Plots de análisis guardados en: $plots_dir"
end

# ---------- CLI entry point ----------
if abspath(PROGRAM_FILE) == @__FILE__
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
        elseif !startswith(arg, "--")
            # Si no es un flag, asumir que es el directorio de derivados
            global deriv_root = arg
        end
    end
    
    main(deriv_root; limit=limit, ch_idx=ch_idx, win_s=win_s)
end
