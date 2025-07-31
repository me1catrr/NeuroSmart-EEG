using DSP
using FileIO
using BSON
using Printf
using Dates

# === CONFIGURACIÓN ===
raw_base = "data/raw"
patients_dir = joinpath(raw_base, "patients")
controls_dir = joinpath(raw_base, "controls")
preprocessed_dir = "data/preprocessed"

fs_target = 512.0               # Hz
lowcut, highcut = 0.5, 45.0     # Hz (bandpass)
notch_freq = 50.0               # Hz

essential_exts = [".vhdr", ".eeg", ".vmrk"]

println("🚀 Iniciando preprocesamiento EEG...")
println("Filtros: $(lowcut)-$(highcut) Hz | Notch: $(notch_freq) Hz | Resample: $(fs_target) Hz\n")

# === FUNCIONES ===

"🧹 Borra el contenido previo de la carpeta preprocessed."
function clear_directory(path::String)
    if isdir(path)
        for f in readdir(path)
            rm(joinpath(path, f); recursive=true, force=true)
        end
    else
        mkpath(path)
    end
end
clear_directory(preprocessed_dir)

"📖 Lee la cabecera .vhdr y extrae parámetros básicos."
function read_vhdr(path::String)
    fs, n_channels, gain = 0.0, 0, 1.0
    open(path, "r") do io
        for line in eachline(io)
            if startswith(line, "SamplingInterval")
                val = split(line, "=")[2]
                fs = 1_000_000.0 / parse(Float64, strip(val)) # µs → Hz
            elseif startswith(line, "NumberOfChannels")
                n_channels = parse(Int, split(line, "=")[2])
            elseif startswith(line, "Ch") && occursin(",", line)
                vals = split(split(line, "=")[2], ",")
                if length(vals) >= 3
                    try
                        gain = parse(Float64, vals[3])
                    catch
                        gain = 1.0
                    end
                end
            end
        end
    end
    return fs, n_channels, gain
end

"📂 Carga datos binarios EEG (.eeg) como matriz (canales x muestras)."
function read_eeg(path::String, n_channels::Int)
    raw = read(path)
    n_samples = div(length(raw) ÷ 4, n_channels) # float32 = 4 bytes
    data = reinterpret(Float32, raw)
    return reshape(data, n_channels, n_samples)
end

"🎛 Aplica filtrado bandpass + notch canal por canal."
function filter_eeg(data::Matrix{Float64}, fs::Float64)
    nyquist = fs / 2
    bp = digitalfilter(Bandpass(lowcut/nyquist, highcut/nyquist), Butterworth(4))
    notch = digitalfilter(Bandstop((notch_freq-1)/nyquist, (notch_freq+1)/nyquist), Butterworth(2))

    filtered = [filt(notch, filt(bp, data[ch, :])) for ch in 1:size(data, 1)]
    result = reduce(vcat, filtered)
    return reshape(result, size(data, 1), length(result) ÷ size(data, 1))
end

"📉 Resamplea la señal a fs_target."
function resample_eeg(data::Matrix{Float64}, fs_original::Float64, fs_new::Float64)
    factor = fs_new / fs_original
    resampled = [DSP.resample(data[ch, :], factor) for ch in 1:size(data, 1)]
    result = reduce(vcat, resampled)
    return reshape(result, size(data, 1), length(result) ÷ size(data, 1))
end

"📝 Genera el log limpio."
function write_log(log_path::String; subject, session, fs, n_channels, gain, size_orig, size_final, status, save_path)
    open(log_path, "w") do log
        println(log, "=== Preprocesamiento EEG ===")
        println(log, "Fecha: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(log, "Sujeto: $subject | Sesión: $session\n")
        println(log, "FS original: $(fs) Hz | Canales: $n_channels | Ganancia: $(gain) µV")
        println(log, "Dimensiones originales: $(size_orig)")
        println(log, "Dimensiones tras resampleo: $(size_final)")
        println(log, "\nEstado: $status")
        if status == "✅ OK"
            println(log, "Guardado en: $save_path")
        end
    end
end

# === PROCESAMIENTO ===
for group_dir in [patients_dir, controls_dir]
    group_name = splitpath(group_dir)[end]
    subjects = sort(readdir(group_dir))

    for subject in subjects
        subject_path = joinpath(group_dir, subject)
        sessions_root = joinpath(subject_path, "EEG")

        println("\n→ Sujeto: $subject ($group_name)")
        if !isdir(sessions_root)
            println("   ⚠ Sin carpeta EEG. Saltando.")
            continue
        end

        for session in sort(readdir(sessions_root))
            session_path = joinpath(sessions_root, session)
            files = readdir(session_path)
            save_dir = joinpath(preprocessed_dir, group_name, subject, session)
            mkpath(save_dir)
            log_file = joinpath(save_dir, "$(subject)_$(session)_preprocessing_log.txt")

            println("   ➤ Sesión: $session")

            # Validación esenciales
            missing = [ext for ext in essential_exts if isempty(filter(f -> endswith(f, ext), files))]
            if !isempty(missing)
                println("     ❌ Faltan esenciales: $(join(missing, ", "))")
                write_log(log_file; subject, session, fs="N/A", n_channels="N/A", gain="N/A", size_orig="N/A", size_final="N/A", status="❌ Faltan archivos", save_path="")
                continue
            end

            vhdr_file = filter(f -> endswith(f, ".vhdr"), files)[1]
            eeg_file = filter(f -> endswith(f, ".eeg"), files)[1]

            vhdr_path = joinpath(session_path, vhdr_file)
            eeg_path = joinpath(session_path, eeg_file)

            try
                println("     🔍 Leyendo cabecera...")
                fs, n_channels, gain = read_vhdr(vhdr_path)
                println("     ✅ FS=$(fs) Hz | Canales=$n_channels | Ganancia=$(gain) µV")

                println("     📂 Cargando datos...")
                data = read_eeg(eeg_path, n_channels)
                println("     ✅ Dimensiones: $(size(data))")

                println("     🔍 Filtrando...")
                data_filt = filter_eeg(Matrix{Float64}(data), fs)

                println("     🔍 Resampleando a $(fs_target) Hz...")
                data_final = resample_eeg(data_filt, fs, fs_target)

                save_path = joinpath(save_dir, "EEG_preprocessed.bson")
                BSON.@save save_path data_final fs_target lowcut highcut notch_freq
                println("     💾 Guardado en $save_path")

                write_log(log_file; subject, session, fs, n_channels, gain,
                        size_orig=size(data), size_final=size(data_final),
                        status="✅ OK", save_path)
            catch e
                short_msg = sprint(showerror, e)
                println("     ❌ Error procesando $subject: $short_msg")
                write_log(log_file; subject, session, fs="Error", n_channels="Error", gain="Error",
                        size_orig="N/A", size_final="N/A", status="❌ ERROR: $short_msg", save_path="")
            end
        end
    end
end

println("\n✅ Preprocesamiento completado.")