using CSV
using DataFrames
using JSON
using Printf
using FilePathsBase

# === CONFIGURACIÓN ===
csv_path = "data/raw/info/demográficos Rafael MIND EM_pac vs cont t1.csv"
base_all_data_dir = "data/raw/all_data"
patients_dir = "data/raw/patients"
controls_dir = "data/raw/controls"
log_path = "data/raw/info/organize/organize_log.txt"

essential_exts = ["vhdr", "eeg", "vmrk"]       # Requeridos
expected_exts = ["vhdr", "eeg", "vmrk", "ehst2", "hfinf2"]  # Todos

# === LIMPIEZA INICIAL ===
println("\n🧹 Limpiando directorios y log previos...")
function clear_directory(path::String)
    if isdir(path)
        for f in readdir(path)
            rm(joinpath(path, f); recursive=true, force=true)
        end
    else
        mkpath(path)
    end
end
clear_directory(patients_dir)
clear_directory(controls_dir)
if isfile(log_path)
    rm(log_path; force=true)
end
println("✅ Limpieza completada.\n")

# Variables de resumen
total_subjects = 0
total_sessions = 0
complete_sessions = 0
partial_sessions = 0
incomplete_sessions = 0

open(log_path, "w") do log
    println("\n📄 Leyendo CSV...")
    df = CSV.read(csv_path, DataFrame; delim=';', ignorerepeated=true)
    rename!(df, names(df) .=> [lowercase(replace(n, " " => "_")) for n in names(df)])
    println("✅ CSV cargado: $(size(df, 1)) filas")

    for row in eachrow(df)
        global total_subjects += 1
        codigo = strip(row[:codigo])
        tipo = strip(row[:tipo])
        subject_base_dir = tipo == "PAC" ? patients_dir : controls_dir
        subject_dir = joinpath(subject_base_dir, codigo)
        mkpath(subject_dir)

        println("\n➡ Procesando sujeto: $codigo ($tipo)")
        println(log, "\n[$codigo]")

        # Buscar archivos
        subject_files = filter(f -> startswith(f, codigo * "_"), readdir(base_all_data_dir))
        if isempty(subject_files)
            println("   ⚠ No hay archivos en all_data")
            println(log, "   ⚠ Sin archivos")
            continue
        end

        # Agrupar por sesión
        sessions = Dict{String, Dict}()
        for file in subject_files
            session = if occursin("_T1_", file)
                "T1"
            elseif occursin("_T2_", file)
                "T2"
            else
                "OTROS"
            end
            session_files = get!(sessions, session, Dict("files" => String[], "status" => ""))
            push!(session_files["files"], file)
        end

        for (session, info) in sessions
            global total_sessions += 1
            session_dir = joinpath(subject_dir, "EEG", session)
            mkpath(session_dir)

            println("   ➤ Sesión $session")
            println(log, "   $session:")

            found_exts = String[]
            for file in info["files"]
                src = joinpath(base_all_data_dir, file)
                condition = if occursin("COG", file)
                    "COG"
                elseif occursin("OJOS CERRADOS", file)
                    "OJOS_CERRADOS"
                elseif occursin("OJOS ABIERTOS", file)
                    "OJOS_ABIERTOS"
                else
                    "UNKNOWN"
                end

                ext = splitext(file)[2]
                push!(found_exts, replace(ext, "." => ""))
                new_name = "$(codigo)_$(session)_$(condition)$(ext)"
                dest = joinpath(session_dir, new_name)

                cp(src, dest; force=true)
                println("      ✅ $file → $new_name")
            end

            # Validación
            missing_essentials = setdiff(essential_exts, found_exts)
            missing_all = setdiff(expected_exts, found_exts)

            if isempty(missing_essentials) && isempty(missing_all)
                info["status"] = "completa"
                global complete_sessions += 1
                println("      ✅ Sesión completa")
                println(log, "      OK")
            elseif isempty(missing_essentials) && !isempty(missing_all)
                info["status"] = "parcial"
                global partial_sessions += 1
                println("      ⚠ Sesión parcial (faltan extras: $(join(missing_all, ", ")))")
                println(log, "      PARCIAL (faltan extras: $(join(missing_all, ", ")))")
            else
                info["status"] = "incompleta"
                global incomplete_sessions += 1
                println("      ❌ Sesión incompleta (faltan esenciales: $(join(missing_essentials, ", ")))")
                println(log, "      INCOMPLETA (faltan esenciales: $(join(missing_essentials, ", ")))")
                # Renombrar carpeta con sufijo _INCOMPLETO
                new_dir = session_dir * "_INCOMPLETO"
                mv(session_dir, new_dir; force=true)
            end
        end

        # Metadata por sujeto
        metadata = Dict(
            "codigo" => codigo,
            "tipo" => tipo,
            "sesiones" => Dict(k => Dict("status" => v["status"], "archivos" => v["files"]) for (k, v) in sessions)
        )
        open(joinpath(subject_dir, "metadata.json"), "w") do io
            JSON.print(io, metadata)
        end
    end

    # === Resumen final ===
    println("\n=== RESUMEN FINAL ===")
    println("Sujetos: $total_subjects")
    println("Sesiones: $total_sessions")
    println("Completas: $complete_sessions")
    println("Parciales: $partial_sessions")
    println("Incompletas: $incomplete_sessions")
    if total_sessions > 0
        println(@sprintf("Completitud: %.2f%%", (complete_sessions / total_sessions) * 100))
    end

    println(log, "\n=== RESUMEN FINAL ===")
    println(log, "Sujetos: $total_subjects")
    println(log, "Sesiones: $total_sessions")
    println(log, "Completas: $complete_sessions")
    println(log, "Parciales: $partial_sessions")
    println(log, "Incompletas: $incomplete_sessions")
end

println("\n✅ Organización completada. Log en: $log_path")