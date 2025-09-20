#!/usr/bin/env julia

using JSON
using DataFrames
using CSV
using Glob

function validate_bids_structure(bids_dir)
    println("🔍 Validando estructura BIDS...")
    
    errors = []
    warnings = []
    
    # Verificar archivos obligatorios en la raíz
    required_files = ["dataset_description.json", "participants.tsv", "participants.json", "README"]
    for file in required_files
        if !isfile(joinpath(bids_dir, file))
            push!(errors, "❌ Archivo obligatorio faltante: $file")
        else
            println("✅ $file encontrado")
        end
    end
    
    # Verificar dataset_description.json
    if isfile(joinpath(bids_dir, "dataset_description.json"))
        try
            dataset_desc = JSON.parsefile(joinpath(bids_dir, "dataset_description.json"))
            required_fields = ["Name", "BIDSVersion", "DatasetType"]
            for field in required_fields
                if !haskey(dataset_desc, field)
                    push!(errors, "❌ Campo obligatorio faltante en dataset_description.json: $field")
                end
            end
            println("✅ dataset_description.json válido")
        catch e
            push!(errors, "❌ Error parseando dataset_description.json: $e")
        end
    end
    
    # Verificar participants.tsv
    if isfile(joinpath(bids_dir, "participants.tsv"))
        try
            participants = CSV.read(joinpath(bids_dir, "participants.tsv"), DataFrame)
            if !("participant_id" in names(participants))
                push!(errors, "❌ Columna 'participant_id' faltante en participants.tsv")
            else
                println("✅ participants.tsv válido ($(nrow(participants)) participantes)")
            end
        catch e
            push!(errors, "❌ Error leyendo participants.tsv: $e")
        end
    end
    
    # Verificar estructura de participantes
    participant_dirs = filter(x -> startswith(x, "sub-"), readdir(bids_dir))
    println("📁 Encontrados $(length(participant_dirs)) directorios de participantes")
    
    for sub_dir in participant_dirs
        sub_path = joinpath(bids_dir, sub_dir)
        
        # Verificar que sea un directorio
        if !isdir(sub_path)
            push!(warnings, "⚠️ $sub_dir no es un directorio")
            continue
        end
        
        # Verificar sesiones
        session_dirs = filter(x -> startswith(x, "ses-"), readdir(sub_path))
        if isempty(session_dirs)
            push!(warnings, "⚠️ $sub_dir no tiene sesiones")
            continue
        end
        
        for ses_dir in session_dirs
            ses_path = joinpath(sub_path, ses_dir)
            eeg_path = joinpath(ses_path, "eeg")
            
            if !isdir(eeg_path)
                push!(warnings, "⚠️ $sub_dir/$ses_dir no tiene directorio eeg/")
                continue
            end
            
            # Verificar archivos EEG
            eeg_files = readdir(eeg_path)
            eeg_json_files = filter(x -> endswith(x, "_eeg.json"), eeg_files)
            eeg_data_files = filter(x -> endswith(x, "_eeg.eeg") || endswith(x, "_eeg.vhdr"), eeg_files)
            
            if isempty(eeg_json_files)
                push!(warnings, "⚠️ $sub_dir/$ses_dir/eeg/ no tiene archivos _eeg.json")
            end
            
            if isempty(eeg_data_files)
                push!(warnings, "⚠️ $sub_dir/$ses_dir/eeg/ no tiene archivos de datos EEG")
            end
            
            # Validar archivos JSON individuales
            for json_file in eeg_json_files
                json_path = joinpath(eeg_path, json_file)
                try
                    eeg_json = JSON.parsefile(json_path)
                    
                    # Verificar campos obligatorios
                    required_eeg_fields = ["TaskName", "SamplingFrequency", "EEGChannelCount"]
                    for field in required_eeg_fields
                        if !haskey(eeg_json, field)
                            push!(warnings, "⚠️ $json_file falta campo: $field")
                        end
                    end
                    
                    # Verificar que el nombre del archivo coincida con el TaskName
                    if haskey(eeg_json, "TaskName")
                        task_name = eeg_json["TaskName"]
                        if !occursin(task_name, json_file)
                            push!(warnings, "⚠️ $json_file: TaskName '$task_name' no coincide con el nombre del archivo")
                        end
                    end
                    
                catch e
                    push!(errors, "❌ Error parseando $json_file: $e")
                end
            end
        end
    end
    
    # Verificar archivos de tareas
    task_files = ["task-eyesopen_eeg.json", "task-eyesclosed_eeg.json"]
    for task_file in task_files
        if isfile(joinpath(bids_dir, task_file))
            println("✅ $task_file encontrado")
        else
            push!(warnings, "⚠️ Archivo de tarea faltante: $task_file")
        end
    end
    
    return errors, warnings
end

function main()
    bids_dir = "bids"
    
    if !isdir(bids_dir)
        println("❌ Directorio BIDS no encontrado: $bids_dir")
        return
    end
    
    println("🚀 Iniciando validación BIDS para: $bids_dir")
    println("="^60)
    
    errors, warnings = validate_bids_structure(bids_dir)
    
    println("\n" * "="^60)
    println("📊 RESUMEN DE VALIDACIÓN")
    println("="^60)
    
    if isempty(errors)
        println("✅ ¡No se encontraron errores críticos!")
    else
        println("❌ ERRORES CRÍTICOS ENCONTRADOS:")
        for error in errors
            println("   $error")
        end
    end
    
    if !isempty(warnings)
        println("\n⚠️ ADVERTENCIAS:")
        for warning in warnings
            println("   $warning")
        end
    end
    
    println("\n📈 ESTADÍSTICAS:")
    println("   • Errores críticos: $(length(errors))")
    println("   • Advertencias: $(length(warnings))")
    
    if isempty(errors)
        println("\n🎉 ¡El dataset BIDS es válido!")
    else
        println("\n🔧 Se requieren correcciones antes de usar el dataset.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
