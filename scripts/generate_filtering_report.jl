#!/usr/bin/env julia
# scripts/generate_filtering_report.jl - Genera reporte resumen del análisis de filtrado
using JSON3, FilePathsBase, Statistics, Dates, Printf

# ---------- Función para extraer métricas de archivos de análisis ----------
function extract_metrics_from_analysis(analysis_dir::AbstractString)
    """Extrae métricas de todos los archivos de análisis"""
    
    metrics_list = NamedTuple[]
    
    for (root, _, files) in walkdir(analysis_dir)
        for f in files
            if endswith(f, ".png") && contains(f, "_analysis")
                # Extraer información del nombre del archivo
                base_name = replace(f, ".png" => "")
                
                # Parsear información del archivo (ej: sub-M10_ses-1_task-eyesopen_ch1_Fz_analysis)
                try
                    parts = split(base_name, "_")
                    if length(parts) >= 6  # sub-M10_ses-1_task-eyesopen_ch1_Fz_analysis
                        subject = parts[1] * "_" * parts[2]  # sub-M10_ses-1
                        task = parts[3]  # task-eyesopen
                        channel_idx_str = parts[4]  # ch1
                        channel_name = parts[5]  # Fz
                        
                        # Extraer número del canal de forma más robusta
                        channel_idx = 1  # Valor por defecto
                        if startswith(channel_idx_str, "ch")
                            try
                                channel_idx = parse(Int, replace(channel_idx_str, "ch" => ""))
                            catch
                                channel_idx = 1  # Valor por defecto si falla el parseo
                            end
                        end
                        
                        push!(metrics_list, (
                            file = f,
                            subject = subject,
                            task = task,
                            channel_idx = channel_idx,
                            channel_name = channel_name,
                            base_name = base_name
                        ))
                    end
                catch e
                    @warn "Error parsing file $f: $e"
                end
            end
        end
    end
    
    return sort(metrics_list, by=x->x.base_name)
end

# ---------- Función para generar reporte HTML ----------
function generate_html_report(metrics_list::Vector, output_path::AbstractString)
    """Genera un reporte HTML con todas las métricas"""
    
    html_content = """
    <!DOCTYPE html>
    <html lang="es">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Reporte de Análisis de Filtrado EEG</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
            .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
            h1 { color: #333; text-align: center; border-bottom: 3px solid #4CAF50; padding-bottom: 10px; }
            h2 { color: #555; margin-top: 30px; }
            .summary { background-color: #e8f5e8; padding: 15px; border-radius: 5px; margin: 20px 0; }
            .metrics-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
            .metrics-table th, .metrics-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            .metrics-table th { background-color: #4CAF50; color: white; }
            .metrics-table tr:nth-child(even) { background-color: #f2f2f2; }
            .metrics-table tr:hover { background-color: #e8f5e8; }
            .channel-name { font-weight: bold; color: #2196F3; }
            .subject { color: #FF9800; }
            .task { color: #9C27B0; }
            .plot-link { color: #4CAF50; text-decoration: none; }
            .plot-link:hover { text-decoration: underline; }
            .footer { text-align: center; margin-top: 30px; color: #666; font-size: 0.9em; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧠 Reporte de Análisis de Filtrado EEG</h1>
            
            <div class="summary">
                <h2>📊 Resumen General</h2>
                <p><strong>Total de archivos analizados:</strong> $(length(metrics_list))</p>
                <p><strong>Fecha de generación:</strong> $(Dates.format(now(), "dd/mm/yyyy HH:MM:SS"))</p>
                <p><strong>Filtros aplicados:</strong> Pasa-alto 0.5 Hz + Notch 50 Hz (fase cero)</p>
            </div>
            
            <h2>📋 Archivos Analizados</h2>
            <table class="metrics-table">
                <thead>
                    <tr>
                        <th>Sujeto</th>
                        <th>Tarea</th>
                        <th>Canal</th>
                        <th>Nombre del Canal</th>
                        <th>Archivo de Análisis</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    for metric in metrics_list
        html_content *= """
                    <tr>
                        <td class="subject">$(metric.subject)</td>
                        <td class="task">$(metric.task)</td>
                        <td>$(metric.channel_idx)</td>
                        <td class="channel-name">$(metric.channel_name)</td>
                        <td><a href="../qc/filtering_analysis/$(metric.file)" class="plot-link" target="_blank">$(metric.file)</a></td>
                    </tr>
        """
    end
    
    html_content *= """
                </tbody>
            </table>
            
            <h2>📁 Estructura de Resultados</h2>
            <ul>
                <li><strong>Datos filtrados:</strong> <code>derivatives/preproc/</code></li>
                <li><strong>Análisis de filtrado:</strong> <code>derivatives/qc/filtering_analysis/</code></li>
                <li><strong>Comparaciones crudo vs filtrado:</strong> <code>derivatives/qc/filtering_comparison/</code></li>
                <li><strong>Previsualizaciones filtradas:</strong> <code>derivatives/qc/filtered_preview/</code></li>
                <li><strong>Previsualizaciones crudas:</strong> <code>derivatives/qc/raw_preview/</code></li>
            </ul>
            
            <h2>🔧 Scripts Utilizados</h2>
            <ul>
                <li><strong>Filtrado:</strong> <code>scripts/filter_continuous.jl</code></li>
                <li><strong>Análisis de filtrado:</strong> <code>scripts/compare_filtering_effects.jl</code></li>
                <li><strong>Comparación crudo vs filtrado:</strong> <code>scripts/compare_raw_vs_filtered.jl</code></li>
                <li><strong>Previsualizaciones filtradas:</strong> <code>scripts/plot_filtered_traces.jl</code></li>
            </ul>
            
            <div class="footer">
                <p>Reporte generado automáticamente por NeuroSmart-EEG</p>
                <p>Para más información, consulte la documentación del proyecto</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    write(output_path, html_content)
    @info "✓ Reporte HTML generado: $output_path"
end

# ---------- Función para generar reporte CSV ----------
function generate_csv_report(metrics_list::Vector, output_path::AbstractString)
    """Genera un reporte CSV con todas las métricas"""
    
    csv_content = "Sujeto,Tarea,Canal,Nombre_Canal,Archivo_Análisis\n"
    
    for metric in metrics_list
        csv_content *= "$(metric.subject),$(metric.task),$(metric.channel_idx),$(metric.channel_name),$(metric.file)\n"
    end
    
    write(output_path, csv_content)
    @info "✓ Reporte CSV generado: $output_path"
end

# ---------- Función principal ----------
function main(analysis_dir::AbstractString="derivatives/qc/filtering_analysis")
    """Función principal para generar reportes"""
    
    project_root = dirname(dirname(abspath(@__FILE__)))
    reports_dir = joinpath(project_root, "derivatives", "reports")
    mkpath(reports_dir)
    
    # Usar ruta absoluta para el directorio de análisis
    if !isabspath(analysis_dir)
        analysis_dir = joinpath(project_root, analysis_dir)
    end
    
    @info "Extrayendo métricas de análisis desde: $analysis_dir"
    metrics_list = extract_metrics_from_analysis(analysis_dir)
    
    if isempty(metrics_list)
        @warn "No se encontraron archivos de análisis en: $analysis_dir"
        return
    end
    
    @info "Encontrados $(length(metrics_list)) archivos de análisis"
    
    # Generar reporte HTML
    html_path = joinpath(reports_dir, "filtering_analysis_report.html")
    generate_html_report(metrics_list, html_path)
    
    # Generar reporte CSV
    csv_path = joinpath(reports_dir, "filtering_analysis_report.csv")
    generate_csv_report(metrics_list, csv_path)
    
    @info "Reportes generados en: $reports_dir"
    @info "Abre el reporte HTML en tu navegador: file://$html_path"
end

# ---------- CLI entry point ----------
if abspath(PROGRAM_FILE) == @__FILE__
    analysis_dir = "derivatives/qc/filtering_analysis"
    
    # Parsear argumentos
    for arg in ARGS
        if startswith(arg, "--analysis-dir=")
            global analysis_dir = split(arg, "=")[2]
        elseif !startswith(arg, "--")
            # Si no es un flag, asumir que es el directorio de análisis
            global analysis_dir = arg
        end
    end
    
    main(analysis_dir)
end
