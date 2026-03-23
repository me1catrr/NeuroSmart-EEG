#!/usr/bin/env julia

# -----------------------------------------------------------------------------
# tools/pages/build.jl
#
# Tooling de publicacion para GitHub Pages (no logica cientifica).
#
# Hace todo en una sola pasada:
#   1) Asegura estructura docs/Pluto/<Modulo>/
#   2) Copia HTML exportado desde un directorio opcional (si existe)
#   3) Si no hay exportado, genera placeholder HTML
#   4) Regenera docs/pages_manifest.json
#   5) Imprime resumen final (published vs placeholder)
#
# Uso:
#   julia tools/pages/build.jl
#   julia tools/pages/build.jl /ruta/a/exportaciones_html
#
# Convencion de busqueda en export_dir para cada modulo:
#   - <export_dir>/<Modulo>/index.html
#   - <export_dir>/<Modulo>.html
#   - <export_dir>/<modulo_minusculas>.html
# -----------------------------------------------------------------------------

using Dates

const MODULES = [
    (
        name = "BIDS",
        source = "Pluto/BIDS/BIDS.jl",
        target = "docs/Pluto/BIDS/index.html",
        description = "Estandarizacion BIDS para datasets EEG",
    ),
    (
        name = "Connectivity",
        source = "Pluto/Connectivity/Connectivity.jl",
        target = "docs/Pluto/Connectivity/index.html",
        description = "Analisis de conectividad funcional en EEG",
    ),
    (
        name = "Preprocessing",
        source = "Pluto/Preprocessing/Preprocessing.jl",
        target = "docs/Pluto/Preprocessing/index.html",
        description = "Limpieza y preparacion de senales EEG",
    ),
    (
        name = "Processing",
        source = "Pluto/Processing/Processing.jl",
        target = "docs/Pluto/Processing/index.html",
        description = "Extraccion y analisis de caracteristicas",
    ),
]

function project_root()::String
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function read_text(path::String)::String
    open(path, "r") do io
        return read(io, String)
    end
end

function write_text(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
end

function json_escape(value::String)::String
    escaped = replace(value, "\\" => "\\\\")
    escaped = replace(escaped, "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n")
    escaped = replace(escaped, "\r" => "\\r")
    escaped = replace(escaped, "\t" => "\\t")
    return escaped
end

function placeholder_html(module_name::String, source_path::String, description::String)::String
    generated_at = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    return """
<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NeuroSmart-EEG | $(module_name)</title>
  <style>
    :root {
      --bg: #f4f7fb;
      --surface: #ffffff;
      --text: #1f2937;
      --muted: #4b5563;
      --primary: #1d4ed8;
      --primary-hover: #1e40af;
      --border: #dbe4f0;
      --shadow: 0 10px 24px rgba(16, 24, 40, 0.08);
      --radius: 14px;
      --max-width: 860px;
    }

    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 0;
      background: var(--bg);
      color: var(--text);
      font-family: "Inter", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      line-height: 1.6;
    }
    main {
      min-height: 100dvh;
      display: grid;
      place-items: center;
      padding: 2rem 1rem;
    }
    .card {
      width: min(100%, var(--max-width));
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      box-shadow: var(--shadow);
      padding: 1.6rem;
    }
    h1 {
      margin-top: 0;
      margin-bottom: 0.6rem;
      font-size: clamp(1.5rem, 3vw, 2rem);
      color: #0f172a;
    }
    p {
      margin: 0.5rem 0;
      color: var(--muted);
    }
    .badge {
      display: inline-block;
      margin-bottom: 0.7rem;
      padding: 0.28rem 0.7rem;
      border-radius: 999px;
      background: #eff6ff;
      border: 1px solid #cfe0ff;
      color: #1e3a8a;
      font-size: 0.86rem;
      font-weight: 600;
    }
    .source {
      margin-top: 0.85rem;
      padding: 0.75rem;
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 10px;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 0.88rem;
      color: #334155;
      overflow-x: auto;
    }
    .actions {
      margin-top: 1.2rem;
      display: flex;
      gap: 0.7rem;
      flex-wrap: wrap;
    }
    .btn {
      display: inline-block;
      padding: 0.55rem 0.9rem;
      border-radius: 10px;
      border: 1px solid transparent;
      background: var(--primary);
      color: #fff;
      text-decoration: none;
      font-weight: 600;
      font-size: 0.92rem;
    }
    .btn:hover {
      background: var(--primary-hover);
      text-decoration: none;
    }
    footer {
      margin-top: 1rem;
      color: #64748b;
      font-size: 0.86rem;
    }
  </style>
</head>
<body>
  <main>
    <article class="card">
      <span class="badge">NeuroSmart-EEG</span>
      <h1>Modulo $(module_name)</h1>
      <p><strong>Contenido pendiente de exportacion.</strong></p>
      <p>$(description)</p>
      <p>
        Esta pagina es un placeholder automatico. Cuando exportes el notebook real
        de Pluto, reemplaza este archivo por el HTML final del modulo.
      </p>
      <p>Notebook fuente esperado:</p>
      <div class="source">$(source_path)</div>
      <div class="actions">
        <a class="btn" href="../../index.html">Volver a la portada</a>
      </div>
      <footer>Placeholder generado automaticamente el $(generated_at).</footer>
    </article>
  </main>
</body>
</html>
"""
end

function candidate_export_paths(export_dir::String, module_name::String)::Vector{String}
    lowercase_name = lowercase(module_name)
    return [
        joinpath(export_dir, module_name, "index.html"),
        joinpath(export_dir, "$(module_name).html"),
        joinpath(export_dir, "$(lowercase_name).html"),
    ]
end

function resolve_exported_html(export_dir::Union{Nothing, String}, module_name::String)::Union{Nothing, String}
    export_dir === nothing && return nothing
    for path in candidate_export_paths(export_dir, module_name)
        isfile(path) && return path
    end
    return nothing
end

function relative_to_root(path::String, root::String)::String
    rel = relpath(path, root)
    return replace(rel, "\\" => "/")
end

function build_manifest_entry(name::String, path::String, source::String, published::Bool)::String
    published_text = published ? "true" : "false"
    return string(
        "  {\n",
        "    \"name\": \"", json_escape(name), "\",\n",
        "    \"path\": \"", json_escape(path), "\",\n",
        "    \"source\": \"", json_escape(source), "\",\n",
        "    \"published\": ", published_text, "\n",
        "  }",
    )
end

function main()
    root = project_root()
    docs_dir = joinpath(root, "docs")
    docs_pluto_dir = joinpath(docs_dir, "Pluto")
    mkpath(docs_pluto_dir)

    export_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : nothing
    if export_dir !== nothing && !isdir(export_dir)
        println("Aviso: no existe el directorio de exportaciones: ", export_dir)
        println("Se usaran placeholders para modulos sin HTML local.")
        export_dir = nothing
    end

    println("== NeuroSmart-EEG Pages Build ==")
    println("Raiz del proyecto: ", root)
    export_dir !== nothing && println("Origen opcional de exportaciones: ", export_dir)

    manifest_entries = String[]
    published_modules = String[]
    placeholder_modules = String[]

    for mod in MODULES
        target_html_abs = joinpath(root, mod.target)
        target_dir_abs = dirname(target_html_abs)
        mkpath(target_dir_abs)

        source_export_abs = resolve_exported_html(export_dir, mod.name)

        used_placeholder = false
        if source_export_abs !== nothing
            write_text(target_html_abs, read_text(source_export_abs))
        elseif !isfile(target_html_abs)
            html = placeholder_html(mod.name, mod.source, mod.description)
            write_text(target_html_abs, html)
            used_placeholder = true
        end

        # Si existia previamente un HTML (real o placeholder), se considera publicado.
        published = isfile(target_html_abs)
        rel_published_path = relative_to_root(target_html_abs, docs_dir) # Pluto/.../index.html
        push!(
            manifest_entries,
            build_manifest_entry(mod.name, rel_published_path, mod.source, published),
        )

        if published
            push!(published_modules, mod.name)
            if used_placeholder
                push!(placeholder_modules, mod.name)
            elseif source_export_abs === nothing
                # Archivo existente y no sobreescrito: estado desconocido, no se marca placeholder.
            end
        end
    end

    manifest_path = joinpath(docs_dir, "pages_manifest.json")
    manifest_content = "[\n" * join(manifest_entries, ",\n") * "\n]\n"
    write_text(manifest_path, manifest_content)

    println("\nResumen:")
    println(" - Manifest generado: ", manifest_path)
    println(" - Modulos con pagina disponible: ", length(published_modules), "/", length(MODULES))
    println(" - Lista publicados: ", isempty(published_modules) ? "-" : join(published_modules, ", "))
    println(
        " - Lista en placeholder (creados en esta ejecucion): ",
        isempty(placeholder_modules) ? "-" : join(placeholder_modules, ", "),
    )
end

main()
