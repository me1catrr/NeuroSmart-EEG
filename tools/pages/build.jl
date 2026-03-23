#!/usr/bin/env julia

# Construye contenido de GitHub Pages desde exportaciones HTML de Pluto.
#
# Flujo:
#   - Fuente cientifica: Pluto/<Modulo>/<Modulo>.jl
#   - Staging HTML manual: exports/Pluto/<Modulo>/index.html
#   - Publicacion web: docs/Pluto/<Modulo>/index.html
#
# Uso:
#   julia tools/pages/build.jl
#   julia tools/pages/build.jl /ruta/custom/exports/Pluto

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

function write_text(path::String, content::String)::Nothing
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return nothing
end

function json_escape(value::String)::String
    escaped = replace(value, "\\" => "\\\\")
    escaped = replace(escaped, "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n")
    escaped = replace(escaped, "\r" => "\\r")
    escaped = replace(escaped, "\t" => "\\t")
    return escaped
end

function placeholder_html(module_name::String, source_path::String)::String
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
      <p><strong>Aun no hay exportacion HTML disponible para este modulo.</strong></p>
      <p>
        Exporta el notebook desde Pluto.jl y coloca el archivo en:
      </p>
      <div class="source">exports/Pluto/$(module_name)/index.html</div>
      <p>Notebook fuente:</p>
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

function default_exports_dir(root::String)::String
    return joinpath(root, "exports", "Pluto")
end

function staged_export_path(exports_pluto_dir::String, module_name::String)::String
    return joinpath(exports_pluto_dir, module_name, "index.html")
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

function copy_file(src::String, dst::String)::Nothing
    mkpath(dirname(dst))
    cp(src, dst; force = true)
    return nothing
end

function main()
    root = project_root()
    docs_dir = joinpath(root, "docs")
    exports_pluto_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : default_exports_dir(root)

    if !isdir(exports_pluto_dir)
        println("Aviso: no existe directorio de staging HTML: ", exports_pluto_dir)
        println("Se generaran placeholders para modulos sin exportacion.")
    end

    println("== NeuroSmart-EEG Pages Build ==")
    println("Raiz del proyecto: ", root)
    println("Staging de exportaciones Pluto: ", exports_pluto_dir)
    println("Destino GitHub Pages: ", joinpath(docs_dir, "Pluto"))

    manifest_entries = String[]
    real_modules = String[]
    placeholder_modules = String[]

    for mod in MODULES
        target_html_abs = joinpath(root, mod.target)
        staged_html_abs = staged_export_path(exports_pluto_dir, mod.name)

        published_real = false
        if isfile(staged_html_abs)
            copy_file(staged_html_abs, target_html_abs)
            published_real = true
        else
            html = placeholder_html(mod.name, mod.source)
            write_text(target_html_abs, html)
        end

        rel_published_path = relative_to_root(target_html_abs, docs_dir)
        push!(
            manifest_entries,
            build_manifest_entry(mod.name, rel_published_path, mod.source, true),
        )

        if published_real
            push!(real_modules, mod.name)
        else
            push!(placeholder_modules, mod.name)
        end
    end

    manifest_path = joinpath(docs_dir, "pages_manifest.json")
    manifest_content = "[\n" * join(manifest_entries, ",\n") * "\n]\n"
    write_text(manifest_path, manifest_content)

    println("\nResumen:")
    println(" - Manifest: ", manifest_path)
    println(" - Modulos con HTML real: ", isempty(real_modules) ? "-" : join(real_modules, ", "))
    println(" - Modulos en placeholder: ", isempty(placeholder_modules) ? "-" : join(placeholder_modules, ", "))
end

main()
