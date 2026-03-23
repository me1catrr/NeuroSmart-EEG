#!/usr/bin/env julia

# -----------------------------------------------------------------------------
# tools/pages/publish.jl
#
# Flujo automatizado para publicar cambios de GitHub Pages:
#   1) Ejecuta tools/pages/build.jl
#   2) Muestra archivos modificados
#   3) git add .
#   4) git commit (mensaje por defecto o personalizado)
#   5) git push
#
# Uso:
#   julia tools/pages/publish.jl
#   julia tools/pages/publish.jl "Mensaje personalizado"
# -----------------------------------------------------------------------------

const DEFAULT_COMMIT_MESSAGE = "Update GitHub Pages (auto build)"

function project_root()::String
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function run_capture(cmd::Cmd; cwd::String = project_root())
    cmd_in_dir = Cmd(cmd; dir = cwd)
    output_buffer = IOBuffer()
    process = run(
        pipeline(ignorestatus(cmd_in_dir), stdout = output_buffer, stderr = output_buffer),
        wait = true,
    )
    output = String(take!(output_buffer))
    return success(process), output
end

function run_or_error(step_name::String, cmd::Cmd; cwd::String = project_root())
    println("→ ", step_name)
    ok, output = run_capture(cmd; cwd = cwd)
    if !isempty(strip(output))
        println(output)
    end
    if !ok
        error("Fallo en paso: $(step_name)")
    end
end

function get_changed_files()::Vector{String}
    ok, output = run_capture(`git status --short`)
    if !ok
        error("No se pudo consultar el estado de git. Verifica que estas en un repositorio valido.")
    end

    lines = split(chomp(output), '\n')
    cleaned = String[]
    for line in lines
        striped = strip(line)
        if !isempty(striped)
            push!(cleaned, striped)
        end
    end
    return cleaned
end

function has_staged_changes()::Bool
    ok, output = run_capture(`git diff --cached --name-only`)
    if !ok
        error("No se pudo verificar el area de staging.")
    end
    return !isempty(strip(output))
end

function main()
    root = project_root()
    commit_message = isempty(ARGS) ? DEFAULT_COMMIT_MESSAGE : join(ARGS, " ")

    println("==============================================")
    println(" NeuroSmart-EEG · Publicacion GitHub Pages")
    println("==============================================")
    println("Repositorio: ", root)
    println("Mensaje de commit: \"", commit_message, "\"")
    println("")

    # Paso 1: build de Pages
    build_script = joinpath(root, "tools", "pages", "build.jl")
    if !isfile(build_script)
        error("No se encontro tools/pages/build.jl")
    end
    run_or_error("Construyendo Pages", `julia $build_script`)

    # Paso 2: mostrar cambios actuales
    println("→ Revisando archivos modificados")
    changed_before_add = get_changed_files()
    if isempty(changed_before_add)
        println("No hay cambios para publicar. Fin del proceso.")
        return
    end
    println("Archivos modificados detectados:")
    for line in changed_before_add
        println("  - ", line)
    end
    println("")

    # Paso 3: add
    run_or_error("Agregando cambios (git add .)", `git add .`)

    # Paso 4: commit (solo si hay cambios staged)
    if !has_staged_changes()
        println("No hay cambios en staging despues de git add. Fin del proceso.")
        return
    end
    run_or_error("Creando commit", `git commit -m $commit_message`)

    # Paso 5: push
    run_or_error("Haciendo push a remoto", `git push`)

    println("✅ Publicacion completada.")
    println("Siguiente verificacion sugerida: revisar el sitio de GitHub Pages en 1-3 minutos.")
end

try
    main()
catch err
    println("❌ Error durante publish: ", err)
    println("Proceso detenido. Revisa el mensaje anterior para corregir el problema.")
    rethrow(err)
end
