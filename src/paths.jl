"""
Utilidades de rutas del proyecto.

Objetivo: evitar strings hardcodeadas dispersas y construir rutas de forma robusta
con `joinpath` y una única noción de "raíz del proyecto".
"""

"""
    project_root() -> String

Devuelve la raíz del proyecto (carpeta que contiene `src/`, `data/`, `results/`, etc.).
"""
project_root() = normpath(joinpath(@__DIR__, ".."))

"""`data_root(cfg=nothing)` devuelve la carpeta base de `data/`."""
function data_root(cfg=nothing)
    if cfg !== nothing && hasproperty(cfg, :data_dir)
        return normpath(getproperty(cfg, :data_dir))
    end
    return joinpath(project_root(), "data")
end

"""`results_root(cfg=nothing)` devuelve la carpeta base de `results/`."""
function results_root(cfg=nothing)
    if cfg !== nothing && hasproperty(cfg, :output_dir)
        return normpath(getproperty(cfg, :output_dir))
    end
    return joinpath(project_root(), "results")
end

"""
    stage_dir(stage::Symbol; kind=:data, cfg=nothing) -> String

Ruta a la carpeta de una etapa, bajo `data/` o bajo `results/`.

Ejemplos:
- `stage_dir(:IO)` → `data/IO`
- `stage_dir(:FFT; kind=:figures)` → `results/figures/FFT`
"""
function stage_dir(stage::Symbol; kind::Symbol = :data, cfg=nothing)
    st = String(stage)
    if kind === :data
        return joinpath(data_root(cfg), st)
    elseif kind === :results
        return joinpath(results_root(cfg), st)
    elseif kind === :figures
        return joinpath(results_root(cfg), "figures", st)
    elseif kind === :tables
        return joinpath(results_root(cfg), "tables", st)
    elseif kind === :logs
        return joinpath(results_root(cfg), "logs", st)
    else
        throw(ArgumentError("kind inválido: $kind (usa :data, :results, :figures, :tables, :logs)"))
    end
end

"""Ruta a `data/raw/`."""
raw_dir(cfg=nothing) = joinpath(data_root(cfg), "raw")

"""Ruta a `data/electrodes/`."""
electrodes_dir(cfg=nothing) = joinpath(data_root(cfg), "electrodes")

