module PlutoTemplateBase

using Dates
using Markdown
using CSV
using DataFrames
using Serialization
using Statistics
using DSP
using Plots

# Rutas del proyecto centralizadas en src/modules/paths.jl
include(joinpath(@__DIR__, "..", "src", "modules", "paths.jl"))

const FS_DEFAULT = 500.0
const NOTEBOOK_AUTHOR = "RAFAEL CASTRO TRIGUERO, ALEJANDRO GALVAO Y JAVIER ESPUNY"

function published_date_str(d::Date = today())
    month_names_es = [
        "enero", "febrero", "marzo", "abril", "mayo", "junio",
        "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre",
    ]
    return "$(Dates.day(d)) de $(month_names_es[Dates.month(d)]) de $(Dates.year(d))"
end

function notebook_intro(title::AbstractString; author::AbstractString = NOTEBOOK_AUTHOR, published::Date = today())
    return Markdown.parse("""
# $(title)

| AUTOR | PUBLICADO |
|---|---|
| $(author) | $(published_date_str(published)) |
""")
end

export FS_DEFAULT
export NOTEBOOK_AUTHOR, published_date_str, notebook_intro
export project_root, data_root, results_root, stage_dir, bids_root, raw_dir, electrodes_dir

end # module PlutoTemplateBase
