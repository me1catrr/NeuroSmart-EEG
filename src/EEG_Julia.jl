module EEG_Julia

"""
Módulo principal del proyecto EEG_Julia.

- Centraliza utilidades reutilizables (rutas, IO binario, logging).
- Expone funciones de alto nivel `run_*` para ejecutar cada etapa de forma individual
  (útil desde REPL y notebooks Pluto.jl).
"""

# Utilidades comunes (deben ir primero)
include("paths.jl")
include("utils.jl")

# Etapas del pipeline
include("IO.jl")
include("filtering.jl")
include("ICA.jl")
include("ICA_cleaning.jl")
include("segmentation.jl")
include("baseline.jl")
include("artifact_rejection.jl")
include("baseline_2st.jl")
include("FFT.jl")

include(joinpath("Connectivity", "CSD.jl"))
include(joinpath("Connectivity", "wPLI.jl"))

# Exports: solo lo "alto nivel" + rutas/IO comunes útiles desde Pluto
export project_root, data_root, results_root, stage_dir, raw_dir, electrodes_dir
export ensure_dir, clean_dir_by_ext, run_id, logmsg, open_log, save_bin, load_bin, pick_first_existing

export run_io, run_filtering, run_ica, run_ica_cleaning, run_segmentation
export run_baseline, run_artifact_rejection, run_baseline_2st, run_fft
export run_csd, run_wpli

end # module EEG_Julia

