using Test

# Importar los tests de cada módulo
include("test_segmentation.jl")
include("test_artifact_rejection.jl")
include("test_fft_analysis.jl")

println("✅ Todos los tests fueron ejecutados correctamente.")