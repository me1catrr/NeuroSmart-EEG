"""
reject_artifacts(epochs::Array{Float64, 3}; threshold=70.0)

Rechaza epochs con amplitud excesiva (artefactos) según un criterio simple de umbral.

# Argumentos
- `epochs::Array{Float64, 3}`: EEG segmentado [canales x muestras x epochs].
- `threshold::Float64`: Umbral de amplitud en µV (por defecto 70 µV).

# Retorno
- `clean_epochs`: Array con los epochs aceptados.
- `removed_epochs`: Vector con índices de los epochs rechazados.
- `report`: Diccionario con estadísticas del proceso.

# Ejemplo
```julia
clean, removed, info = reject_artifacts(epochs; threshold=80.0)
println(info["porcentaje_rechazado"])
"""

module ArtifactRejection

export reject_artifacts

function reject_artifacts(epochs::Array{Float64, 3}; threshold=70.0)

n_epochs = size(epochs, 3)
keep = trues(n_epochs)
removed_epochs = Int[]

for i in 1:n_epochs
    if maximum(abs.(epochs[:, :, i])) > threshold
        keep[i] = false
        push!(removed_epochs, i)
    end
end

clean_epochs = epochs[:, :, keep]

report = Dict(
    "total_epochs" => n_epochs,
    "rechazados" => length(removed_epochs),
    "aceptados" => sum(keep),
    "porcentaje_rechazado" => round(100 * length(removed_epochs) / n_epochs, digits=2)
)

return clean_epochs, removed_epochs, report