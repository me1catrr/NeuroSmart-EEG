module Segmentation

using Statistics

export segment_eeg, baseline_correct

"""
    segment_eeg(data::Matrix{Float64}, fs::Float64; epoch_len=1.0)

Divide la señal EEG continua en segmentos (*epochs*) de longitud fija.

# Argumentos
- `data::Matrix{Float64}`: Matriz EEG [canales x muestras].
- `fs::Float64`: Frecuencia de muestreo en Hz.
- `epoch_len::Float64`: Duración de cada epoch en segundos (default = 1.0).

# Retorno
- Array de 3D: [canales x muestras_por_epoch x n_epochs].

# Notas
- Si el número total de muestras no es múltiplo exacto del tamaño del epoch, se descartan las muestras sobrantes.
"""
function segment_eeg(data::Matrix{Float64}, fs::Float64; epoch_len=1.0)
    samples_per_epoch = Int(round(epoch_len * fs))
    n_channels, n_samples = size(data)
    n_epochs = div(n_samples, samples_per_epoch)

    if n_epochs == 0
        error("❌ La señal es demasiado corta para segmentar en epochs de $(epoch_len)s.")
    end

    epochs = Array{Float64}(undef, n_channels, samples_per_epoch, n_epochs)
    for i in 1:n_epochs
        idx = ((i - 1) * samples_per_epoch + 1):(i * samples_per_epoch)
        epochs[:, :, i] = data[:, idx]
    end

    return epochs
end

"""
    baseline_correct(epoch::Matrix{Float64}, fs::Float64; baseline=(0.0, 0.1))

Corrige el baseline de un epoch restando la media en un intervalo inicial.

# Argumentos
- `epoch::Matrix{Float64}`: Señal de un epoch [canales x muestras].
- `fs::Float64`: Frecuencia de muestreo (Hz).
- `baseline::Tuple{Float64,Float64}`: Ventana para baseline en segundos (default = (0.0, 0.1)).

# Retorno
- Epoch corregido.

# Notas
- Si el baseline excede la longitud del epoch, se limita al rango disponible.
"""
function baseline_correct(epoch::Matrix{Float64}, fs::Float64; baseline=(0.0, 0.1))
    n_samples = size(epoch, 2)
    idx_start = clamp(Int(round(baseline[1] * fs)) + 1, 1, n_samples)
    idx_end   = clamp(Int(round(baseline[2] * fs)), 1, n_samples)

    correction = Statistics.mean(epoch[:, idx_start:idx_end], dims=2)
    return epoch .- correction
end

end # module