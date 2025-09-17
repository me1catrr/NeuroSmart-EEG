module FFTAnalysis

using FFTW
using Statistics
using DataFrames

export compute_band_powers

"""
    compute_band_powers(epochs::Array{Float64, 3}, fs::Float64)

Calcula potencia en bandas (delta, theta, alpha, beta, gamma) para cada canal y epoch.

# Argumentos
- `epochs`: Array [canales × muestras × epochs]
- `fs`: Frecuencia de muestreo (Hz)

# Retorno
- `df`: DataFrame con columnas:
    :epoch, :channel, :delta, :theta, :alpha, :beta, :gamma
"""
function compute_band_powers(epochs::Array{Float64, 3}, fs::Float64)
    n_channels, n_samples, n_epochs = size(epochs)
    freqs = collect(0:(fs/n_samples):(fs/2))
    bands = Dict(
        "delta" => (0.5, 4),
        "theta" => (4, 8),
        "alpha" => (8, 12),
        "beta"  => (12, 30),
        "gamma" => (30, 45)
    )

    rows = []

    for e in 1:n_epochs
        for ch in 1:n_channels
            signal = epochs[ch, :, e]
            # FFT
            fft_vals = abs.(fft(signal))[1:div(n_samples,2)+1]
            psd = (fft_vals .^ 2) / (fs * n_samples)

            # Band powers
            band_powers = Dict()
            for (band, (fmin, fmax)) in bands
                idx = findall(f -> f ≥ fmin && f ≤ fmax, freqs)
                band_powers[band] = Statistics.mean(psd[idx])
            end

            push!(rows, (
                epoch = e,
                channel = ch,
                delta = band_powers["delta"],
                theta = band_powers["theta"],
                alpha = band_powers["alpha"],
                beta = band_powers["beta"],
                gamma = band_powers["gamma"]
            ))
        end
    end

    return DataFrame(rows)
end

end # module