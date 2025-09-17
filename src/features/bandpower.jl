"""
compute_band_power(power, freqs, bands)

Calcula la potencia promedio por banda EEG.

# Argumentos
- `power::Matrix`: Potencia espectral [canales x frecuencias].
- `freqs::Vector`: Frecuencias asociadas.
- `bands::Dict`: Diccionario con bandas, e.g., Dict("alpha" => (8,13)).

# Retorna
- Dict con potencia media por banda y canal.

# Uso
bands = Dict("delta" => (0.5,4), "theta" => (4,8), "alpha" => (8,13))
bandpower = compute_band_power(power, freqs, bands)
"""

module BandPower

using Statistics

export compute_band_power

function compute_band_power(power::Matrix{Float64}, freqs::Vector{Float64}, bands::Dict)
    band_powers = Dict()
    for (band, range) in bands
        idx = findall(f -> f ≥ range[1] && f ≤ range[2], freqs)
        band_powers[band] = mean(power[:, idx], dims=2)
    end
    return band_powers
end

end