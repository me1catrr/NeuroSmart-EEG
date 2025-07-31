using Test, FFTW
include("../src/spectral/fft_analysis.jl")

@testset "Análisis FFT" begin
    fs = 256.0
    t = 0:1/fs:1-1/fs
    signal = sin.(2π*10 .* t)  # 10 Hz
    data = repeat(reshape(signal, 1, :), 32, 1)  # 32 canales iguales

    freqs, psd = FFTAnalysis.compute_fft(data, fs)
    peak_freq = freqs[argmax(psd[:, 1])]

    @test abs(peak_freq - 10) < 0.5
end