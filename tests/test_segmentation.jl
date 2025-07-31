using Test
include("../src/preprocessing/segmentation.jl")

@testset "Segmentación EEG" begin
    fs = 256.0
    data = rand(32, 2560)  # 32 canales, 10 segundos
    epochs = Segmentation.segment_eeg(data, fs; epoch_len=1.0)

    @test size(epochs, 3) == 10  # 10 epochs
    @test size(epochs)[1:2] == (32, 256)  # Canales x Muestras por epoch
end