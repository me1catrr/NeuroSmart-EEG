using Test
include("../src/preprocessing/artifact_rejection.jl")

@testset "Rechazo de artefactos" begin
    # Simulamos datos con un artefacto extremo
    epochs = rand(32, 256, 5) .* 20.0  # µV normal
    epochs[:, :, 3] .= 200.0  # Artefacto en epoch 3 (>70 µV)

    clean, removed, report = ArtifactRejection.reject_artifacts(epochs; threshold=70.0)

    @test length(removed) == 1
    @test report["rechazados"] == 1
    @test size(clean, 3) == 4
end