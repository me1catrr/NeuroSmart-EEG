include("../src/preprocessing/segmentation.jl")
include("../src/preprocessing/artifact_rejection.jl")
include("../src/spectral/fft_analysis.jl")

using BSON, DataFrames, CSV, Statistics

raw_dir = "data/preprocessed"
features_dir = "data/features"

mkpath(features_dir)

for group in ["patients", "controls"]
    group_path = joinpath(raw_dir, group)
    for subject in readdir(group_path)
        subj_path = joinpath(group_path, subject)
        for session in readdir(subj_path)
            println("Procesando: $subject / $session")

            # Load preprocessed EEG
            file_path = joinpath(subj_path, session, "EEG_preprocessed.bson")
            data_dict = BSON.load(file_path)
            data = data_dict[:data_final]
            fs = data_dict[:fs_target]

            # Segmentación
            epochs = segment_eeg(data, fs; epoch_len=1.0)

            # Baseline correction
            for i in 1:size(epochs, 3)
                epochs[:, :, i] = baseline_correct(epochs[:, :, i], fs)
            end

            # Rechazo de artefactos
            clean_epochs, removed, report = reject_artifacts(epochs; threshold=70.0)
            println("   Rechazados: $(report["rechazados"]) / $(report["total_epochs"])")

            # FFT y band power
            df = compute_band_powers(clean_epochs, fs)

            # Guardar
            save_dir = joinpath(features_dir, group, subject, session)
            mkpath(save_dir)
            CSV.write(joinpath(save_dir, "fft_features.csv"), df)

            # Log
            open(joinpath(save_dir, "spectral_log.txt"), "w") do log
                println(log, "=== Análisis Espectral ===")
                println(log, "Sujeto: $subject | Sesión: $session")
                println(log, "Epochs totales: $(report["total_epochs"])")
                println(log, "Epochs rechazados: $(report["rechazados"])")
                println(log, "Guardado en: fft_features.csv")
            end
        end
    end
end