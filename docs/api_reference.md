# API Reference - NeuroSmart-EEG

## Módulos principales:

### 1. Preprocessing
- **preprocess_eeg.jl**
  - `read_vhdr(path)`: Lee cabecera .vhdr.
  - `read_eeg(path, n_channels)`: Carga datos EEG.
  - `filter_eeg(data, fs)`: Filtrado bandpass y notch.
  - `resample_eeg(data, fs_original, fs_target)`: Resampleo.

### 2. Segmentation
- **segmentation.jl**
  - `segment_eeg(data, fs; epoch_len)`: Divide en epochs.
  - `baseline_correct(epoch, fs; baseline)`: Corrección de línea base.

### 3. Artifact Rejection
- **artifact_rejection.jl**
  - `reject_artifacts(epochs; threshold)`: Elimina epochs con amplitud excesiva.

### 4. Spectral Analysis
- **fft_analysis.jl**
  - `compute_fft(epochs, fs)`: Calcula espectro.
  - `band_power(fft_data, bands)`: Potencia por bandas.

---

## Pipeline Scripts
- `run_preprocessing.jl`: Ejecuta preprocesamiento.
- `run_spectral.jl`: Segmenta, limpia y calcula FFT.

---

### Salida esperada:
- Preprocesado: `data/preprocessed/...`
- Features: `data/features/...`