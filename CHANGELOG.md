# Changelog
Todos los cambios notables en este proyecto se documentarán en este archivo.

El formato sigue [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) y este proyecto se adhiere a [Semantic Versioning](https://semver.org/).

## [Unreleased]
- Implementación futura: métricas de conectividad funcional.
- Mejora del pipeline para análisis avanzado (EMD, CEEMD, PLV).

## [0.2.0] - 2025-07-31
### Added
- **Pipeline de análisis espectral**:
  - `segment_eeg` (segmentación en epochs).
  - `baseline_correct` (corrección de línea base).
  - `reject_artifacts` (rechazo de artefactos por amplitud).
  - `fft_analysis` (análisis FFT y cálculo de potencia por bandas).
- Documentación técnica:
  - `docs/api_reference.md` y `docs/methodology.md`.

### Changed
- Mejoras en `preprocess_eeg.jl`:
  - Filtrado Bandpass (0.5–45 Hz) y Notch 50 Hz.
  - Resampleo uniforme a 512 Hz.
  - Conversión segura a Float64.
  - Logs más limpios y detallados.

### Fixed
- Corrección en la lectura de `gain` en archivos `.vhdr`.
- Ajuste en reshaping tras filtrado y resampleo para evitar errores de dimensión.

## [0.1.0] - 2025-07-29
### Added
- Estructura inicial del proyecto:
  - Carpetas: `src/`, `scripts/`, `data/`, `tests/`.
  - Archivo principal `preprocess_eeg.jl`.
- Primera versión del pipeline de preprocesamiento:
  - Lectura de cabeceras `.vhdr`.
  - Carga de datos binarios `.eeg`.
  - Aplicación de filtros básicos.
  - Resampleo a 512 Hz.
- Inclusión de `README.md` con estructura y flujo de trabajo.