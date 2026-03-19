# EEG_Julia

Pipeline de análisis de señales EEG en **Julia** para el proyecto BRAIN (Investigación universitaria). Incluye carga de datos, preprocesamiento, filtrado, ICA, segmentación, corrección de baseline, rechazo de artefactos, análisis espectral (FFT) y medidas de conectividad (CSD, wPLI).

## Requisitos

- **Julia** ≥ 1.12 (recomendado; el `Manifest.toml` está generado con 1.12.3)
- Dependencias del proyecto (se instalan con `Pkg.instantiate()`)

## Instalación

1. Clona o abre el proyecto y entra en su directorio:

   ```bash
   cd EEG_Julia
   ```

2. Instala dependencias en el entorno del proyecto:

   - Opción A (una línea):

     ```bash
     julia --project=. -e 'using Pkg; Pkg.instantiate()'
     ```

   - Opción B (REPL):

     ```bash
     julia --project=.
     ]
     instantiate
     backspace
     ```

   - Opción C (desde un script Julia):

     ```julia
     using Pkg
     Pkg.activate(".")
     Pkg.instantiate()
     ```

## Estructura del proyecto

```text
EEG_Julia/
├── config/
│   └── default_config.jl     # Configuración tipada del pipeline (PipelineConfig)
├── Pluto/
│   ├── EEG_Pipeline.jl       # Notebook Pluto (pipeline completo interactivo)
│   └── Notebook.jl           # Notebook auxiliar
├── script/
│   └── EEG.jl                # Punto de entrada (pipeline completo, secuencial)
├── src/
│   ├── BIDS/
│   │   ├── build_participants.jl
│   │   ├── build_eeg_bids.jl
│   │   └── validate_bids.jl
│   ├── Preprocessing/
│   │   ├── IO.jl             # Carga inicial/preprocesamiento de datos EEG
│   │   └── filtering.jl      # Filtrado (notch, bandreject, highpass, lowpass)
│   ├── Processing/
│   │   ├── ICA.jl
│   │   ├── ICA_cleaning.jl
│   │   ├── segmentation.jl
│   │   ├── baseline.jl
│   │   ├── artifact_rejection.jl
│   │   ├── baseline_2st.jl
│   │   └── FFT.jl
│   ├── Connectivity/
│       ├── CSD.jl            # Current Source Density (spline esférico Perrin)
│       └── wPLI.jl           # Weighted Phase Lag Index por bandas
│   └── modules/
│       ├── EEG_Julia.jl      # Módulo principal: incluye etapas + exporta `run_*`
│       ├── paths.jl          # Utilidades de rutas (project_root, stage_dir, etc.)
│       └── utils.jl          # Utilidades comunes (logs, binarios, limpieza)
├── data/
│   ├── raw/                          # Entrada principal EEG (TSV/derivados iniciales)
│   ├── electrodes/                   # Posiciones de electrodos (TSV)
│   ├── BIDS/                         # Salidas de conversión/validación BIDS
│   ├── Preprocessing/
│   │   ├── IO/                       # Salida de carga inicial
│   │   └── filtering/                # Salidas de filtrado
│   ├── Processing/
│   │   ├── ICA/
│   │   ├── segmentation/
│   │   ├── baseline/
│   │   ├── artifact_rejection/
│   │   └── FFT/
│   └── Connectivity/
│       ├── CSD/
│       └── wPLI/
├── results/
│   ├── BIDS/
│   │   ├── figures/
│   │   ├── logs/
│   │   └── tables/
│   ├── Preprocessing/
│   │   ├── figures/
│   │   ├── logs/
│   │   └── tables/
│   ├── Processing/
│   │   ├── figures/
│   │   ├── logs/
│   │   └── tables/
│   └── Connectivity/
│       ├── figures/
│       ├── logs/
│       └── tables/
├── docs/
│   ├── index.html            # Sitio estático publicado en GitHub Pages
│   └── .nojekyll             # Evita procesamiento Jekyll en Pages
├── Project.toml
├── Manifest.toml
└── README.md
```

## Privacidad y datos (importante)

Por tratarse de un proyecto de **EEG/medicina**, este repositorio **NO incluye** datos ni derivados por defecto.

- **No se suben a Git**: `data/`, `results/`, `Javier_results/`
- **Qué sí se sube**: código (`src/`, `script/`), configuración (`config/`), notebooks (`Pluto/`) y documentación.

Para ejecutar el pipeline necesitarás proporcionar tus propios archivos de entrada en `data/raw/` (y, si aplica, `data/electrodes/`) en tu entorno local.

## Publicación en GitHub Pages

La carpeta `docs/` se mantiene porque se publica en GitHub Pages.

- `docs/index.html`: página estática a servir.
- `docs/.nojekyll`: evita que GitHub Pages aplique Jekyll sobre el contenido.
- En la configuración del repositorio, la fuente de Pages debe apuntar a la carpeta `docs/` de la rama publicada.

## Flujo del pipeline

El procesamiento sigue un orden secuencial; cada etapa lee la salida de la anterior. El módulo `src/modules/EEG_Julia.jl` expone funciones de alto nivel `run_*` para ejecutar etapas individualmente (REPL/Pluto) o en un flujo completo (vía `script/EEG.jl`).

Antes del pipeline principal, el proyecto incluye utilidades en `src/BIDS/` para convertir datos originales de BrainVision a una estructura BIDS (`build_participants.jl`, `build_eeg_bids.jl`, `validate_bids.jl`).

1. **BIDS (preparación opcional)**: Conversión/validación desde BrainVision hacia estructura BIDS.
2. **IO**: Carga TSV raw, organiza canales (diccionario), PSD, calidad de canales.
3. **filtering**: Notch 50 Hz, bandreject ~100 Hz, highpass 0.5 Hz, lowpass 150 Hz (Butterworth, filtfilt).
4. **ICA**: FastICA simétrico sobre datos filtrados; guarda componentes y matrices de mezcla.
5. **ICA_cleaning**: Evaluación automática de ICs (features/scores), eliminación de artefactos, reconstrucción.
6. **segmentation**: División en segmentos de longitud fija (ej. 1 s), sin solapamiento; hipermatriz 3D.
7. **baseline**: 1ª corrección de baseline por segmento (intervalo 0–0.1 s).
8. **artifact_rejection**: Rechazo de épocas por umbral de amplitud (±70 µV).
9. **baseline_2st**: 2ª corrección de baseline (mismo intervalo).
10. **FFT**: FFT, ventana Hamming, zero-padding; potencia espectral y por bandas (Delta–Gamma).
11. **CSD**: Current Source Density (Laplaciano esférico Perrin) sobre datos segmentados.
12. **wPLI**: Conectividad wPLI por bandas de frecuencia (entrada: datos CSD).

Los datos intermedios se guardan en `data/` (p. ej. `.bin` serializados); figuras y tablas en `results/figures/` y `results/tables/`.

## Configuración

La configuración central está en **`config/default_config.jl`** y es **tipada** mediante `struct`s (p. ej. `FilterConfig`, `ICAConfig`, etc.) agregadas en `PipelineConfig`.

- **Rutas base** (constantes): `ROOT_DIR`, `DATA_DIR`, `RESULTS_DIR`, `FIGURES_DIR`, `LOGS_DIR`, `TABLES_DIR`
- **`PipelineConfig`**:
  - `data_dir`, `output_dir`
  - `filter` (`fs`, `notch_freq`, `bandreject_center`, `bandreject_width`, `hp_cutoff`, `lp_cutoff`)
  - `ica` (`method`, `max_steps`, `tol`, `components_to_zero`)
  - `segmentation` (`length_s`, `overlap_s`)
  - `artifact` (`amp_min`, `amp_max`)
  - `fft` (`pad_to`, `window`, `full_spectrum`)

Además, `src/modules/paths.jl` proporciona utilidades para construir rutas sin hardcodear strings por todo el proyecto:

- `project_root()`, `data_root(cfg)`, `results_root(cfg)`
- `stage_dir(:IO; kind=:data, cfg)` -> `data/Preprocessing/IO`
- `stage_dir(:FFT; kind=:figures, cfg)` -> `results/Processing/figures/FFT`
- `raw_dir(cfg)`, `electrodes_dir(cfg)`

## Cómo ejecutar

- **Pipeline completo (recomendado)**  
  Desde la raíz del proyecto:

  ```bash
  julia script/EEG.jl
  ```

  Este script:
  - activa el entorno del proyecto (`Pkg.activate(...)`)
  - instancia dependencias (`Pkg.instantiate()`)
  - carga `config/default_config.jl`
  - carga el módulo `src/modules/EEG_Julia.jl`
  - ejecuta secuencialmente `run_io(cfg)`, `run_filtering(cfg)`, ..., `run_wpli(cfg)`

- **Ejecutar una etapa suelta (REPL / Pluto / `-e`)**  
  La forma robusta es cargar la config y el módulo, y luego llamar a la etapa:

  ```bash
  julia --project=. -e 'include("config/default_config.jl"); using .DefaultConfig; include("src/modules/EEG_Julia.jl"); using .EEG_Julia; run_io(DEFAULT_CONFIG)'
  ```

  Cambia `run_io` por `run_filtering`, `run_ica`, `run_fft`, etc.

## Datos de entrada

- **EEG raw**: se espera un archivo TSV en `data/raw/` (ej. `sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv`) con:
  - Primera columna: nombres de canales (`Channel`)
  - Resto de columnas: muestras temporales (una columna por punto de tiempo), en µV
- **Metadata** y **electrodos**: opcionales en `data/raw/` y `data/electrodes/` (p. ej. para CSD y visualizaciones).

Frecuencia de muestreo por defecto en el pipeline: **500 Hz**.

## Dependencias principales (Project.toml)

- **CSV**, **DataFrames**: lectura de TSV y tablas  
- **DSP**, **FFTW**: filtros y FFT  
- **Plots**, **GR**, **CairoMakie**, **StatsPlots**: gráficos  
- **MultivariateStats**, **Statistics**, **StatsBase**: estadísticas e ICA  
- **Serialization**, **Dates**: guardado binario y logs  

## Bandas de frecuencia (FFT / wPLI)

Definiciones típicas en el proyecto:

- **Delta (Δ)**: 0.5–4 Hz  
- **Theta (θ)**: 4–8 Hz  
- **Alpha (α)**: 8–12 Hz  
- **Beta low**: 12–15 Hz  
- **Beta mid**: 15–18 Hz  
- **Beta high**: 18–30 Hz  
- **Gamma (γ)**: 30–50 Hz  

## Resultados típicos

- **results/figures/**: PSD, topografías ICA, mapas CSD, heatmaps wPLI, etc.  
- **results/tables/**: potencia por banda (FFT), matrices y listas de conexiones (wPLI), estadísticas de segmentos, baseline y artefactos.  
- **results/logs/**: logs con fecha para CSD, wPLI y otros módulos.  

## Referencias metodológicas

- **Filtrado**: Butterworth, zero-phase (`filtfilt`).  
- **ICA**: FastICA simétrico (blanqueo PCA + maximización de no-gaussianidad).  
- **CSD**: Método de Perrin (spline esférico, Laplaciano), compatible con BrainVision Analyzer.  
- **wPLI**: Filtrado bandpass por banda, Hilbert, agregación sobre muestras y segmentos (estilo BrainVision).  
- **FFT**: Ventana Hamming, zero-padding (ej. 512), potencia en µV² y resolución ~0.98 Hz.  

## Licencia y autoría

Proyecto de investigación BRAIN (Universidad). Para uso interno y académico.
