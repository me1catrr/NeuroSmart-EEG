### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 60620496-6fc4-45df-a2fc-eca64e76a52f
begin
	using PyMNE
	using Plots
	using CSV
	using DataFrames
	using Serialization
	using Dates 
	using DSP
	using Serialization
	using Statistics
	using StatsBase
	using StatsPlots
	using FFTW
	using UnfoldMakie;
	using CairoMakie
	using DataFramesMeta
	using UnfoldSim
	using Unfold
	using MakieThemes
	using TopoPlots;
end

# ╔═╡ 4dc6ee64-1d24-11f1-b5d2-1b41a361acaa
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido")
end

# ╔═╡ f36b8ac2-d883-4d32-a7dc-2072e987165d
#=
begin
# 1) Montaje base
montage = PyMNE.channels.make_standard_montage("standard_1020")

# 2) Labels y posiciones 3D
labels_all = pyconvert(Vector{String}, montage.ch_names)
ch_pos = pyconvert(Dict{String,Any}, montage.get_positions()["ch_pos"])

# 3) Elige tus 32 canales
wanted = [
    "Fp1","Fp2",
    "F7","F3","Fz","F4","F8",
    "FC5","FC1","FC2","FC6",
    "T7","C3","Cz","C4","T8",
    "CP5","CP1","CP2","CP6",
    "P7","P3","Pz","P4","P8",
    "PO9","O1","Oz","O2","PO10","TP9","TP10"
]

labels = [l for l in labels_all if l in wanted]

# 4) Posiciones 3D -> 2D
pos3d = hcat([ch_pos[l] for l in labels]...)
positions = to_positions(pos3d)

# 5) Tus datos: un valor por canal
# reemplaza esto con tus betas/ERP/voltajes/etc.
values = randn(length(labels))

# 6) Topoplot estilo UnfoldMakie
plot_topoplot(
    values;
    labels = labels,
    positions = Point2f.(positions),
    visual = (; label_text = true, label_scatter = false),
    axis = (; xlabel = "")
)
end
=#

# ╔═╡ a73289ad-5d2c-4452-b19e-8a0ef3e42836
#=
begin
	x = 0:0.01:10
	y = sin.(x)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel="x", ylabel="sin(x)")

	lines!(ax, x, y)

	fig
end
=#


# ╔═╡ 078ecfa6-0ca4-42d3-8ffb-867e45c65c1d
md"
### Pipeline fase 1 (Preproceso)

Después del filtrado, las señales se almacenan tras cada etapa (por ejemplo,   
`dict_EEG_Notch.bin`, `dict_EEG_Lowpass.bin`, etc...) para facilitar la **trazabilidad**.

Procesos adicionales como:

- **re-referenciación** (por ejemplo referencia promedio),
- **rechazo de artefactos** (por ejemplo mediante ICA),
- **segmentación en épocas de longitud fija**

pueden aplicarse en etapas posteriores del pipeline antes de la **estimación de conectividad**.
"

# ╔═╡ 1e3c7d2a-3f68-4107-81f9-8ef9f2245802
md"""
!!! info "Pipeline Fase 1: Preprocesamiento"

    **Objetivo:**
    - Eliminar interferencias de red eléctrica y del hardware, deriva lenta y ruido 	fuera de banda del **EEG crudo** generado en la **Fase 0**.
    - Almacenar las señales filtradas para análisis posteriores específicos por 		banda.

    **Entrada.**
    - **EEG** serializado procedente de la **Fase** ``0`` (por ejemplo,  			`data/IO/dict_EEG.bin`).
    - Frecuencia de muestreo ``f_s = 500`` Hz.

    **Salida.**
    - Señales filtradas almacenadas en `data/filtering/` tras cada etapa (Notch, 		Bandreject, Highpass, Lowpass).
    - Gráficas de respuesta del **filtro** (magnitud y fase) y comparaciones de 		**PSD promedio** (antes/después de cada etapa) para control de calidad.

    **Pasos de procesamiento (alineados con `src/filtering.jl`).**

    1. **Carga de datos EEG.**  
       Cargar el diccionario EEG desde `data/IO/`; definir ``f_s`` y la lista de canales.

    2. **Filtro Notch 50 Hz.**  
       Aplicar un filtro Butterworth *bandstop* (orden 4, ancho 1 Hz).  
       Guardar `dict_EEG_Notch.bin`.  
       Representar la respuesta del filtro y comparar la PSD promedio (señal original vs señal con notch).

    3. **Rechazo de banda 100 Hz.**  
       Aplicar un filtro Butterworth *bandstop* (orden 4, ancho de banda 1 Hz) sobre la salida del notch.  
       Guardar `dict_EEG_Bandreject.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    4. **Filtro paso alto 0.5 Hz.**  
       Aplicar un filtro Butterworth *high-pass* (orden 4, `filtfilt`).  
       Guardar `dict_EEG_Highpass.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    5. **Filtro paso bajo 150 Hz.**  
       Aplicar un filtro Butterworth *low-pass* (orden 4, `filtfilt`).  
       Guardar `dict_EEG_Lowpass.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    6. **Salida del preprocesamiento.**  
       Devolver el diccionario completamente filtrado (o la ruta a `dict_EEG_Lowpass.bin`) para la **Fase 2a (ICA)**.  

	Opcionalmente, el re-referenciado, la eliminación de artefactos y la segmentación en épocas pueden realizarse en un paso separado antes del cálculo de conectividad.
"""

# ╔═╡ bc054369-e61b-4a71-92a1-a5d08ff69224
md"
### Pipeline fase 2a (ICA decompostion)
"

# ╔═╡ 3cacaf2c-57df-4eba-a29e-f5b82feeb79c
md"""
!!! info "Pipeline Fase 2a: Descomposición ICA"

    **Objetivo.**
    - Descomponer el **EEG filtrado** en componentes estadísticamente independientes 	utilizando **FastICA simétrico**.
    - Generar las matrices ``S``, ``W_{\text{total}}`` y ``A`` para su posterior 		evaluación y limpieza de artefactos.

    **Entrada.**
    - **EEG filtrado** procedente de la Fase 1: `data/filtering/dict_EEG_Lowpass.bin`
      (diccionario indexado por canal, ``C \times N`` muestras).
    - **Frecuencia de muestreo** ``f_s = 500`` Hz.
    - **Canales** ``C = 32`` 
	- **Número de componentes** ``k = C``.

    **Salida:** (`data/ICA/dict_EEG_ICA.bin`)
    - ``S`` (ICs × muestras), 
	- ``W_{\text{total}}``
	- ``A``, lista de canales
	- `max_iter`, `tol`.

    **Pasos de procesamiento (alineados con `src/ICA.jl`).**

    1. **Carga del EEG filtrado.**  
       - Cargar el **EEG** desde `data/filtering/`
       - construir la matriz ``X \in \mathbb{R}^{C \times N}`` (canales × muestras) 		en un orden fijo de canales (por ejemplo, orden alfabético).

    2. **Centrado de la señal.**  
       Restar la media de cada canal en ``X``.

    3. **Blanqueamiento (*whitening*).**  
       - Calcular la matriz de covarianza de ``X``.  
       - Realizar descomposición en autovalores mediante **PCA** y construir la 			señal blanqueada ``Z`` con ``k = C`` componentes.

    4. **Algoritmo FastICA.**  
       - Inicializar ``W``
	   - Aplicar iteraciones hasta convergencia:
	       - Proyectar ``Y = WZ``.  
	       - Aplicar la función de contraste `tanh`  
	       - Actualizar ``W`` con decorrelación simétrica  

       El proceso se repite hasta que

       ```math
       \max_i \left|1 - |M_{ii}|\right| < \text{tol} = 10^{-7}
       ```

       o hasta alcanzar el número máximo de iteraciones ``\mathrm{max_iter} = 512``

    5. **Resultados.**
       Calcular:

       ```math
       S = WZ
       ```

       ```math
       W_{\text{total}} = W \cdot V_{\text{whit}}
       ```

       ```math
       A = W_{\text{total}}^{-1}
       ```

       Serializar los resultados en `data/ICA/dict_EEG_ICA.bin`.
"""



# ╔═╡ 89849536-05c6-41a9-bb94-7eabe4057549
md"
## Segmentación

After ICA cleaning (Phase 2b), the continuous EEG is free of the main artefactual sources and is ready to be divided into fixed-length epochs for spectral and connectivity analysis.

Segmentation produces a 3D hypermatrix

```math
\mathbf{E} \in \mathbb{R}^{C \times L \times K}
```

channels × samples per segment × segments.

Each segment is a contiguous block of

```math
L = F_s \cdot T_{\mathrm{seg}}
```

samples (e.g. $T_{\mathrm{seg}} = 1\,\text{s}$ at $F_s = 500\,\text{Hz}$ gives $L = 500$).

The step between consecutive segment onsets is

```math
\text{step} = L - \text{overlap}
```

With zero overlap, $\text{step} = L$ and segments are statistically independent, which simplifies FFT-based spectral estimation and later averaging.

The number of full segments is

```math
K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
```

where $N_{\mathrm{total}}$ is the length of the continuous signal; any trailing samples that do not form a complete segment are discarded.

The 1 s segment length balances spectral resolution ($\approx F_s/L \approx 1\,\text{Hz}$, adequate for band-based connectivity) with the assumption of weak stationarity within each segment, which is required for stable cross-spectral and wPLI estimation.

The implementation is in `src/segmentation.jl`, which loads `data/ICA/dict_EEG_ICA_clean.bin`, applies the splitting, computes per-segment statistics (mean, std, min, max, RMS), generates example plots (overlaid and stacked segments), and saves the hypermatrix and metadata to `data/segmentation/`.

Algorithm 1 and Table 1 summarise the step; the box below formalises Pipeline Phase 3.
"

# ╔═╡ 4c757bb9-70a7-4fc2-b9c3-d17a019ba75e
md"""
!!! info "Algoritmo Segmentación en épocas de longitud fija (`segmentation.jl`)"

    **Entrada:**  
    `dict_EEG` (diccionario de canales con EEG limpio tras ICA),  
    `F_s`, `T_seg`, `overlap`

    1. Cálculo de parámetros de segmentación

       ```text
       L      ← floor(F_s · T_seg)
       step   ← L − overlap
       N_total ← length(dict_EEG[channels[1]])
       ```

    2. Número de segmentos

       ```text
       K ← floor((N_total − L) / step) + 1
       ```

    3. Inicialización de la hipermatriz

       ```text
       Allocate E ∈ ℝ^{C × L × K}
       ```

    4. Segmentación por canal

       ```text
       for c = 1 to C

           signal ← dict_EEG[channels[c]]

           for k = 1 to K

               start ← (k − 1) · step + 1
               end   ← start + L − 1

               if end ≤ N_total
                   E[c,:,k] ← signal[start : end]
               else
                   copy available samples
                   pad or truncate if needed
               end

           end

       end
       ```

    **Salida:**  
    `E` (hipermatriz de segmentos),  
    `K` (número de segmentos),  
    `L` (muestras por segmento),  
    `step` (desplazamiento entre segmentos)

    Opcionalmente: estadísticas por segmento  
    `(mean, std, min, max, RMS)`
"""

# ╔═╡ 5839954c-87c2-4333-81a3-9a5f64021000
md"
### Pipeline fase 3 (Segmentation)
"

# ╔═╡ 9ea266c8-0ffb-49c9-94f8-b77d206b087f
md"""
!!! info "Pipeline Fase 3: Segmentación"

    **Objetivo.**
    - Dividir el EEG continuo limpio tras ICA en segmentos de longitud fija sin solapamiento.
    - Construir una hipermatriz tridimensional (canales × muestras × segmentos).
    - Generar metadatos y estadísticas necesarias para las fases posteriores de corrección de baseline y estimación de conectividad.

    **Entrada.**
    - EEG limpio tras ICA: `data/ICA/dict_EEG_ICA_clean.bin`  
      (señales continuas indexadas por canal).

    - Parámetros de segmentación:

      ```math
      F_s = 500\,\text{Hz}
      ```

      ```math
      T_{\mathrm{seg}} = 1\,\text{s}
      ```

      ```math
      \text{overlap} = 0\,\text{s}
      ```

    **Salida.**
    - `data/segmentation/eeg_segmented.bin`: hipermatriz

      ```math
      \mathbf{E} \in \mathbb{R}^{C \times L \times K}
      ```

      donde

      ```math
      L = 500
      ```

      y

      ```math
      K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
      ```

    - `data/segmentation/dict_segmentation_info.bin`: diccionario con:
        - `eeg_segmented`
        - `channels`
        - `fs`
        - `segment_length_s`
        - `segment_length_samples`
        - `segment_overlap_s`
        - `segment_step_samples`
        - `n_segments`
        - `n_channels`
        - `n_samples_total`

    - `data/segmentation/segment_statistics.csv`: estadísticas por segmento  
      (media, desviación estándar, mínimo, máximo y RMS) para verificación.

    - Opcional: gráficos de control de calidad  
      (primer canal con segmentos superpuestos o segmentos apilados).

    **Pasos de procesamiento (alineados con `src/segmentation.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_EEG_ICA_clean.bin` y obtener:

       - `channels`
       - número total de muestras

       ```math
       N_{\mathrm{total}}
       ```

    2. **Definición de parámetros de segmentación.**  
       Calcular la longitud del segmento:

       ```math
       L = F_s \cdot T_{\mathrm{seg}}
       ```

       Definir el paso entre segmentos (sin solapamiento):

       ```math
       \text{step} = L
       ```

       Calcular el número total de segmentos:

       ```math
       K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
       ```

    3. **Construcción de la hipermatriz.**  
       Reservar memoria para

       ```math
       \mathbf{E} \in \mathbb{R}^{C \times L \times K}
       ```

       Para cada canal \(c\) y segmento \(k\), copiar el bloque correspondiente de señal:

       ```math
       \mathbf{E}[c,:,k] =
       \text{signal}[(k-1)\text{step}+1 : (k-1)\text{step}+L]
       ```

    4. **Cálculo de estadísticas.**  
       Para cada segmento calcular:

       - media  
       - desviación estándar  
       - mínimo  
       - máximo  
       - RMS  

       Guardar los resultados en `segment_statistics.csv`.

    5. **Visualización (opcional).**  
       Generar gráficos de control de calidad, por ejemplo:

       - primeros segmentos superpuestos del primer canal  
       - todos los segmentos apilados  

    6. **Guardado de resultados.**  
       - Serializar la hipermatriz en `eeg_segmented.bin`.
       - Construir el diccionario `dict_segmentation_info` con los metadatos y guardarlo en `dict_segmentation_info.bin`.

       Los datos resultantes se utilizan como entrada para la **Fase 4: corrección de baseline**.
"""

# ╔═╡ 2ba2934f-e388-462e-a632-8bc906cc6a4c
md"""
!!! info "Algoritmo Corrección de baseline por segmento (`baseline.jl`)"

    **Entrada:**  
    ``E ∈ ℝ^{C × L × K}`` (EEG segmentado),  
    `F_s`, `t_start`, `t_end`

    1. Cálculo de índices del intervalo de baseline

       ```text
       n_start ← floor(t_start · F_s) + 1
       n_end   ← floor(t_end   · F_s)
       ```

    2. Inicialización de la señal corregida

       ```text
       E_corr ← copy(E)
       ```

    3. Corrección por canal y segmento

       ```text
       for c = 1 to C
           for k = 1 to K

               segment ← E_corr[c,:,k]

               baseline_interval ← segment[n_start : n_end]

               μ ← mean(baseline_interval)

               E_corr[c,:,k] ← segment − μ

           end
       end
       ```

    **Salida:**  
    `E_corr` (EEG segmentado con corrección de baseline aplicada)
"""

# ╔═╡ 4450bc71-a1f1-4f8b-b390-99a8221b9b7d
md"""
### Pipeline 4 (Baseline before artifacts)
"""

# ╔═╡ 4f237b4c-6598-4fcc-8dd3-2b9b07e141bd
md"""
!!! info "Pipeline Fase 4: Corrección de baseline"

    **Objetivo.**
    - Eliminar el desplazamiento DC de cada segmento restando la media de un intervalo de baseline.
    - Utilizar como baseline los primeros 100 ms de cada segmento.
    - Generar EEG segmentado corregido y estadísticas de verificación para los análisis posteriores de conectividad y espectro.

    **Entrada.**
    - EEG segmentado: `data/segmentation/dict_segmentation_info.bin`, que contiene:
        - `eeg_segmented` (hipermatriz \(C \times L \times K\))
        - `channels`
        - `fs`
        - `segment_length_samples`
        - `n_segments`
        - `n_channels`

    - Ventana de baseline:

      ```math
      t_{\mathrm{start}} = 0\,\text{s}, \quad t_{\mathrm{end}} = 0.1\,\text{s}
      ```

      (primeros **100 ms** de cada segmento).

    **Salida.**
    - `data/baseline/eeg_1st_baseline_correction.bin`: hipermatriz corregida

      ```math
      \mathbf{E}^{\mathrm{corr}}
      ```

      (canales × muestras × segmentos).

    - `data/baseline/dict_1st_baseline_correction.bin`: diccionario con:
        - `eeg_baseline_corrected`
        - `eeg_segmented_original`
        - `channels`
        - `fs`
        - parámetros de segmentación y baseline
        - `baseline_type = "segment-based"`
        - `correction_method = "subtract_baseline_mean"`

    - `data/baseline/baseline_correction_statistics.csv`: estadísticas por segmento  
      (media, desviación estándar, media del baseline antes y después de la corrección y diferencia).

    **Pasos de procesamiento (alineados con `src/baseline.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_segmentation_info.bin` y extraer:

       - `eeg_segmented`
       - `channels`
       - `fs`
       - `segment_length_samples`
       - `n_segments`
       - `n_channels`

    2. **Cálculo de índices de baseline.**  
       Determinar los índices correspondientes al intervalo de baseline:

       ```math
       n_{\mathrm{start}} = \lfloor t_{\mathrm{start}} F_s \rfloor + 1
       ```

       ```math
       n_{\mathrm{end}} = \lfloor t_{\mathrm{end}} F_s \rfloor
       ```

       Por ejemplo, con

       ```math
       F_s = 500\,\text{Hz}
       ```

       el intervalo \([0, 0.1]\) s corresponde a las muestras **1–50**.

    3. **Corrección de baseline.**  
       Para cada canal \(c\) y segmento \(k\), calcular la media del baseline:

       ```math
       \mu_{c,k} = \text{mean}(\mathbf{E}[c, n_{\mathrm{start}}:n_{\mathrm{end}}, k])
       ```

       y restarla a todo el segmento:

       ```math
       \mathbf{E}^{\mathrm{corr}}[c,:,k] = \mathbf{E}[c,:,k] - \mu_{c,k}
       ```

    4. **Verificación.**  
       Calcular por segmento:

       - media  
       - desviación estándar  
       - media del baseline antes de la corrección  
       - media del baseline después de la corrección  

       Confirmar que la media del baseline posterior sea aproximadamente

       ```math
       0
       ```

       Opcionalmente generar un gráfico comparativo (por ejemplo, primer canal y primeros segmentos) con el intervalo de baseline marcado.

    5. **Guardado de resultados.**  
       - Serializar `eeg_baseline_corrected` en `eeg_1st_baseline_correction.bin`.
       - Construir el diccionario `dict_1st_baseline_correction` y guardarlo en `dict_1st_baseline_correction.bin`.
       - Exportar las estadísticas a `baseline_correction_statistics.csv`.

       Los datos resultantes se utilizan como entrada para la **Fase 5: rechazo de artefactos**.
"""

# ╔═╡ 97c5d3c1-41f1-42cb-ad01-192e7507c528
md"""
!!! info "Algoritmo Rechazo de artefactos basado en amplitud (`artifact_rejection.jl`)"

    **Entrada:**  
    ``E_bl ∈ ℝ^{C × L × K}`` (segmentos corregidos por baseline),  
    lista de canales, `N_used`, `A_min`, `A_max`

    1. Selección de canales a inspeccionar

       ```text
       C_used ← {1, … , min(N_used, C)}
       ```

    2. Inicialización

       ```text
       segment_is_good[1:K] ← true
       initialise arrays for:
           per-segment min amplitude
           per-segment max amplitude
           violating channel lists
       ```

    3. Evaluación de cada segmento

       ```text
       for k = 1 to K

           segment_used ← E_bl[C_used,:,k]

           m_k ← min(segment_used)
           M_k ← max(segment_used)

           if (m_k < A_min) or (M_k > A_max)

               segment_is_good[k] ← false

               V_k ← { c ∈ C_used :
                       min_n(E_bl[c,n,k]) < A_min
                       or
                       max_n(E_bl[c,n,k]) > A_max }

           else

               V_k ← ∅

           end

       end
       ```

    4. Identificación de segmentos válidos

       ```text
       good_indices ← { k : segment_is_good[k] = true }
       ```

    5. Construcción de la hipermatriz filtrada

       ```text
       E_AR ← E_bl[:,:,good_indices]
       ```

    **Salida:**  
    `E_AR` (EEG con segmentos artefactuales eliminados)  
    `segment_is_good` (vector booleano por segmento)  
    estadísticas por segmento y listas de canales que violan los umbrales
"""

# ╔═╡ 95d47537-1faf-47d0-ac9a-4be72a90c832
md"""
### Pipeline fase 5 (Artifact rejection)
"""

# ╔═╡ ece3c6d2-723c-4489-83eb-ff36cf9bc10e
md"""
!!! info "Pipeline Fase 5: Rechazo de artefactos"

    **Objetivo.**
    - Aplicar un criterio conservador basado en amplitud para descartar segmentos con artefactos evidentes.
    - Detectar saturaciones o ráfagas EMG extremas que puedan contaminar el análisis.
    - Conservar únicamente los segmentos limpios de EEG en reposo para los análisis posteriores de conectividad y espectro.

    **Entrada.**
    - EEG segmentado y corregido por baseline (Fase 4):  
      `data/baseline/dict_1st_baseline_correction.bin`, que contiene:
        - `eeg_baseline_corrected` (hipermatriz \(C \times L \times K\))
        - `channels`
        - `fs`
        - `segment_length_samples`
        - `n_segments`
        - `n_channels`

    - Configuración de rechazo de artefactos:

      ```math
      N_{\mathrm{used}} = 30
      ```

      canales inspeccionados.

      ```math
      A_{\min} = -70\,\mu\text{V}, \quad A_{\max} = +70\,\mu\text{V}
      ```

      Umbrales de amplitud.

      Los parámetros `before_event_ms` y `after_event_ms` se almacenan para posible uso futuro, pero no se aplican en esta implementación.

    **Salida.**
    - `data/artifact_rejection/eeg_artifact_rejected.bin`: hipermatriz filtrada con solo segmentos válidos

      ```math
      C \times L \times K_{\mathrm{good}}
      ```

    - `data/artifact_rejection/dict_artifact_rejection.bin`: diccionario con:
        - hipermatriz original y filtrada
        - lista de canales
        - umbrales de amplitud
        - índices de segmentos válidos
        - máscara booleana de segmentos buenos
        - conteos por segmento

    - `data/artifact_rejection/artifact_rejection_statistics.csv`: estadísticas por segmento (mínimos, máximos, flags y canales que violan los umbrales).

    **Pasos de procesamiento (alineados con `src/artifact_rejection.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_1st_baseline_correction.bin` y extraer:

       - `eeg_baseline_corrected`
       - `channels`
       - `fs`
       - `segment_length_samples`
       - `n_segments`
       - `n_channels`

    2. **Configuración de parámetros.**  
       Definir:

       ```math
       N_{\mathrm{used}},\; A_{\min},\; A_{\max}
       ```

       Calcular equivalentes en muestras para `before_event_ms` y `after_event_ms` (para uso futuro).  
       Seleccionar `channels_used` como los primeros

       ```math
       N_{\mathrm{used}}
       ```

       canales.

    3. **Detección de artefactos.**  
       Para cada segmento \(k\):

       - calcular

       ```math
       m_k = \min(E_{c,n,k})
       ```

       ```math
       M_k = \max(E_{c,n,k})
       ```

       considerando únicamente `channels_used`.

       El segmento se marca como inválido si

       ```math
       m_k < A_{\min} \quad \text{o} \quad M_k > A_{\max}
       ```

       Registrar los canales que violan los umbrales.

    4. **Cálculo de estadísticas.**  
       Construir un `DataFrame` con:

       - amplitud mínima por segmento  
       - amplitud máxima por segmento  
       - flag de segmento válido/inválido  
       - canales que superan los umbrales  

       Calcular también:

       - número de segmentos conservados  
       - número de segmentos rechazados  
       - porcentaje de rechazo  
       - mínimos y máximos globales.

    5. **Filtrado de segmentos.**  
       Construir `good_segment_indices` a partir de `segment_is_good`.

       Crear la hipermatriz filtrada:

       ```math
       \mathbf{E}^{\mathrm{AR}}
       ```

       seleccionando únicamente los segmentos válidos a lo largo de la tercera dimensión.

    6. **Visualización (opcional).**  
       Generar:

       - gráficos de segmentos de ejemplo con umbrales superpuestos  
       - histogramas de amplitudes máximas y mínimas por segmento con los umbrales marcados.

    7. **Guardado de resultados.**  
       - Serializar `eeg_artifact_rejected` en `eeg_artifact_rejected.bin`.
       - Construir el diccionario `dict_artifact_rejection` con todos los metadatos y guardarlo en `dict_artifact_rejection.bin`.
       - Exportar las estadísticas a `artifact_rejection_statistics.csv`.

       Los datos resultantes se utilizan como entrada para la **Fase 6: segunda corrección opcional de baseline**.
"""
# ╔═╡ 20569852-d89f-49f0-8d73-30ad7d2261d9
md"""
### Pipeline fase 6 (Baseline after artifacts)

Idem fase 4
"""

# ╔═╡ e79ba77a-2e08-47db-a1f9-91efe34c2510
md"""
!!! info "Algoritmo Análisis espectral mediante FFT (`FFT.jl`)"

    **Entrada:**  
    ``E_bl2 ∈ ℝ^{C × L × K}`` (segmentos con segunda corrección de baseline),  
    `F_s`, `N_fft = 512`,  
    parámetros de ventana (`P = 0.10`, Hamming `α = 0.54`, `β = 0.46`)

    1. Eliminación del componente DC

       ```text
       for c = 1 to C
           for k = 1 to K

               μ_ck ← mean(E_bl2[c,:,k])

               E_dc[c,:,k] ← E_bl2[c,:,k] − μ_ck

               store μ_ck for later DC reintroduction

           end
       end
       ```

    2. Aplicación de ventana

       ```text
       build taper window w[n]
           (Tukey/Hamming with length P · L at both ends)

       w2_mean ← mean(w^2)

       for c = 1 to C
           for k = 1 to K
               E_win[c,:,k] ← E_dc[c,:,k] ⊙ w
           end
       end
       ```

    3. Zero-padding

       ```text
       allocate E_pad ∈ ℝ^{C × N_fft × K}

       E_pad[:,1:L,:] ← E_win
       remaining samples ← 0

       Δf ← F_s / N_fft
       n_bins ← N_fft/2 + 1
       ```

    4. Cálculo de la rFFT

       ```text
       for c = 1 to C
           for k = 1 to K
               X[c,:,k] ← rFFT(E_pad[c,:,k])
           end
       end
       ```

    5. Reintroducción del componente DC

       ```text
       X[c,1,k] ← X[c,1,k] + μ_ck · N_fft
       ```

    6. Cálculo de potencia espectral

       ```text
       P[c,:,k] ← |X[c,:,k]|^2

       P[c,:,k] ← P[c,:,k] / w2_mean

       P[c,2:n_bins−1,k] ← 2 · P[c,2:n_bins−1,k]

       P[c,:,k] ← P[c,:,k] / N_fft^2
       ```

    7. Cálculo de potencia por bandas de frecuencia

       ```text
       for each frequency band b with range [f_min, f_max)

           I_b ← {ℓ : f_ℓ ∈ [f_min, f_max)}

           band_power[c,b] ← mean(P[c,I_b,:])

       end
       ```

    8. Promediado entre segmentos

       ```text
       P_mean ← mean(P, segments)
       ```

    **Salida:**  
    `P` (potencia por bin de frecuencia)  
    `P_mean` (potencia media por canal)  
    potencia por bandas de frecuencia  

    además de los archivos:

    - `data/FFT/dict_FFT.bin`  
    - `data/FFT/dict_FFT_power.bin`  
    - figuras y tablas generadas del análisis espectral
"""


# ╔═╡ c9854ffc-344b-4f87-bb93-2e71e71aa504
md"""
!!! info "Pipeline Fase 7: Análisis espectral (FFT)"

    **Objetivo.**
    - Estimar el espectro de potencia del EEG segmentado (tras la Fase 6) siguiendo un procedimiento compatible con BrainVision Analyzer (BVA).
    - Aplicar eliminación de DC, ventana tipo Hamming, zero-padding a 512 puntos, rFFT, corrección de varianza y normalización de potencia.
    - Calcular la potencia media por canal y la potencia por bandas de frecuencia para su análisis y visualización.

    **Entrada.**
    - EEG segmentado con segunda corrección de baseline:  
      `data/baseline/dict_2nd_baseline_correction.bin`, que contiene:
        - `eeg_2nd_baseline_corrected` (hipermatriz \(C \times L \times K\))
        - `channels`
        - `fs`
        - `segment_length_samples`
        - `n_segments`
        - `n_channels`

    **Salida.**
    - `data/FFT/dict_FFT.bin`: datos intermedios del cálculo espectral:
        - segmentos con ventana aplicada
        - zero-padding
        - ventana utilizada
        - valor medio de la ventana

          ```math
          \overline{w^2}
          ```

        - valores DC
        - vectores de frecuencia
        - tamaño FFT

          ```math
          N_{\mathrm{fft}}
          ```

        - resolución espectral

          ```math
          \Delta f
          ```

    - `data/FFT/dict_FFT_power.bin`: resultados espectrales:
        - potencia espectral

          ```math
          P
          ```

        - `P_mean`
        - `power_by_band`
        - `band_idx`
        - `freqs_bva`
        - definiciones de bandas

    - `results/figures/FFT/`: figuras generadas:
        - ventana Hamming
        - efecto del windowing
        - ejemplo de zero-padding
        - espectro FFT
        - espectro de potencia
        - rejilla de espectros por canal
        - espectro estilo BVA
        - topografía de potencia (por ejemplo banda alfa)

    - `results/tables/FFT/`:
        - `channel_summary.csv`
        - `power_by_band.csv`

    **Pasos de procesamiento (alineados con `src/FFT.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_2nd_baseline_correction.bin` y extraer:

       - EEG segmentado
       - `channels`
       - frecuencia de muestreo

       ```math
       F_s
       ```

       - número de muestras por segmento

       ```math
       L
       ```

       - número de segmentos

       ```math
       K
       ```

    2. **Eliminación de DC (estilo BVA).**  
       Para cada canal y segmento, restar la media del segmento.  
       Guardar los valores DC para su posterior reintroducción en el bin de frecuencia

       ```math
       0\,\text{Hz}
       ```

    3. **Aplicación de ventana.**  
       Construir una ventana tipo **Hamming** (Tukey con longitud de ventana del 10%).  
       Calcular

       ```math
       \overline{w^2}
       ```

       y aplicar la ventana a los segmentos sin DC.

    4. **Zero-padding.**  
       Extender cada segmento hasta

       ```math
       N_{\mathrm{fft}} = 512
       ```

       muestras.

       La resolución espectral resultante es

       ```math
       \Delta f = \frac{F_s}{N_{\mathrm{fft}}}
       ```

       Definir los bins de frecuencia desde

       ```math
       0
       ```

       hasta la frecuencia de Nyquist.

    5. **Cálculo de FFT y potencia.**  
       Calcular la **FFT real (rFFT)** de cada segmento.  
       Reintroducir el valor DC en el primer bin.

       Calcular la potencia espectral:

       ```math
       P = |X|^2
       ```

       Aplicar:

       - corrección de varianza

         ```math
         P \leftarrow P / \overline{w^2}
         ```

       - plegado del espectro completo (duplicar bins interiores)
       - normalización

         ```math
         P \leftarrow P / N_{\mathrm{fft}}^2
         ```

       La unidad final de potencia es

       ```math
       \mu\text{V}^2
       ```

    6. **Cálculo por bandas de frecuencia.**  
       Definir bandas:

       - Delta
       - Theta
       - Alpha
       - Beta (Low / Mid / High)
       - Gamma

       Para cada banda, promediar la potencia de los bins dentro del rango correspondiente y calcular `power_by_band`.

    7. **Promediado entre segmentos.**  
       Calcular la potencia media:

       ```math
       P_{\mathrm{mean}} = \mathrm{mean}(P, \text{segmentos})
       ```

       Opcionalmente calcular SNR (media / desviación estándar entre segmentos).

       Guardar resultados en:

       - `channel_summary.csv`
       - `power_by_band.csv`

    8. **Visualización.**  
       Generar figuras como:

       - forma de la ventana
       - ejemplo de señal con ventana aplicada
       - ejemplo de zero-padding
       - magnitud FFT
       - espectro de potencia (por ejemplo canal Oz)
       - espectros por canal en rejilla
       - espectro estilo BrainVision Analyzer
       - topografía de potencia alfa

       Guardar en `results/figures/FFT/`.

    9. **Guardado final.**  
       Serializar:

       - `dict_FFT.bin`
       - `dict_FFT_power.bin`
"""

# ╔═╡ e28cf442-e820-47ef-947c-85398eadadd7


# ╔═╡ 2c8d9444-4e8f-4234-b31d-c7e3425f8d13


# ╔═╡ 1c2bd8c6-3309-4a92-8a10-5437f5036c7a
md"""
!!! info "Algoritmo CSD mediante Laplaciano esférico de Perrin (`CSD.jl`)"

    **Entrada:**  
    ``E ∈ ℝ^{C × L × K}`` (segmentos EEG),  
    coordenadas de electrodos `(x_c, y_c, z_c)` desde TSV BIDS,  
    parámetros `m`, `degree`, `λ`

    1. Comprobación de consistencia canal–electrodo

       ```text
       load electrode TSV

       keep rows where type = EEG

       if any EEG channel missing in TSV
           abort "Channel without electrode position"
       end

       reorder electrode rows to match order of channels in E
       ```

    2. Normalización de posiciones en la esfera unitaria

       ```text
       for c = 1 to C

           r_c ← sqrt(x_c² + y_c² + z_c²)

           R_c ← (x_c / r_c, y_c / r_c, z_c / r_c)

       end
       ```

    3. Precomputación de coeficientes de Legendre

       ```text
       for k = 1 to degree

           c_g[k] ← (2k + 1) / (k^m · (k + 1)^m)

           c_h[k] ← (2k + 1) / (k^(m−1) · (k + 1)^(m−1))

       end
       ```

    4. Construcción de las matrices G y H (simétricas)

       ```text
       for i = 1 to C

           G[i,i] ← Σ_k c_g[k]
           H[i,i] ← Σ_k c_h[k]

           for j = i+1 to C

               cosγ ← R_i · R_j
               clamp cosγ to [−1, 1]

               G[i,j] ← Σ_k c_g[k] · P_k(cosγ)
               G[j,i] ← G[i,j]

               H[i,j] ← Σ_k c_h[k] · P_k(cosγ)
               H[j,i] ← H[i,j]

           end

       end
       ```

    5. Regularización y restricción de referencia libre

       ```text
       G_λ ← G + λ · I

       1_vec ← (1,1,…,1)ᵀ

       C ← G_λ^{-1} − (G_λ^{-1} · 1_vec · 1_vecᵀ · G_λ^{-1})
                        / (1_vecᵀ · G_λ^{-1} · 1_vec)

       L ← H · C
       ```

       `L` es el **operador CSD** de dimensión `C × C`.

    6. Aplicación del operador a cada segmento

       ```text
       for s = 1 to K

           E_CSD[:,:,s] ← L · E[:,:,s]

       end
       ```

    **Salida:**  
    `E_CSD` (EEG transformado a densidad de corriente),  

    archivos guardados:

    - `data/CSD/eeg_csd.bin`  
    - `data/CSD/dict_csd.bin`

    opcionalmente: controles de calidad  
    (RMS por canal, topografías CSD).
"""

# ╔═╡ 5b99601e-2fb5-4a21-9ab9-ead4df41c9e6
md"""
!!! info "Algoritmo Cálculo de wPLI por banda"

    **Entrada:** ``E_{CSD} \in \mathbb{R}^{(C \times L \times K)}``, banda `(f1,f2)`, `fs`

    1. Precalcular señales analíticas por canal

       ```text
       for c = 1 to C
           for s = 1 to K
               x        ← E_CSD[c,:,s]
               x_f      ← bandpass_filtfilt(x, fs, f1, f2; order=8)
               Z_c[:,s] ← analytic_signal(x_f)
           end
       end
       ```

    2. Construir matriz simétrica de wPLI

       ```text
       for i = 1 to C
           for j = i+1 to C
               Im_ij  ← Im(Z_i ⊙ conj(Z_j))
               num    ← |Σ Im_ij|
               den    ← Σ |Im_ij| + ε
               W[i,j] ← num / den
               W[j,i] ← W[i,j]
           end
       end
       ```

    **Salida:** `W` simétrica, diagonal 0
"""

# ╔═╡ 35f3cbd4-8a04-4f72-9c43-d73a23c6de8c
md"""
!!! info "Weighted Phase Lag Index (wPLI)"

    **Objetivo.** Calcular la conectividad funcional estática específica por banda utilizando wPLI a partir de segmentos transformados mediante CSD, agregando sobre todos los *samples* y segmentos (estilo BrainVision Analyzer), y guardar matrices de conectividad, listas de aristas y mapas de calor por banda.

    **Entrada**

    - Datos CSD: `data/CSD/dict_csd.bin` (`eeg_csd C × L × K`, canales, fs)
    - Definición de bandas: `data/FFT/dict_FFT_power.bin` o `dict_FFT.bin` (`bands_hz`); en caso de ausencia se utilizan bandas estándar (DELTA, THETA, ALPHA, BETA_LOW, BETA_MID, BETA_HIGH, GAMMA)

    **Salida**

    - `data/wPLI/dict_wpli.bin`: wpli, bandas, canales, fs, space = CSD  
    - `results/tables/wPLI/wPLI_{band}_matrix.csv`, `wPLI_{band}_edges.csv`: matriz de conectividad y lista de aristas por banda  
    - `results/figures/wPLI/wPLI_{band}_heatmap.png`: mapa de calor por banda  
    - `results/logs/wPLI/wPLI_{timestamp}.log`: registro de ejecución

    ---

    **Pasos de procesamiento** (alineados con `src/Connectivity/wPLI.jl`)

    1. Cargar `dict_csd.bin`.

    2. Cargar la definición de bandas desde `dict_FFT_power.bin` o `dict_FFT.bin`.  
       Si no existe, usar las bandas estándar: DELTA, THETA, ALPHA, BETA_LOW, BETA_MID, BETA_HIGH, GAMMA.

    3. Para cada banda `b` con rango `(f₁, f₂)`:

       **(a) Señales analíticas por canal**

       ```text
       for c = 1:C
           for s = 1:K
               x        ← eeg_csd[c,:,s]
               x_f      ← bandpass_filtfilt(x, fs, f1, f2; order=8)
               Z_c[:,s] ← analytic_signal(x_f)
           end
       end
       ```

       donde `Z_c ∈ ℂ^{L×K}`.

       **(b) Cálculo de wPLI por pares de canales**

       ```text
       for i = 1:C
           for j = i+1:C
               Im_ij = Im(Z_i ⊙ conj(Z_j))
               num   = |Σ Im_ij|
               den   = Σ |Im_ij| + ε
               W[i,j] = num / den
               W[j,i] = W[i,j]
           end
       end
       ```

       Matemáticamente:

       ```math
       \mathrm{wPLI}_{ij} =
       \frac{\left|\sum \mathrm{Im}_{ij}\right|}
            {\sum |\mathrm{Im}_{ij}| + \epsilon}
       ```

       **(c) Guardado de resultados**

       - Guardar matriz de conectividad CSV (`wPLI_{band}_matrix.csv`)
       - Guardar lista de aristas CSV (`from, to, wpli`)
       - Generar mapa de calor (`wPLI_{band}_heatmap.png`)

    4. Guardar `dict_wpli.bin` con resultados y metadatos.

    5. Ejecutar verificación:

       - valores en rango `[0,1]`
       - matriz simétrica
       - diagonal `0`
       - ausencia de `NaN` o `Inf`
"""

# ╔═╡ efc6b5dc-f489-48dc-bdbb-a79cdd744e86
md"
# Modelo surrogado

Para evaluar la significación estadística, se generaron señales EEG sustitutas de forma independiente para cada canal mediante **aleatorización de fase en el dominio de Fourier**. Si $X(f)$ representa la transformada de Fourier de la señal original, cada surrogate se construye como

```math
X^{(s)}(f) = |X(f)|\, e^{i\phi^{(s)}(f)},
\qquad
\phi^{(s)}(f) \sim \mathcal{U}(0, 2\pi)
```

imponiendo **simetría conjugada** para garantizar que la transformada inversa produzca una señal real. Con ello se preserva la **densidad espectral de potencia** de cada canal,

```math
|X^{(s)}(f)|^2 = |X(f)|^2
\tag{2}
```

mientras se destruyen las dependencias de fase entre canales y a lo largo del tiempo. En este *pipeline* se utilizaron $N_s = 200$ realizaciones por sujeto y condición.
"

# ╔═╡ 8713cdf3-d880-48f6-8575-a39bb16acb01
md"
## Implementación operativa

1. Calcular la FFT de cada canal y almacenar ``|X(f)|``.
2. Muestrear fases aleatorias ``\phi^{(s)}(f)`` de una distribución uniforme en $[0,2\pi]$.
3. Imponer simetría conjugada.
4. Reconstruir la señal sustituta mediante

```math
x^{(s)}(t) = \mathrm{IFFT}\!\left(|X(f)|\, e^{i\phi^{(s)}(f)}\right)
\tag{3}
```

5. Repetir el cálculo de wPLI sobre cada realización sustituta.

Finalmente, para cada conexión ``(i,j)``, el valor-p empírico se calculó como

```math
p_{ij}
=
\frac{
1 + \#\left\{\mathrm{wPLI}_{ij}^{(s)} \ge \mathrm{wPLI}_{ij}^{\mathrm{real}}\right\}
}{
1 + N_s
}
\tag{4}
```

de modo que, bajo la hipótesis nula, $p_{ij}$ es aproximadamente uniforme sobre el soporte discreto correspondiente.
"

# ╔═╡ fd21bfc3-e448-4b2c-a723-257a8d1aa6ca
md"""
!!! info "Generación de señales sintéticas mediante aleatorización de fase"

    **Objetivo.** Generar señales EEG sintéticas (*surrogates*) para construir distribuciones nulas empíricas.

    **Entrada:** Señales EEG preprocesadas.

    **Salida:** Matrices de conectividad sustitutas para la inferencia estadística.

    ---

    **Pasos**

    1. **Para cada canal** `i`:


         - Aleatorizar la fase de forma uniforme en `[0, 2\pi]` 
         - Imponer simetría conjugada  
         - Calcular la **FFT inversa** para obtener la señal sustituta

    2. Generar ``N_{\text{surr}}`` realizaciones (p. ej., `200`).

    3. Para cada conjunto de datos sustituto, repetir **Fase 7** (cálculo de `wPLI`).

    4. Almacenar las distribuciones de conectividad sustituta por banda y condición.
"""