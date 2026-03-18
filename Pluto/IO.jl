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

# ╔═╡ 1c2137c1-9117-405f-97dd-d2392fcc5911
md"""
**Packages input**
"""

# ╔═╡ f36b8ac2-d883-4d32-a7dc-2072e987165d
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

# ╔═╡ a73289ad-5d2c-4452-b19e-8a0ef3e42836
begin
	x = 0:0.01:10
	y = sin.(x)
	
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="x", ylabel="sin(x)")
	
	lines!(ax, x, y)
	
	fig
end

# ╔═╡ fb6880f3-5181-4299-af69-2675820284ca
md"""
# Conectividad Funcional EEG

En este **notebook** describimos un proceso reproducible, implementado en el lenguaje *Julia*, basado en señales registradas mediante EEG para estimar la **conectividad funcional estática** a través del **weighted Phase Lag Index (wPLI)**, que puede traducirse al castellano como **Índice de Retardo de Fase Ponderado**.

Las señales EEG fueron registradas en **estado de reposo**, tanto con **ojos abiertos (OA)** como con **ojos cerrados (OC)**, en **controles sanos** y en **pacientes con Esclerosis Múltiple (EM)**.

Este documento describe el diseño del estudio, las características de los participantes y el pipeline de análisis empleado para estimar la conectividad funcional.

---

## Características de los participantes

- **Duración promedio de la enfermedad:**
  - Total: **8.52 ± 5.74 años**
  - Mujeres: **8.23 ± 5.84 años**
  - Hombres: **9.14 ± 5.68 años**

- **Tratamientos modificadores de la enfermedad:**
  - **Ocrelizumab:** 38.63%
  - **Natalizumab:** 31.81%
  - **Alemtuzumab:** 13.63%
  - **Cladribina:** 13.63%
  - **Rituximab:** 2.27%
  - **Interferón:** 0%

---

## Diseño del estudio

El estudio combina **dos componentes complementarios**: uno **transversal (caso-control)** y otro **longitudinal**.

### Componente transversal (caso-control)

En el momento **T1** se realiza una comparación entre:

- **Pacientes con Esclerosis Múltiple (EM):** N = 44  
- **Grupo control sano:** N = 40  

Este análisis permite evaluar diferencias en la conectividad funcional entre pacientes y controles.

### Componente longitudinal

Dentro del grupo de pacientes con EM se realizó un **seguimiento longitudinal** entre dos momentos temporales:

- **T1:** N = 44  
- **Pérdida de seguimiento:** N = 14  
  - Deserción
  - Brote de la enfermedad
- **T2:** N = 30  

Este componente permite analizar **cambios intra-sujeto a lo largo del tiempo** en la conectividad funcional.

---

## Consideraciones para el análisis

El diseño del estudio implica que el pipeline de análisis debe permitir dos tipos de análisis complementarios:

- **Análisis entre grupos:** comparación **EM vs controles** en T1.
- **Análisis longitudinal intra-sujeto:** evaluación de **cambios entre T1 y T2** dentro del grupo de pacientes con EM.

Por tanto, la estructura del dataset y del pipeline se ha concebido para **incorporar explícitamente la dimensión longitudinal**, permitiendo integrar ambos tipos de análisis dentro de un mismo marco reproducible.

---

El marco metodológico está diseñado para garantizar **reproducibilidad completa**, integrando cálculos y análisis implementados en *Julia* y visualizaciones generadas mediante el paquete *Makie*.
"""

# ╔═╡ b8e9ae7c-942f-4f7c-9ed8-02dbbd0cb657
md"""
# Esclerosis Múltiple (breve intro)

La **Esclerosis Múltiple (EM)** es una enfermedad **crónica, inflamatoria y desmielinizante** del **sistema nervioso central (SNC)** que afecta a la **vaina de mielina** que recubre las fibras nerviosas. Este daño altera la transmisión normal de los impulsos nerviosos y puede provocar una amplia variedad de **síntomas neurológicos y discapacidades**.

La etiología exacta de la EM aún no se conoce completamente, aunque se considera generalmente una **enfermedad autoinmune** en la que el sistema inmunitario ataca de forma errónea la mielina del SNC.

## Síntomas principales

Los síntomas de la EM pueden variar considerablemente entre pacientes, pero algunos de los más frecuentes incluyen:

- **Fatiga**
- **Debilidad muscular**
- **Problemas de visión**
- **Deterioro cognitivo**

## Diagnóstico

El diagnóstico de la EM se basa en una combinación de:

- **síntomas clínicos**
- hallazgos en **resonancia magnética (RM o MRI)**
- **pruebas de laboratorio**

Estos criterios se aplican siguiendo **guías diagnósticas internacionalmente aceptadas**, que permiten confirmar la presencia de lesiones desmielinizantes en el sistema nervioso central.

## Tratamiento

El tratamiento de la EM suele combinar varias estrategias terapéuticas, entre ellas:

- **terapias modificadoras de la enfermedad**
- **medicación sintomática**
- **rehabilitación física**
- **adaptaciones del estilo de vida**

## Pronóstico

El curso de la enfermedad es **altamente variable entre pacientes**. En muchos casos, la EM sigue un curso **progresivo**, que puede conducir a discapacidades significativas, incluyendo dificultades para **caminar**, **hablar** o realizar **actividades de la vida diaria**.

## Recursos y datos en España

Para información y datos sobre la realidad de la EM desde la perspectiva de los pacientes en España —incluyendo epidemiología, comorbilidades y necesidades de cuidadores— la plataforma [**EMData**](https://emdata.esclerosismultiple.com) **(Esclerosis Múltiple España)** ofrece un recurso centralizado: 
"""

# ╔═╡ 3b1ca46d-580a-4a4c-83fc-d5a8249532a0
md"
# Electroencefalografía (EEG)

La **electroencefalografía (EEG)** es una técnica no invasiva que registra la actividad eléctrica del cerebro mediante electrodos colocados sobre el cuero cabelludo. Estos electrodos se distribuyen siguiendo sistemas estandarizados, como el **sistema internacional 10–20**, que garantiza un muestreo espacial reproducible de la actividad cortical.

El conjunto de datos comprende registros de **EEG** en **estado de reposo** procedentes de dos grupos: **controles sanos** y **pacientes diagnosticados con esclerosis múltiple (EM)**. La pertenencia a cada grupo se definió de acuerdo con los criterios clínicos proporcionados en los metadatos asociados.

Las variables demográficas (**edad y sexo**) se registraron para todos los participantes y se incorporaron a la estructura de metadatos a nivel de sujeto, con el fin de permitir análisis estadísticos posteriores.

En este sistema, los electrodos se posicionan a distancias equivalentes al **10 % o 20 % de determinadas medidas craneales**, lo que permite una cobertura sistemática de la superficie del cráneo. Las etiquetas de los electrodos reflejan aproximadamente la región cortical subyacente:

- **F**: frontal  
- **C**: central  
- **T**: temporal  
- **P**: parietal  
- **O**: occipital  

Además, los **números impares** indican posiciones en el **hemisferio izquierdo**, los **números pares** posiciones en el **hemisferio derecho**, y la letra **z** denota posiciones en la **línea media**.

## Origen fisiológico de la señal EEG

Las señales registradas mediante EEG reflejan principalmente la **suma de potenciales postsinápticos** generados por grandes poblaciones de neuronas corticales que se activan de forma relativamente sincronizada. En particular, las contribuciones más relevantes provienen de **neuronas piramidales** orientadas perpendicularmente a la superficie cortical.

Estas neuronas generan **campos eléctricos dipolares** que se propagan a través del tejido cerebral, el líquido cefalorraquídeo, el cráneo y el cuero cabelludo mediante un proceso conocido como **conducción de volumen** (*volume conduction*).

Como consecuencia, el **EEG de superficie** no mide potenciales de acción individuales, sino una **actividad sináptica integrada espacial y temporalmente** que proviene de áreas corticales relativamente amplias. Esto implica:

- **Resolución espacial limitada**
- **Posible mezcla de señales entre canales**

## EEG en estado de reposo

El **EEG en estado de reposo** (*resting-state EEG*, **RS-EEG**) es un paradigma en el que la actividad cerebral se registra mientras el sujeto permanece relajado y no realiza ninguna tarea estructurada. Este enfoque captura las dinámicas oscilatorias basales y las fluctuaciones espontáneas de la actividad neuronal, evitando los factores de confusión asociados al diseño de tareas y a la carga motora o cognitiva.

En neurofisiología clínica, los registros de **RS-EEG** se utilizan ampliamente para:

- evaluar **ritmos de fondo**
- analizar la **reactividad cerebral**
- detectar **asimetrías**
- identificar **anomalías focales o difusas**

## Relevancia en esclerosis múltiple

En el contexto de la **esclerosis múltiple (EM)**, el RS-EEG resulta particularmente interesante porque la enfermedad afecta a **redes cerebrales distribuidas**. Estas alteraciones pueden manifestarse como cambios en los **patrones oscilatorios** y en la **conectividad funcional**, incluso en ausencia de tareas cognitivas explícitas.

Además, el EEG constituye una herramienta **relativamente simple, portátil y de bajo coste**, lo que lo convierte en un complemento útil a técnicas de neuroimagen como la **resonancia magnética (RM)** en estudios longitudinales o comparativos.

## Condiciones experimentales

Los registros suelen obtenerse en dos condiciones principales:

- **Ojos abiertos (EO)**  
  mantiene un nivel moderado de activación cortical y reduce la dominancia del ritmo alfa.

- **Ojos cerrados (EC)**  
  potencia típicamente el **ritmo alfa**, siendo sensible a cambios dependientes del estado cerebral y potencialmente relacionados con la enfermedad.

La comparación entre **EO** y **EC** permite evaluar tanto efectos asociados al **nivel de alerta** como posibles **alteraciones neurofisiológicas** entre controles sanos y pacientes con EM.
"

# ╔═╡ 7b611338-d3c2-4e05-bf4e-bb87b1b56e5f
md"
## Parámetros de adquisición EEG

Las señales de **EEG** se registraron a partir de $C = 32$ electrodos de cuero cabelludo colocados según el **sistema internacional (SI) 10–20**. La disposición incluye **31 canales de EEG**, además del **electrodo de referencia REF (FCz)** y el **electrodo de tierra GRND (Fpz)**.

Las señales se muestrearon a una frecuencia de $f_s = 500\,\text{Hz}$, lo que produce

$N = f_s \cdot T = 50\,000$

muestras por canal y por condición.

Los datos se representan como matrices

$X_{s,k} \in \mathbb{R}^{C \times N},$

donde $s$ indexa a los sujetos y $k \in \{\text{EO}, \text{EC}\}$ denota la condición de registro. Cada fila corresponde a la serie temporal de un electrodo, expresada en microvoltios ($\mu V$).
"

# ╔═╡ dadc3838-9f52-4397-a36d-e74626f873cb
md"
## Software BrainVision

**BrainVision Analyzer** es un software comercial ampliamente utilizado para el procesamiento y análisis de señales de **EEG**, desarrollado por **Brain Products**.

### Estructura típica de archivos BrainVision

Está diseñado para trabajar de forma nativa con el **formato BrainVision**, el cual separa la información del registro en tres archivos complementarios:

```
recording.vhdr
recording.eeg
recording.vmrk
```

- **`recording.eeg`**  
  Contiene los datos binarios, las señales EEG crudas, almacenadas como datos binarios.  
  Cada canal representa la evolución temporal del potencial eléctrico medido en el cuero cabelludo.

- **`recording.vhdr`**  
  Archivo de texto estructurado en secciones que incluye:
  - información del sistema de adquisición
  - número de canales
  - frecuencia de muestreo
  - configuración de referencia
  - formato del archivo binario
  - **referencias a los archivos `.eeg` y `.vmrk`**

- **`recording.vmrk`**  
  Archivo de texto con la lista de eventos experimentales (por ejemplo estímulos o inicio de ensayo).  
  Cada marcador suele incluir:
  - tipo de evento
  - posición temporal
  - duración
  - canal asociado (si aplica)


| Archivo | Tipo | Función |
|-------|------|--------|
| `.eeg` | Binario | Contiene la **serie temporal cruda del EEG** (voltaje continuo en $\mu V$). |
| `.vhdr` | Texto | Archivo de **cabecera** con información sobre el registro: número de canales, frecuencia de muestreo, referencias, formato de datos y enlaces a los otros archivos. |
| `.vmrk` | Texto | Archivo de **marcadores de eventos** que describe los tiempos y tipos de eventos experimentales. |

Esta estructura modular permite una gestión clara de la información y facilita su integración en flujos de análisis reproducibles.

Uno de los aspectos clave de este formato es la representación de los datos en el archivo `.eeg`, donde las señales pueden almacenarse en distintos tipos numéricos. Entre ellos, el formato **`IEEE_FLOAT_32`** (coma flotante de 32 bits) es especialmente relevante, ya que permite una alta precisión en la representación de la señal **EEG**, evitando pérdidas de información durante el registro y el procesamiento. Este nivel de precisión es importante para análisis sensibles, como la detección de componentes, el filtrado o el modelado de fuentes.

El archivo `.vhdr` especifica, entre otros parámetros, el tipo de datos (`BinaryFormat=IEEE_FLOAT_32`), la frecuencia de muestreo, el número de canales y las referencias a los archivos asociados. Por su parte, el archivo `.vmrk` contiene la información temporal de los eventos experimentales, esencial para segmentar y analizar la señal en función del diseño del experimento.

En conjunto, BrainVision Analyzer no solo proporciona herramientas avanzadas de procesamiento (filtrado, ICA, análisis espectral, entre otros), sino que también se apoya en un formato de datos transparente y bien documentado. Esto facilita su compatibilidad con estándares como **EEG-BIDS** y su integración con otras herramientas de análisis dentro del ecosistema de neurociencia computacional.

Estos tres archivos forman un **conjunto inseparable**. El archivo `.vhdr` actúa como punto de entrada para la mayoría del software de análisis, ya que describe cómo interpretar el archivo binario `.eeg` y dónde encontrar los eventos en `.vmrk`.

---

### Enlaces internos entre archivos

Dentro de los archivos `.vhdr` y `.vmrk` existen **campos que apuntan a los otros archivos del conjunto**.

Ejemplo típico dentro de un archivo `.vhdr`:

```
DataFile=recording.eeg
MarkerFile=recording.vmrk
```

En el archivo `.vmrk` también aparece una referencia al archivo de datos:

```
DataFile=recording.eeg
```

Estos enlaces permiten que los programas de análisis localicen correctamente los datos y los eventos.
"

# ╔═╡ aa3253d3-8a64-4e0f-9e68-dd9b1c07db2d
md"
## BIDS

El **Brain Imaging Data Structure (BIDS)** es un estándar desarrollado por la comunidad de neuroimagen para organizar, describir y compartir datos experimentales de forma consistente entre laboratorios y herramientas de análisis. 

Su objetivo principal es reducir la heterogeneidad en los formatos y estructuras de datos, facilitando su **reutilización, interoperabilidad y reproducibilidad** en estudios científicos.

En el caso del **EEG**, la extensión **EEG-BIDS** adapta este estándar a las particularidades de los registros electroencefalográficos, definiendo convenciones claras para nombres de archivos, organización por sujeto/sesión y el uso de archivos de metadatos (`.json`, `.tsv`). Además, promueve el uso de formatos ampliamente soportados como **BrainVision** (`.vhdr`, `.eeg`, `.vmrk`), lo que permite integrar datos provenientes de distintos sistemas en un flujo de trabajo común.

---

### Validación de datasets con BIDS Validator

El **BIDS Validator** es una herramienta diseñada para comprobar automáticamente si un dataset cumple con la especificación **BIDS**. Su función principal es detectar errores en la estructura de carpetas, nombres de archivos y metadatos, así como identificar información faltante o inconsistente, lo que ayuda a garantizar que los datos sean interoperables y reutilizables.

El validador puede ejecutarse tanto desde la línea de comandos como directamente en el navegador, lo que facilita su uso en diferentes entornos de trabajo. Al utilizarlo durante el proceso de organización de datos, es posible corregir problemas de forma temprana y asegurar la compatibilidad con herramientas de análisis como **BrainVision Analyzer**, **EEGLAB** o **MNE-Python**.

Puedes acceder al validador en el siguiente enlace:  
[BIDS Validator](https://bids-standard.github.io/bids-validator/)

---

### Introducción a BIDS y su uso con BrainVision Analyzer

Al adaptar un dataset al estándar **BIDS (Brain Imaging Data Structure)**, no basta con renombrar los archivos siguiendo la convención establecida. En el caso del formato **BrainVision** (`.vhdr`, `.eeg`, `.vmrk`), es necesario también actualizar las **referencias internas** contenidas en los archivos de texto, ya que estos incluyen enlaces explícitos entre sí. Si los nombres cambian pero estas referencias no se actualizan, el software de análisis no podrá localizar correctamente los datos asociados, generando errores en la carga o interpretación del dataset.

La integración de **BIDS** con herramientas como **BrainVision Analyzer** es especialmente natural, dado que este software utiliza el formato BrainVision, uno de los formatos recomendados dentro de EEG-BIDS. No obstante, convertir un dataset a BIDS implica más que reorganizar archivos: supone estructurar la información de forma coherente, estandarizar metadatos y asegurar la consistencia entre todos los componentes del registro.

De este modo, **BIDS** no solo facilita el uso de los datos en BrainVision Analyzer, sino que también permite su interoperabilidad con otras herramientas ampliamente utilizadas como EEGLAB, FieldTrip o MNE. En conjunto, BIDS actúa como una capa de estandarización que favorece la reproducibilidad, mientras que BrainVision Analyzer se integra como una herramienta de procesamiento dentro de un ecosistema más amplio y compatible.

---

### Problema al convertir a BIDS

En el estándar **BIDS-EEG**, los archivos deben renombrarse siguiendo una convención estructurada, por ejemplo:

```
sub-01_ses-01_task-rest_run-01_eeg.vhdr
sub-01_ses-01_task-rest_run-01_eeg.eeg
sub-01_ses-01_task-rest_run-01_eeg.vmrk
```

Sin embargo, cuando se cambian los nombres de archivo:

- los **enlaces internos dentro de `.vhdr` y `.vmrk` siguen apuntando al nombre antiguo**
- esto provoca que el software no encuentre los archivos correspondientes.

---

### Actualización necesaria de enlaces

Después de renombrar los archivos, se deben actualizar las referencias internas.

En el archivo `.vhdr`:

```
DataFile=sub-01_ses-01_task-rest_run-01_eeg.eeg
MarkerFile=sub-01_ses-01_task-rest_run-01_eeg.vmrk
```

En el archivo `.vmrk`:

```
DataFile=sub-01_ses-01_task-rest_run-01_eeg.eeg
```

---

### Importancia para el análisis de EEG

Actualizar estos enlaces garantiza que:

- el software compatible con **BrainVision** (por ejemplo EEGLAB, FieldTrip o MNE)
  pueda cargar correctamente los datos
- los **marcadores de eventos** se alineen correctamente con la señal EEG
- los **validadores BIDS** puedan interpretar el dataset sin errores

En resumen, la actualización de los enlaces internos asegura la **consistencia estructural del dataset** y permite que las herramientas de análisis y validación localicen correctamente tanto la señal EEG como los eventos experimentales asociados.

Para cada sujeto y sesión, **BIDS** proporciona un archivo `electrodes.tsv` (por ejemplo, `sub-M05_ses-T2_electrodes.tsv`) que contiene las columnas:

(`name`, `x`, `y`, `z`, `type`)

Las coordenadas $(x, y, z)$ están **normalizadas** (por ejemplo, sobre una **esfera unitaria** o proporcionalmente a las dimensiones de la cabeza) y se utilizan para:

- generar **mapas topográficos**
- realizar **modelado de fuentes**

La **Tabla** enumera los nombres de los electrodos, sus coordenadas y el tipo (`EEG`, `REF` o `GND`) tal como aparecen en el conjunto de datos.

En esta fase de Input-Ouput **Fase (IO)**, donde las señales *EEG* registradas mediante el equipo comercial (en este caso **BrainVision**) se han organizado como un diccionario o matriz indexada por canales, junto con metadatos temporales y de canal, y estadísticas opcionales de calidad obtenidas durante la inspección inicial.
"

# ╔═╡ 43e0934e-d7e9-421d-a678-ee7766bc90bc
md"### Tabla de electrodos (name, x, y, z, type)"

# ╔═╡ c62af5d5-c30d-44b7-a067-e01bbdc5e835
begin
	electrodes_path = joinpath(@__DIR__, "..", "data", "electrodes", "sub-M05_ses-T2_electrodes.tsv")
	electrodes = if isfile(electrodes_path)
		CSV.read(electrodes_path, DataFrame; delim = '\t')
	else
		@warn "No se encontró electrodes.tsv; usando tabla vacía (modo público/export)." electrodes_path
		DataFrame(name = String[], x = Float64[], y = Float64[], z = Float64[], type = String[])
	end

	electrodes
end

# ╔═╡ d4f94266-9f5f-469c-abc9-90ed2995b95a
md"
Tipos de archivos **BIDS**: ubicación, nombre y contenido principal

| Localización | Archivo | Descripción |
|----------|------|-------------|
| Root | `dataset_description.json` | Metadatos generales del dataset: nombre, autores, versión de BIDS y licencia. |
| Root | `participants.tsv` | Demografía tabular: `participant_id` (p. ej., `sub-M05`), `age`, `sex`, `group` (p. ej., `control`, `MS`). |
| Root | `participants.json` | Descripción y unidades para cada columna en `participants.tsv`. |
| Root | `README.md` | Documentación del dataset legible por humanos. |
| `sub-*/ses-*/eeg/` | `sub-*_ses-*_task-*_run-*_eeg.eeg` | Serie temporal EEG binaria cruda (voltaje continuo, $\mu V$). |
| `sub-*/ses-*/eeg/` | `*_eeg.vhdr` | Cabecera BrainVision: referencias a `.eeg` y `.vmrk`, información de canales y frecuencia de muestreo. Tras renombrar en BIDS, `DataFile=` y `MarkerFile=` deben referenciar los nuevos nombres base. |
| `sub-*/ses-*/eeg/` | `*_eeg.vmrk` | Archivo de marcadores BrainVision: tiempos de eventos y tipos. `DataFile=` debe apuntar al nuevo nombre del archivo `.eeg`. |
| `sub-*/ses-*/eeg/` | `*_eeg.json` | *Sidecar* BIDS-EEG con metadatos: `SamplingFrequency` (p. ej., 500), `PowerLineFrequency` (50), `EEGReference` (p. ej., `TP9-TP10` linked), `EEGGround` (p. ej., `Fpz`), `EEGChannelCount` (31), `RecordingType` (continuous), entre otros campos requeridos por BIDS-EEG. |
| `sub-*/ses-*/eeg/` | `*_events.tsv` | Eventos del *run*: columnas `onset`, `duration`, `trial_type` (p. ej., `rest`); tiempos en segundos (decimal). |
| `sub-*/ses-*/eeg/` | `*_events.json` | Descripción y unidades de las columnas presentes en `events.tsv`. |
"

# ╔═╡ b15557a3-4733-4bc7-85c7-313b08aa1039
md"""
# Pipeline

El **Pipeline** de análisis de señales **EEG** en **Julia** para el proyecto **BRAIN** (investigación universitaria), incluye carga de datos, preprocesamiento, filtrado, ICA, segmentación, corrección de baseline, rechazo de artefactos, análisis espectral (FFT) y medidas de conectividad (CSD, wPLI).

## Requisitos

- **Julia** ≥ 1.12 (recomendado; el `Manifest.toml` está generado con 1.12.3)
- Dependencias del proyecto (se instalan con `Pkg.instantiate()`)

## Instalación

1. Clona o abre el proyecto y entra en su directorio:
   ```bash
   cd EEG_Julia
   ```

2. Activa el entorno del proyecto e instala dependencias:
   ```bash
   julia --project=.
   ]
   instantiate
   backspace
   ```

   O desde un script:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Estructura del proyecto

```
EEG_Julia/
├── config/
│   └── default_config.jl    # Configuración: rutas, filtros, ICA, segmentación, FFT
├── script/
│   └── EEG.jl               # Punto de entrada (carga config + IO)
├── src/
│   ├── IO.jl                # Carga de datos raw y preprocesamiento inicial
│   ├── filtering.jl         # Filtrado (notch, bandreject, highpass, lowpass)
│   ├── ICA.jl               # Descomposición ICA (FastICA)
│   ├── ICA_cleaning.jl      # Evaluación y eliminación de componentes artefactuales
│   ├── segmentation.jl      # Segmentación en épocas
│   ├── baseline.jl          # 1ª corrección de baseline
│   ├── artifact_rejection.jl # Rechazo de artefactos por amplitud
│   ├── baseline_2st.jl      # 2ª corrección de baseline
│   ├── FFT.jl               # Análisis espectral (FFT, potencia por bandas)
│   └── Connectivity/
│       ├── CSD.jl           # Current Source Density (spline esférico Perrin)
│       └── wPLI.jl          # Weighted Phase Lag Index por bandas
├── data/
│   ├── raw/                 # Datos EEG en TSV + metadata + electrodos
│   ├── IO/                  # Salida del paso de carga
│   ├── filtering/           # Salidas de cada filtro
│   ├── ICA/                 # Resultados ICA y datos limpios
│   ├── segmentation/        # Hipermatriz segmentada
│   ├── baseline/            # Datos tras 1ª y 2ª baseline
│   ├── artifact_rejection/  # Segmentos tras rechazo de artefactos
│   └── (CSD, etc. según módulos)
├── results/
│   ├── figures/             # Gráficos (PSD, ICA, CSD, wPLI, etc.)
│   ├── tables/              # CSV y tablas (FFT, wPLI, estadísticas)
│   └── logs/                # Logs por fecha (CSD, wPLI, etc.)
├── Project.toml
├── Manifest.toml
└── README.md
```

## Flujo del pipeline

El procesamiento sigue un orden secuencial; cada paso lee la salida del anterior:

| Orden | Módulo               | Descripción breve |
|------:|----------------------|-------------------|
| 1     | **IO**               | Carga TSV raw, organiza canales (diccionario), PSD, calidad de canales. |
| 2     | **filtering**        | Notch 50 Hz, bandreject ~100 Hz, highpass 0.5 Hz, lowpass 150 Hz (Butterworth, filtfilt). |
| 3     | **ICA**              | FastICA simétrico sobre datos filtrados; guarda componentes y matrices de mezcla. |
| 4     | **ICA_cleaning**     | Evaluación automática de ICs (features/scores), eliminación de artefactos, reconstrucción. |
| 5     | **segmentation**     | División en segmentos de longitud fija (ej. 1 s), sin solapamiento; hipermatriz 3D. |
| 6     | **baseline**         | 1ª corrección de baseline por segmento (intervalo 0–0.1 s). |
| 7     | **artifact_rejection** | Rechazo de épocas por umbral de amplitud (±70 µV). |
| 8     | **baseline_2st**     | 2ª corrección de baseline (mismo intervalo). |
| 9     | **FFT**              | FFT, ventana Hamming, zero-padding; potencia espectral y por bandas (Delta–Gamma). |
| 10    | **CSD**              | Current Source Density (Laplaciano esférico Perrin) sobre datos segmentados. |
| 11    | **wPLI**             | Conectividad wPLI por bandas de frecuencia (entrada: datos CSD). |

Los datos intermedios se guardan en `data/` (p. ej. `.bin` serializados); figuras y tablas en `results/figures/` y `results/tables/`.
"""

# ╔═╡ 15966fd1-e39d-4c1c-ba38-f10a09e598dd
md"
# Raw data (datos crudos)

Antes de cualquier etapa de preprocesamiento es recomendable realizar una **inspección básica del EEG crudo**.  
Esto permite detectar problemas de adquisición, artefactos evidentes o canales defectuosos que podrían afectar al análisis posterior.

Las siguientes representaciones gráficas y estadísticas son útiles para realizar un **control de calidad inicial (QC)** de los datos.
"

# ╔═╡ 8d4caa60-3407-4510-9db4-8132946572fb
md"
## Carga datos raw`.tsv`

Se carga el archivo `.tsv` que contiene los datos **EEG raw**, en este caso solamente para un sujeto denominado **M05**, en condiciones de Ojos Cerrados **(EC)**: 

`sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv`

El formato esperado es:

- Primera columna: nombres de canales **(Channel)**
- Columnas siguientes: muestras temporales (una columna por punto de tiempo)
- Los datos están en unidades de microvoltios **(µV)**

Este archivo es la entrada principal del pipeline y contiene los datos sin ningún procesamiento previo **(raw data)**. El `.tsv` es el resultado de transformar los archivos de **BrainVision** al formato **BIDS**.
"

# ╔═╡ 4001ae18-0212-4445-8950-75c3c21f40d3
begin
println("=" ^ 60)
println("📊 CARGA DE DATOS EEG")
println("=" ^ 60)

# Ruta al archivo de datos raw
# Se usa @__DIR__ para obtener la ruta relativa al archivo del script
dir_raw = joinpath(@__DIR__, "..", "data", "raw", "sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv")

# Leer el archivo TSV como DataFrame
# CSV.read carga el archivo y lo convierte en una estructura tabular
data_raw = if isfile(dir_raw)
	CSV.read(dir_raw, DataFrame)
else
	@warn "No se encontró el TSV raw; usando datos demo (modo público/export)." dir_raw
	n = 2000
	t = range(0, step = 1 / 500, length = n)
	DataFrame(
		Channel = ["Cz", "Pz", "Fz", "Oz"],
		[Symbol("s$(i)") => [sin(2π * 10 * t[i]) + 0.05randn(), sin(2π * 8 * t[i]) + 0.05randn(),
		                    sin(2π * 12 * t[i]) + 0.05randn(), sin(2π * 6 * t[i]) + 0.05randn()]
		 for i in 1:n]...,
	)
end

println("✓ Archivo cargado: $(basename(dir_raw))")
println("✓ Dimensiones: $(size(data_raw, 1)) canales × $(size(data_raw, 2) - 1) muestras")
println("✓ Tabla de datos:")
display(data_raw)
end

# ╔═╡ 3ca6e623-ba24-413d-b383-4f82bbe1de87
md"
**Información temporal del registro EEG**

En este caso se extrae información sobre las características temporales del registro **EEG**, como son:

- Frecuencia de muestreo **(fs)**: número de muestras por segundo
- Número de muestras: longitud temporal del registro
- Duración total **(t seg)**: tiempo total de registro en segundos y minutos
"

# ╔═╡ f6f00a09-de36-4617-b669-a633ed9c9664
begin
println("=" ^ 60)
println("⏱️  INFORMACIÓN TEMPORAL")
println("=" ^ 60)

# Frecuencia de muestreo (Hz)
# Define la resolución temporal: cuántas muestras se toman por segundo
fs = 500

# Número de muestras temporales
# Se resta 1 porque la primera columna es "Channel" (nombres de canales)
n_muestras = size(data_raw, 2) - 1

# Vector de tiempo en segundos
# Se genera desde 0 hasta (n_muestras-1) y se divide por fs para obtener segundos
tiempo_seg = collect(0:(n_muestras-1)) ./ fs

# Duración total del registro
duracion_total = tiempo_seg[end]

println("  Frecuencia de muestreo: $(fs) Hz")
println("  Número de muestras: $(n_muestras)")
println("  Duración total: $(round(duracion_total, digits=2)) segundos ($(round(duracion_total/60, digits=2)) minutos)")
end

# ╔═╡ 42756cc1-b01a-43bc-847a-9263233984fd
md"
**Extracción y organización de canales**

Se extraen los nombres de los canales y se organizan los datos en estructuras
convenientes para el procesamiento, de forma siguiente:

- **channels**: vector con nombres de canales
- **EEG_matrix**: matriz (canales × muestras) con todos los datos
- **dict_EEG**: diccionario que mapea nombre de canal → señal (vector)

El diccionario es especialmente útil porque permite acceso por nombre de canal
y facilita el procesamiento posterior en el pipeline.
"

# ╔═╡ 44d68fb3-8b55-45ae-a1dc-834acde8556e
begin
println("=" ^ 60)
println("📡 INFORMACIÓN DE CANALES")
println("=" ^ 60)

# Nombres de canales (primera columna del DataFrame)
channels = data_raw.Channel

# Matriz de datos EEG (canales × muestras)
# Se excluye la columna "Channel" usando Not(:Channel)
# Cada fila corresponde a un canal, cada columna a una muestra temporal
EEG_matrix = Matrix(data_raw[:, Not(:Channel)])

println("  Número de canales: $(length(channels))")
println("  Canales: $(join(channels, ", "))")
println()

# Diccionario: canal → señal (vector de muestras)
# Esta estructura facilita el acceso por nombre de canal y es compatible
# con el resto del pipeline que espera datos organizados de esta manera
dict_EEG = Dict(data_raw.Channel[i] => EEG_matrix[i, :]
                for i in 1:size(data_raw, 1))
end

# ╔═╡ 937f5801-e2d0-411c-a7e9-2d8358e8ac95
md"
**Guardado de datos organizados**

Se guarda el diccionario de canales en formato binario nativo de *Julia* (`dict_EEG.bin`). Este formato es eficiente para matrices grandes y permite serialización rápida de estructuras complejas. El archivo guardado será la entrada para el siguiente paso del pipeline (filtrado, corrección de línea base, etc.).
"

# ╔═╡ 64f4f3c1-8b34-48fb-8a5f-3c414f09ad17
begin
# Directorio de salida para datos IO
dir_io = joinpath(@__DIR__, "..", "data", "IO")
path_dict = joinpath(dir_io, "dict_EEG.bin")

# Asegurar que el directorio existe
isdir(dir_io) || mkpath(dir_io)

# Serializar y guardar el diccionario
Serialization.serialize(path_dict, dict_EEG)

println("💾 Diccionario EEG raw guardado en: $(abspath(path_dict))")
end

# ╔═╡ 80bfe9b0-af73-4cb8-be31-65455652ad59
md"
## Control de calidad señales EEG 

Es necesario generar visualizaciones exploratorias en el dominio temporal/frecuencial. Estos gráficos se generan a partir de los **datos raw**. Estas visualizaciones ayudan a identificar:

- Artefactos obvios (parpadeos, movimiento, ruido de línea).
- Diferencias entre canales.
- Calidad general del registro.

Algunas notas prácticas para control de calidad de las señales **EEG**:

- Un **canal plano** suele indicar desconexión del electrodo o fallo del amplificador.
- **Picos estrechos** en **50/60 Hz** indican interferencia de red eléctrica.
- **Exceso de potencia <1 Hz** puede deberse a mal contacto electrodo–piel o movimiento.
- En registros con **ojos cerrados (EC)**, suele observarse aumento de potencia en **banda alfa (8–12 Hz)** en regiones occipitales.
- Es recomendable marcar o excluir **canales sospechosos** antes de realizar análisis espectral o conectividad.
"

# ╔═╡ 96c6f9db-a65b-4769-8ea7-dcd4b26a4378
md"
Tabla x: Gráficas y estadísticas para la inspección y control de calidad de **EEG crudo**

| Representación o estadístico | Propósito y uso |
|---|---|
| **Trazas crudas** (EEG vs. tiempo, vista multicanal apilada o desplazable) | Visualizar el comportamiento temporal de los canales, identificar artefactos o canales planos, comprobar la continuidad de la señal y localizar eventos; primera verificación de la integridad y sincronía de la grabación. |
| **Mapa de impedancias** | Refleja la calidad del contacto electrodo–piel. Impedancias altas (>10–20 kΩ) suelen anticipar ruido o inestabilidad. Un mapa topográfico ayuda a detectar áreas problemáticas (p. ej., electrodos frontales secos). |
| **Histograma de amplitud por canal** | Permite comprobar la amplitud típica por canal y detectar saturación, *clipping* o ruido excesivo. Una distribución simétrica centrada en cero y sin colas largas suele indicar buena calidad. |
| **Topografía de varianza** | Resume la actividad media en el dominio espacial. Los canales con varianza (o desviación estándar) mucho mayor que el resto suelen estar contaminados por ruido o desconectados. |
| **PSD media cruda por canal** | El espectro de potencia sin filtrar revela interferencias eléctricas (picos estrechos a 50/100 Hz) o mal acoplamiento (exceso de potencia en frecuencias muy bajas). Útil para planificar el filtrado. |
| **Amplitud media, varianza o desviación estándar por canal** | Caracteriza el nivel de actividad basal y su variabilidad; permite comparar condiciones (p. ej., ojos abiertos vs. ojos cerrados) y detectar canales ruidosos o planos. |
| **Potencia en bandas absoluta o relativa** | Cuantifica la potencia en bandas canónicas (δ, θ, α, β, γ); caracteriza el perfil espectral y facilita la comparación entre condiciones y la detección de artefactos. |
| **Implementación** | Los pasos anteriores (carga de TSV, extracción temporal y de canales, gráficos temporales crudos y apilados, PSD y potencia por bandas, estadísticas por canal y marcado de canales sospechosos) están implementados en la rutina de Julia `src/IO.jl`. |
"

# ╔═╡ aa0047fa-a8e4-4e59-ae92-4e9ee1344811
md"
**Gráfico de un canal individual**, como ejemplo, en este caso el **canal Cz** (central) suele ser representativo de la actividad **EEG** general.
"

# ╔═╡ a1e97a97-3e0b-4352-8211-6ca4a1daaa7a
begin
    p_cz = Plots.plot(
        tiempo_seg, dict_EEG["Cz"],
        xlabel = "Tiempo (s)",
        ylabel = "Amplitud (µV)",
        title = "EEG - Canal Cz",
        legend = false
    )
    p_cz
end

# ╔═╡ 1f9ebe5a-ad63-4be0-8e92-7aa1ee2b5bea
md"
En la representación de todos los **canales apilados**, cada canal se desplaza verticalmente para facilitar la comparación visual. Esta visualización permite detectar artefactos que afectan a múltiples canales. Una separación vertical entre canales (en µV) mayor aumenta la separación visual pero puede hacer que los detalles
de amplitud sean menos visibles.
"

# ╔═╡ eef886e7-89e9-4fa9-8e5a-819a5396aab2
begin
# Separación vertical entre canales (en µV)
# Un valor mayor aumenta la separación visual pero puede hacer que los detalles
# de amplitud sean menos visibles
sep = 80.0

# Calcular offsets (posiciones Y) para cada canal
# Los canales se apilan de arriba a abajo, con el primer canal en la parte superior
offsets = [(length(channels) - i) * sep for i in 1:length(channels)]

# Gráfico apilado
# Se transpone EEG_matrix para que cada columna sea un canal
# Se suman los offsets para desplazar verticalmente cada señal
pila = Plots.plot(
    tiempo_seg,
    EEG_matrix' .+ offsets',  # transpuesta + offsets para apilar
    xlabel = "Tiempo (s)",
    ylabel = "",
    yticks = (offsets, channels),  # etiquetas Y con nombres de canales
    legend = false,
    grid = false,
    size = (1000, 600)
)

pila
end

# ╔═╡ 1ff72c58-70fb-4f55-aa0f-a883c79f72bc
md"
En el **análisis estadístico** de control de calidad de las señales de los canales **EEG** se calculan estadísticas descriptivas para cada canal que permiten evaluar la calidad de los datos y detectar canales problemáticos:

- **Mean**: valor medio (indica offset DC)
- **Min/Max**: rango de amplitudes
- **Std**: desviación estándar (variabilidad, ruido)
- **RMS**: raíz cuadrada de la media de cuadrados (potencia promedio)
- **Kurtosis**: medida de picos (valores altos indican outliers/artefactos)
- **Skewness**: asimetría de la distribución

De manera que los canales con **kurtosis** muy alta (>6) o desviación estándar muy baja (<5 µV) pueden indicar problemas (artefactos, canales muertos, etc.).
"

# ╔═╡ 86335a85-f9da-4ebf-91b0-84e6212048d3
begin
# Calcular estadísticas para cada canal
# eachrow itera sobre cada fila de la matriz (cada canal)
stats_channels = DataFrame(
    :Channel                 => channels,
    Symbol("Mean (µV)")     => mean.(eachrow(EEG_matrix)),      # Media
    Symbol("Min (µV)")      => minimum.(eachrow(EEG_matrix)),   # Mínimo
    Symbol("Max (µV)")      => maximum.(eachrow(EEG_matrix)),  # Máximo
    Symbol("Std (µV)")      => std.(eachrow(EEG_matrix)),       # Desviación estándar
    Symbol("RMS (µV)")      => rms.(eachrow(EEG_matrix)),       # RMS (root mean square)
    :Kurtosis                => kurtosis.(eachrow(EEG_matrix)), # Curtosis (medida de picos)
    :Skewness                => skewness.(eachrow(EEG_matrix)), # Asimetría
)
end

# ╔═╡ 009b77aa-a435-47d6-8c01-c547427f7155
md"
En un **histograma de la desviación estándar** **Std(μV)** por canales, se puede interpretar la variabilidad de los mismos. Valores muy bajos **(<5 µV)** pueden indicar canales muertos o con muy poca señal.
"

# ╔═╡ 6343e583-9f06-42bf-974a-ee722805c7f5
histogram(
    stats_channels[!, Symbol("Std (µV)")],
    bins   = 20,
    xlabel = "Std (µV)",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de desviación típica por canal"
)  

# ╔═╡ 51ea38ab-dc1d-4cad-bacd-1922f7dddbd2
md"
En un **histograma de kurtosis** por canales, se puede ver qué tan **picuda** es la distribución. Como se ha dicho antes, valores altos **(>6)** pueden indicar presencia de artefactos (picos, spikes) o canales problemáticos.
"

# ╔═╡ fb1e4bb0-51db-41d6-8086-dc625b2e45d6
histogram(
    stats_channels.Kurtosis,
    bins = 20,
    xlabel = "Kurtosis",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de kurtosis por canal"
)

# ╔═╡ fb669a4f-0320-4f2c-9589-d91484eac516
begin
# Extraer variables para el scatter plot
stds   = stats_channels[!, Symbol("Std (µV)")]
kurts  = stats_channels.Kurtosis
names_ = stats_channels.Channel

# Criterio de canales sospechosos:
# - Kurtosis > 6: distribución muy picuda (posibles artefactos)
# - Std < 5 µV: variabilidad muy baja (posible canal muerto)
idx_raros = (kurts .> 6) .| (stds .< 5)

# Colores: rojo para sospechosos, azul para normales
colors = ifelse.(idx_raros, :red, :blue)

# Crear scatter plot
scatter_p = Plots.scatter(
    stds, kurts;
    color  = colors,
    xlabel = "Std (µV)",
    ylabel = "Kurtosis",
    title  = "Calidad de canales EEG (Std vs Kurtosis)",
    legend = false,
)

# Añadir etiquetas SOLO a los canales sospechosos
# Se desplazan ligeramente para mejor legibilidad
offsetx = 1.5   # desplazamiento en X
offsety = 0.2   # desplazamiento en Y

for i in eachindex(names_)
    if idx_raros[i]
        annotate!(
            scatter_p,
            stds[i] + offsetx,
            kurts[i] + offsety,
            Plots.text(string(names_[i]), 6, :red)
        )
    end
end

scatter_p

end

# ╔═╡ 7a6de37e-6370-4b16-8da9-f857fc79c01f
md"
## Análisis Espectral

El **análisis espectral** constituye una de las herramientas fundamentales en el estudio de señales **EEG**, ya que permite descomponer la actividad cerebral en sus distintas componentes de frecuencia. A diferencia del análisis en el dominio temporal, el enfoque espectral facilita la identificación de patrones rítmicos que están directamente relacionados con **estados cognitivos**, **procesos fisiológicos** y posibles **alteraciones patológicas**.

En este contexto, la estimación de la **densidad de potencia espectral (PSD)** proporciona una medida cuantitativa de cómo se distribuye la **energía de la señal** a lo largo de las **frecuencias**. Esto resulta clave para caracterizar la dinámica cerebral, comparar condiciones experimentales y evaluar cambios en la actividad neuronal de forma objetiva.

El **análisis espectral** es ampliamente utilizado en múltiples aplicaciones, como el estudio de estados de vigilia y sueño, la detección de anomalías (por ejemplo, actividad epiléptica), el análisis de ritmos asociados a funciones cognitivas (atención, memoria) y el desarrollo de interfaces cerebro-computador (BCI). Además, el **cálculo de potencia en bandas específicas** permite resumir la información de la señal en métricas interpretables y comparables entre sujetos o condiciones.

En conjunto, estas herramientas convierten al **análisis espectral** en un componente esencial dentro de cualquier pipeline de procesamiento de **EEG**, proporcionando una base sólida para la interpretación funcional de la actividad cerebral.
"

# ╔═╡ 424c9b90-8d93-407a-9fc7-840e2c1eec52
md"
### Bandas espectrales

Los análisis espectrales y de conectividad se realizaron dentro de **siete bandas de frecuencia predefinidas**. Estas bandas se seleccionaron en función de las convenciones neurofisiológicas estándar y de su relevancia para la dinámica del **EEG en estado de reposo**:

| Banda | Rango de frecuencia | Interpretación neurofisiológica |
|------|---------------------|--------------------------------|
| $\delta$ | $0.5$–$4\ \text{Hz}$ | Actividad cortical lenta asociada al sueño profundo; incrementos anómalos en vigilia pueden indicar disfunción cortical difusa. |
| $\theta$ | $4$–$8\ \text{Hz}$ | Relacionada con somnolencia, estados de baja vigilancia y procesos de memoria. |
| $\alpha$ | $8$–$12\ \text{Hz}$ | Ritmo dominante en regiones occipitales durante reposo con **ojos cerrados**; marcador clásico del estado basal del EEG. |
| $\beta_{\text{low}}$ | $12$–$15\ \text{Hz}$ | Asociada a actividad sensoriomotora y a estados de alerta moderada. |
| $\beta_{\text{mid}}$ | $15$–$18\ \text{Hz}$ | Relacionada con procesos atencionales y actividad cortical sostenida. |
| $\beta_{\text{high}}$ | $18$–$30\ \text{Hz}$ | Vinculada a procesamiento cognitivo activo y actividad motora; puede aumentar con tensión muscular. |
| $\gamma$ | $30$–$50\ \text{Hz}$ | Oscilaciones rápidas asociadas al procesamiento cortical de orden superior e integración neuronal; susceptible a artefactos musculares. |

La **banda δ** refleja actividad cortical lenta y puede indicar disfunción difusa cuando aparece anormalmente elevada.  
El rango **θ** se asocia con somnolencia y con procesos relacionados con la memoria.  
El **ritmo α** domina en regiones posteriores durante el reposo con **ojos cerrados** y constituye un marcador clásico del estado de reposo.

Las **subbandas β** capturan procesos sensoriomotores y atencionales a frecuencias progresivamente más altas.  
Por último, la **banda γ** refleja actividad oscilatoria rápida vinculada al procesamiento cortical de orden superior, aunque es más susceptible a artefactos musculares.
"

# ╔═╡ 94f78fd4-f8f8-4dfb-a035-4e2b9b227602
begin
# Límites de cada banda en Hz (frecuencia mínima, frecuencia máxima)
Δ = (0.5, 4.0)      # Delta
δ = (4.0, 8.0)      # Theta
α = (8.0, 12.0)     # Alpha
β_low = (12.0, 15.0)      # Beta bajo
β_medium = (15.0, 18.0)   # Beta medio
β_high = (18.0, 30.0)     # Beta alto
γ = (30.0, 50.0)    # Gamma

# Diccionario de límites de bandas
band_limits = Dict(
    :Δ => Δ,
    :δ => δ,
    :α => α,
    :β_low => β_low,
    :β_medium => β_medium,
    :β_high => β_high,
    :γ => γ,
)

# Orden de bandas para visualización y tablas
ordered_bands = [:Δ, :δ, :α, :β_low, :β_medium, :β_high, :γ]

# Etiquetas legibles para cada banda
band_labels = Dict(
    :Δ => "DELTA (0.5-4)",
    :δ => "THETA (4-8)",
    :α => "ALPHA (8-12)",
    :β_low => "BETA LOW (12-15)",
    :β_medium => "BETA MID (15-18)",
    :β_high => "BETA HIGH (18-30)",
    :γ => "GAMMA (30-50)",
)

# Colores para visualización de bandas en gráficos
band_colors = [:lightcyan, :lavender, :lightgoldenrod, :lightgreen, :lightsalmon, :lightpink, :lightgray]
end

# ╔═╡ 042b4a34-d0ed-42f2-9595-f533807b7218
md"
### Densidad de Potencia Espectral (PSD)

Se calcula la **densidad de potencia espectral (PSD)** usando el **método de Welch**, que divide la señal en segmentos superpuestos, aplica ventanas (**Hamming**) y promedia los periodogramas.

La **PSD** muestra cómo se distribuye la potencia de la señal en el dominio frecuencial, lo cual es esencial para:

- Identificar bandas de frecuencia dominantes (α, β, etc.)
- Detectar ruido de línea (50/60 Hz)
- Evaluar la calidad espectral de cada canal
- Comparar características espectrales entre canales

Los parámetros fundamentales para calcular la **PSD** son los siguientes:

- **fs**: frecuencia de muestreo (necesaria para escalar frecuencias)
- **window**: tipo de ventana (**Hamming** reduce leakage espectral)
- **nfft**: tamaño de **FFT** (usar n_muestras para máxima resolución)
"

# ╔═╡ 8a99df75-275d-48ab-be48-1aee6194b9ad
begin
# Calcular PSD para cada canal usando método de Welch
# welch_pgram: divide la señal en segmentos, aplica ventana Hamming,
# calcula periodograma de cada segmento y promedia

PSD = Dict(channel => begin
    p = welch_pgram(dict_EEG[channel]; fs=fs, window=hamming, nfft=n_muestras)
    # Extraer frecuencias y potencias del periodograma
    (; freq = DSP.freq(p), power = DSP.power(p))
end for channel in channels)
	
display(PSD)
end

# ╔═╡ a31e013d-243e-430b-9d5b-ceaeb84eb378
md"
### Visualización de la PSD

Se superponen todas las **PSD**, canales superpuestos, para comparar características espectrales entre canales. La escala logarítmica facilita visualizar rangos amplios de potencia.
"

# ╔═╡ 1f95dd62-aa74-4cf4-8096-9ad6faa9dcdb
begin
# Calcular límites Y amigables para escala logarítmica
# Se buscan potencias de 10 que enmarquen todos los datos
all_powers = vcat([PSD[ch].power for ch in channels]...)
power_min_all, power_max_all = extrema(all_powers)
y_min_log_all = floor(log10(power_min_all))
y_max_log_all = ceil(log10(power_max_all))
ylim_log_all = (10.0^y_min_log_all, 10.0^y_max_log_all)

# Gráfico de PSD superpuestas
p_psd = Plots.plot(
    xlabel = "Frecuencia (Hz)",
    ylabel = "Potencia (µV²/Hz)",
    title = "PSD por canal (Welch)",
    legend = :right,
    xlim = (0, fs / 2),  # Frecuencia de Nyquist = fs/2
    yscale = :log10,     # Escala logarítmica para mejor visualización
    ylim = ylim_log_all,
)

# Añadir cada canal al gráfico
for channel in channels
    Plots.plot!(p_psd, PSD[channel].freq, PSD[channel].power; label = channel, lw = 1)
end

p_psd
end

# ╔═╡ 8ffd90e5-8aa9-4bca-a68a-c7bb88c3b349
md"
### Cálculo de potencia por banda

El **cálculo de la potencia integrada por banda** se realiza a partir del cálculo de la potencia total en cada banda de frecuencia para cada canal.
La potencia se integra sumando las densidades espectrales dentro de cada banda y multiplicando por la **resolución espectral (df)**.

Esto permite cuantificar:

- Qué bandas dominan en cada canal
- Distribución relativa de potencia entre bandas
- Comparaciones entre canales
"

# ╔═╡ 05e66c38-32c9-4e5f-81bd-f066e0e9fcd9
begin
# Calcular potencia integrada por banda para cada canal
# La integración se hace sumando las densidades de potencia en cada banda
# y multiplicando por la resolución espectral (df = diferencia entre frecuencias)
band_powers = Dict(channel => begin
    freqs = PSD[channel].freq
    powers = PSD[channel].power
    
    # Resolución espectral: diferencia entre frecuencias adyacentes
    # Esto es necesario para convertir densidad (µV²/Hz) en potencia total (µV²)
    df = freqs[2] - freqs[1]
    
    # Para cada banda, sumar las potencias dentro del rango de frecuencias
    Dict(
        band => begin
            # Índices de frecuencias dentro de la banda [f_min, f_max)
            idx = (freqs .>= limits[1]) .& (freqs .< limits[2])
            # Potencia integrada = suma de densidades × resolución
            sum(powers[idx]) * df
        end
        for (band, limits) in band_limits
    )
end for channel in channels)
end

# ╔═╡ 9f7c42a7-dc71-4db5-9b51-07c38cea79ba
begin
# Crear DataFrame con potencias por banda y canal
band_powers_df = DataFrame(Band = [band_labels[band] for band in ordered_bands])
channel_totals = Dict{String, Float64}()
# Añadir columnas para cada canal (potencia absoluta y porcentaje)
for channel in channels
    values = [band_powers[channel][band] for band in ordered_bands]
    total_power = sum(values)
    channel_totals[channel] = total_power
    
    # Potencia absoluta en µV²
    band_powers_df[!, Symbol(channel)] = values
    
    # Potencia relativa (porcentaje del total)
    # Evita división por cero si total_power es muy pequeño
    band_powers_df[!, Symbol(channel * "_pct")] = total_power ≈ 0 ? zeros(length(values)) : (values ./ total_power)
end

# Añadir fila de totales
total_row_pairs = Any[:Band => "TOTAL POWER"]
for channel in channels
    push!(total_row_pairs, Symbol(channel) => channel_totals[channel])
    push!(total_row_pairs, Symbol(channel * "_pct") => 1.0)
end
total_row = (; total_row_pairs...)

# Combinar fila de totales con el resto del DataFrame
band_powers_df = vcat(DataFrame([total_row]), band_powers_df; cols = :union)

display(band_powers_df)
end

# ╔═╡ 800c5880-f5de-45bd-8458-3666931ec0ff
md"
## Pipeline fase 0 (Raw data)
"

# ╔═╡ ac015a68-1c14-46e8-9ef2-b00962add288
md"""
!!! info "Pipeline Fase 0: Ingesta de datos y estandarización de metadatos"

    **Objetivo:**
    - Cargar datos **EEG crudos** desde archivos `.tsv`.
    - Validar metadatos y la consistencia dimensional entre sujetos y condiciones.
    - Construir una representación estructurada por sujeto (matriz y diccionario opcional indexado por canal).
    - Generar visualizaciones exploratorias y estadísticas de calidad por canal para inspección previa al preprocesamiento.

    **Entrada:**
    - Directorio de datos (`data/raw/`) con archivos `(.tsv, .json)` por sujeto y condición.  
      Formato compatible con **BrainVision/EEGLAB**; primera columna `Channel`; columnas restantes = muestras temporales (unidades: µV).
    - Registro de sujetos `participants.csv` (opcional, para ejecuciones por lotes).
    - Parámetros fijos de adquisición: ``f_s = 500`` Hz, ``T = 100`` s, ``C = 32`` canales  (``N = F_s \cdot T = 50000`` muestras por canal).

    **Salida:**
    - Estructura por sujeto con metadatos validados (ID, grupo, edad, sexo),  
      matrices **EEG** específicas por condición ``X_{s,k} \in \mathbb{R}^{C \times N}``, nombres de canales, ``f_s``, ``T``  
      y rutas de salida (por ejemplo `results/connectivity/id/`).
    - Opcional: binario serializado (p. ej. `dict_EEG.bin`) para recarga rápida; metadatos `.json` para trazabilidad.
    - Salidas exploratorias:  
      - gráficas temporales (canal individual y canales apilados),  
      - PSD por canal y PSD promedio con sombreado por bandas,  
      - potencia integrada por banda por canal,  
      - estadísticas de calidad de canal (media, desviación estándar, mínimo, máximo, RMS, curtosis, asimetría)  
        y marcado de canales sospechosos (p. ej. kurtosis > 6, desviación estándar < 5 µV).

    **Pasos de procesamiento (rutina `src/IO.jl`).**

    1. **Carga de datos crudos.**  
       - Leer el archivo `.tsv` (p. ej. `CSV.read`) en una tabla; la primera columna 		contiene los nombres de los canales y las columnas restantes corresponden a 		``N`` muestras.  
       - Construir la lista de canales y la matriz **EEG** ``\in \mathbb{R}^{C 				\times N}``; opcionalmente crear un diccionario `canal ↦ señal` para uso 			posterior.

    2. **Metadatos temporales y de canal.**  
       - Definir ``f_s`` y calcular ``N = columnas - 1`` y ``T = N/f_s``.  
       - Validar dimensiones y etiquetas de canal respecto al número esperado de 			canales y la convención de nombres.

    3. **Persistencia opcional.**  
       - Serializar la matriz o el diccionario de canales a un archivo binario  
       (por ejemplo mediante `Serialization.serialize`) en `data/IO/`  
       o en una ruta específica del sujeto para la Fase 1.

    4. **Visualizaciones exploratorias.**  
       - Dominio temporal: señal de un canal (p. ej. Cz) y vista multicanal apilada. 
       - Dominio espectral: PSD por canal (Welch, ventana de Hamming).  
       - PSD promedio con regiones de bandas ``β`` sombreadas.  
       - Opcional: tabla de potencia integrada por banda (absoluta y relativa).

    5. **Calidad de canales.**  
       - Calcular estadísticas por canal: media, mínimo, máximo, desviación 				estándar, RMS, curtosis y asimetría.  
       - Generar histogramas (p. ej. curtosis, desviación estándar) y diagramas de 			dispersión (desviación estándar vs kurtosis).  
       - Marcar como sospechosos los canales con curtosis > 6 o desviación estándar 		< 5 µV (posibles artefactos o mal contacto).

    6. **Agregación a nivel de sujeto (procesamiento por lotes).**  
       - Si se utiliza una lista de sujetos: para cada ``s \in S`` cargar cada 				condición ``k \in K``, construir ``X_{s,k}``, crear 		`results/connectivity/id/` y escribir metadatos.  
       - Devolver la colección preparada para la Fase 1 (Preprocesamiento y control 		de calidad).
"""

# ╔═╡ 188a16df-19e8-4455-9bf1-15bedb31efea
md"""
# Preprocesamiento EEG

El **preprocesamiento** de las señales **EEG crudas** tiene como objetivo eliminar:

- ruido de línea eléctrica,
- interferencias relacionadas con el hardware,
- deriva lenta del baseline,

y limitar el ancho de banda antes del análisis de conectividad específico por banda.

La implementación se proporciona en `src/filtering.jl`, que:

1. carga los datos serializados desde `data/IO/dict_EEG.bin`,
2. aplica una **cascada de filtros digitales**,
3. guarda en cada paso el resultado en `data/filtering/`,
4. genera gráficos de la respuesta del filtro (magnitud y fase),
5. compara la **PSD media** antes y después del filtrado.

Todos los filtros son del tipo **Butterworth**. Las etapas **pasa-altos** y **pasa-bajos** utilizan `filtfilt` (filtrado de fase cero), de modo que:

- el **orden efectivo es el doble** del orden de diseño
- se preserva la **alineación temporal** de la señal.

La cascada de filtrado se aplica en el siguiente orden:

- Filtro **Notch** (50 Hz): elimina interferencia de la red eléctrica.
- Filtro **Bandreject** (100 Hz): elimina posible interferencia de hardware.
- Filtro **Highpass** (0.5 Hz): elimina deriva (drift) preservando actividad lenta.
- Filtro **Lowpass** (150 Hz): limita alias y ruido de alta frecuencia.

---
"""

# ╔═╡ d38f84c5-f99d-4d21-9113-7b7697b9ac59
md"
### Carga de datos (dict_EEG.bin)

Se carga el diccionario **EEG** desde el paso de IO (`dict_EEG.bin`) y se configuran las rutas de salida para los datos filtrados.
"

# ╔═╡ 2eb0aada-37e4-4b30-87d1-76822d782adc
begin
# Directorio base para datos filtrados
dir_filtering = joinpath(@__DIR__, "..", "data", "filtering")
end

# ╔═╡ 611b3e3b-9931-4d8b-b9d3-3ad5d048f60f
md"
### Filtro Notch 

Filtro **band-stop** ($50\,\mathrm{Hz}$) para eliminar la interferencia de la red eléctrica (Europa).

**Diseño:**

- orden: $4$
- ancho de banda: $1\,\mathrm{Hz}$ alrededor de $50\,\mathrm{Hz}$
---
"

# ╔═╡ 830fd9d7-7916-4f56-bb3c-75102e4a26fd
"""Aplica un filtro Butterworth notch_filter a un vector x"""
function notch_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 5, width::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del notch
    freq_low = freq - width/2
    freq_high = freq + width/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end

# ╔═╡ fd19d3b2-e552-457a-aa8c-dd170cb05d7c
begin
"""Plotea la respuesta en frecuencia de un filtro (magnitud y fase)"""
function plot_filter_response(flt, fs::Real; title::String = "Respuesta del Filtro", n_points::Int = 1024)
    # Frecuencias en rad/muestra donde evaluar (0..π)
    ω = range(0, stop = π, length = n_points)
    # Respuesta compleja del filtro en esas frecuencias
    H = freqresp(flt, ω)
    # Pasamos a Hz reales
    f = ω .* fs ./ (2π)
    # Magnitud en dB y fase (convertida a grados sexagesimales)
    mag_db   = 20 .* log10.(abs.(H))
    fase_rad = angle.(H)
    fase_deg = fase_rad .* (180 / π)
    # Aseguramos que la fase esté en el rango [-180°, 180°]
    fase_deg = mod.(fase_deg .+ 180, 360) .- 180
    # Plot de magnitud
    p_mag = Plots.plot(
        f, mag_db,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Magnitud [dB]",
        title  = "$(title) - Magnitud",
        grid   = true,
        linewidth = 2,
        legend = false,
    )
    # Plot de fase
    p_phase = Plots.plot(
        f, fase_deg,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Fase [degrees]",
        title  = "$(title) - Fase",
        grid   = true,
        linewidth = 2,
        ylim = (-180, 180),
        yticks = -180:90:180,
        legend = false,
    )
    # Combinamos ambos plots
    p_combined = Plots.plot(p_mag, p_phase, layout = (2, 1), size = (800, 600))
    
    return p_combined, p_mag, p_phase
end
end

# ╔═╡ 016a1c47-ed53-4dee-9322-274e599d8e88
begin
"""Guarda un diccionario de señales filtradas en un archivo binario"""
function save_filtered_signals(dict_EEG_filtered::Dict{String, Vector{Float64}}, filename::String, base_dir::String = dir_filtering)
    # Crear el directorio si no existe
    if !isdir(base_dir)
        mkpath(base_dir)
        println("📁 Directorio creado: $(abspath(base_dir))")
    end
    # Construir el path completo
    save_path = joinpath(base_dir, filename)
    
    # Guardar el archivo
    Serialization.serialize(save_path, dict_EEG_filtered)
    println("✓ Diccionario EEG filtrado guardado en: $(abspath(save_path))")
    println()
end
end

# ╔═╡ b1337098-3e65-4b6b-8c88-f5becac7fff3
begin
"""Calcula la PSD promedio de un diccionario de señales EEG (sin plotear)"""
function calculate_PSD_average(dict_EEG_data::Dict, fs::Real)
    channels = collect(keys(dict_EEG_data))
    n_muestras = length(dict_EEG_data[first(channels)])
    
    # Calcular PSD para cada canal
    PSD = Dict(channel => begin
        p = welch_pgram(dict_EEG_data[channel]; fs=fs, window=hamming, nfft=n_muestras)
        (; freq = DSP.freq(p), power = DSP.power(p))
    end for channel in channels)
    
    # Calcular PSD promedio
    freqs_avg = PSD[first(channels)].freq
    avg_power = mapreduce(ch -> PSD[ch].power, +, channels) ./ length(channels)
    
    return freqs_avg, avg_power, PSD
end
end

# ╔═╡ 668144ad-6b12-4c89-8280-0a1ba5ed3609
begin
"""Compara y grafica superpuestas las PSD promedio de dos diccionarios de señales EEG"""
function compare_PSD_averages(dict_EEG_1::Dict, dict_EEG_2::Dict, fs::Real;
                              label_1::String = "Original",
                              label_2::String = "Filtrado",
                              title::String = "Comparación PSD Promedio")
    println("📈 COMPARACIÓN DE PSD PROMEDIO")
    println()
    
    # Calcular PSD promedio de ambos diccionarios
    freqs_1, avg_power_1, _ = calculate_PSD_average(dict_EEG_1, fs)
    freqs_2, avg_power_2, _ = calculate_PSD_average(dict_EEG_2, fs)
    
    # Calcular límites Y para ambos conjuntos de datos
    all_powers = vcat(avg_power_1, avg_power_2)
    power_min, power_max = extrema(all_powers)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)
    
    # Crear plot con ambos promedios superpuestos
    p_comparison = Plots.plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title = title,
        legend = :topright,
        xlim = (0, fs / 2),
        yscale = :log10,
        ylim = ylim_log,
    )
    
    Plots.plot!(p_comparison, freqs_1, avg_power_1; label = label_1, lw = 2, color = :blue)
    Plots.plot!(p_comparison, freqs_2, avg_power_2; label = label_2, lw = 2, color = :red)
    
    display(p_comparison)
    println()
    
    return p_comparison
end
end

# ╔═╡ 8db31d27-4aa6-40bc-8c83-fb70727981fa
begin
println("=" ^ 30)
println("📊 Filtro Notch (50 Hz)")
println("=" ^ 30)

Notch_cutoff = 50
Notch_order = 4
Notch_width = 1.0

println()
println("Frecuencia de corte: $(Notch_cutoff) Hz")
println("Orden del filtro: $(Notch_order)")
println("Ancho de banda: $(Notch_width) Hz")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Notch = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_notch = notch_filter(dict_EEG[first(keys(dict_EEG))], Notch_cutoff, fs, Notch_order, Notch_width)

# Aplicamos el filtro Notch a todos los canales
for (channel, signal) in dict_EEG
    y_notch, _ = notch_filter(signal, Notch_cutoff, fs, Notch_order, Notch_width)
    dict_EEG_Notch[channel] = y_notch
end
end

# ╔═╡ 359e26ad-2713-4aff-873f-bd2e765be8d2
begin
# Ploteamos la respuesta del filtro Notch
p_notch, _, _ = plot_filter_response(flt_notch, fs; title = "Filtro Notch (50 Hz)")
p_notch
end

# ╔═╡ b96bf1b3-2e15-49c6-b7d0-877f41569e8a
begin
# Comparar PSD promedio del original vs Notch filtrado (superpuestos)
p_psd_comparison = compare_PSD_averages(
    dict_EEG, 
    dict_EEG_Notch, 
    fs;
    label_1 = "Original",
    label_2 = "Notch (50 Hz)",
    title = "Comparación PSD Promedio: Original vs Notch"
)
end

# ╔═╡ a55e367e-be5a-43da-8bce-279ebaca8ca6
save_filtered_signals(dict_EEG_Notch, "dict_EEG_Notch.bin")

# ╔═╡ b59da6b5-0b1a-4448-86fa-f5fd8196a5b6
md"

### Rechazo de banda

Filtro **band-stop estrecho** ($100\,\mathrm{Hz}$) para eliminar posibles interferencias relacionadas con el hardware alrededor de $100\,\mathrm{Hz}$ (y múltiplos de $50\,\mathrm{Hz}$) sin atenuar ampliamente las frecuencias altas.

**Diseño:**

- orden: $4$
- ancho de banda: $1\,\mathrm{Hz}$
"

# ╔═╡ 91d53d89-df18-42d5-8943-af8a7d6cb74b
begin
"""Aplica un filtro Butterworth bandreject (bandstop) a un vector x"""
function bandreject_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 4, bandwidth::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del bandreject
    freq_low = freq - bandwidth/2
    freq_high = freq + bandwidth/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end
end

# ╔═╡ bc98c842-7ce3-4a31-bfa0-bec06f3997c3
begin
println("=" ^ 40)
println("📊 Filtro Bandreject (100 Hz)")
println("=" ^ 40)
Bandreject_freq = 100
Bandreject_order = 4
Bandreject_bandwidth = 1.0
println()
println("Frecuencia central: $(Bandreject_freq) Hz")
println("Orden del filtro: $(Bandreject_order)")
println("Ancho de banda: $(Bandreject_bandwidth) Hz")
println()
end

# ╔═╡ 94d21c47-149b-4168-b090-410c45e16f4e
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Bandreject = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_bandreject = bandreject_filter(dict_EEG_Notch[first(keys(dict_EEG_Notch))], Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)

# Aplicamos el filtro Bandreject a todos los canales (usando dict_EEG_Notch como entrada)
println("Aplicando filtro Bandreject a todos los canales...")
for (channel, signal) in dict_EEG_Notch
    y_bandreject, _ = bandreject_filter(signal, Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)
    dict_EEG_Bandreject[channel] = y_bandreject
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Bandreject aplicado a $(length(dict_EEG_Bandreject)) canales")
println()
end

# ╔═╡ 0bbe4118-3e47-4854-98e8-73b79aecf55f
begin
# Ploteamos la respuesta del filtro Bandreject
p_bandreject, _, _ = plot_filter_response(flt_bandreject, fs; title = "Filtro Bandreject (100 Hz)")
p_bandreject
end

# ╔═╡ c780bd60-ed60-4fc0-877a-36a85c328a1f
begin
# Comparar PSD promedio del Notch vs Bandreject filtrado (superpuestos)
p_psd_comparison_br = compare_PSD_averages(
    dict_EEG_Notch, 
    dict_EEG_Bandreject, 
    fs;
    label_1 = "Notch (50 Hz)",
    label_2 = "Notch + Bandreject (100 Hz)",
    title = "Comparación PSD Promedio: Notch vs Notch+Bandreject"
)
end

# ╔═╡ a753bee5-a981-4a42-8dda-2b1e10d0f003
save_filtered_signals(dict_EEG_Bandreject, "dict_EEG_Bandreject.bin")

# ╔═╡ 47511ab1-2099-4385-93b0-32f1772713b8
md"
### Filtro Pasa-Alto

Elimina la deriva lenta preservando actividad de baja frecuencia ($0.5\,\mathrm{Hz}$)

**Diseño:**

- orden: $4$
- orden efectivo: $8$ con `filtfilt`

El corte en $0.5\,\mathrm{Hz}$ se eligió para mantener actividad **infra-lenta** que puede ser relevante en **EM (Esclerosis Múltiple)**, mientras se reduce la deriva del baseline.
"

# ╔═╡ 7b988af1-3aab-464e-b3e6-dc9a07128877
begin
"""Aplica un filtro Butterworth pasa-altos a un vector x con zero phase shift (filtfilt)
   
Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo. Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function highpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Highpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end 
end

# ╔═╡ bd16c06c-16bf-4368-8389-5bb413325649
begin
println("=" ^ 30)
println("📊 Filtro Highpass (0.5 Hz)")
println("=" ^ 30)
Highpass_cutoff = 0.5
Highpass_order_design = 4          # Orden del diseño del filtro
Highpass_order_effective = 8       # Orden efectivo (4 × 2 con filtfilt)
println()
println("Frecuencia de corte: $(Highpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Highpass_order_design)")
println("Orden efectivo (con filtfilt): $(Highpass_order_effective)")
println()
end

# ╔═╡ add5d704-abbd-4f5c-a626-05600a057bd4
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Highpass = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_highpass = highpass_filter(dict_EEG_Bandreject[first(keys(dict_EEG_Bandreject))], Highpass_cutoff, fs, Highpass_order_design)

# Aplicamos el filtro Highpass a todos los canales (usando dict_EEG_Bandreject como entrada)
println("Aplicando filtro Highpass a todos los canales...")
for (channel, signal) in dict_EEG_Bandreject
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_highpass, _ = highpass_filter(signal, Highpass_cutoff, fs, Highpass_order_design)
    dict_EEG_Highpass[channel] = y_highpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Highpass aplicado a $(length(dict_EEG_Highpass)) canales")
println()
end

# ╔═╡ 06f567b9-feda-4381-a032-ffe2176926e6
begin
# Ploteamos la respuesta del filtro Highpass
p_highpass, _, _ = plot_filter_response(flt_highpass, fs; title = "Filtro Highpass (0.5 Hz)")
p_highpass
end

# ╔═╡ 64e75122-d2b9-425d-ad99-0d20c891d910
save_filtered_signals(dict_EEG_Highpass, "dict_EEG_Highpass.bin")

# ╔═╡ 2a43830c-9ffe-4dae-9cd2-e23decfdf230
begin
# Comparar PSD promedio del Bandreject vs Highpass filtrado (superpuestos)
p_psd_comparison_hp = compare_PSD_averages(
    dict_EEG_Bandreject, 
    dict_EEG_Highpass, 
    fs;
    label_1 = "Notch + Bandreject",
    label_2 = "Notch + Bandreject + Highpass (0.5 Hz)",
    title = "Comparación PSD Promedio: Bandreject vs Bandreject+Highpass"
)
end

# ╔═╡ d6a7a4c0-5680-4fb4-bc91-4d57d5d3adc5
md"
### Filtro Paso-Bajo 

Limita el ancho de banda superior ($150\,\mathrm{Hz}$).

**Diseño:**

- orden: $4$
- orden efectivo: $8$ con `filtfilt`

Un corte conservador en $150\,\mathrm{Hz}$ preserva la banda **gamma** para inspecciones posteriores; la señal puede restringirse a bandas más bajas (por ejemplo $<50\,\mathrm{Hz}$) en pasos posteriores si es necesario.
"

# ╔═╡ 4f6f1059-01e1-4213-af37-5fd304474212
begin
"""Aplica un filtro Butterworth pasa-bajos a un vector x con zero phase shift (filtfilt)
   
   Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo.
   Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
   Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function lowpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Lowpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end
end

# ╔═╡ 712f89de-087f-4e88-83de-40d58c886896
begin
println("=" ^ 30)
println("📊 Filtro Lowpass (150 Hz)")
println("=" ^ 30)
Lowpass_cutoff = 150
Lowpass_order_design = 4          # Orden del diseño del filtro
Lowpass_order_effective = 8       # Orden efectivo (4 × 2 con filtfilt)
println()
println("Frecuencia de corte: $(Lowpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Lowpass_order_design)")
println("Orden efectivo (con filtfilt): $(Lowpass_order_effective)")
println()
end

# ╔═╡ c4a28722-e0b9-4726-ba15-2aa8bb28ebea
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Lowpass = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_lowpass = lowpass_filter(dict_EEG_Highpass[first(keys(dict_EEG_Highpass))], Lowpass_cutoff, fs, Lowpass_order_design)

# Aplicamos el filtro Lowpass a todos los canales (usando dict_EEG_Highpass como entrada)
println("Aplicando filtro Lowpass a todos los canales...")
for (channel, signal) in dict_EEG_Highpass
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_lowpass, _ = lowpass_filter(signal, Lowpass_cutoff, fs, Lowpass_order_design)
    dict_EEG_Lowpass[channel] = y_lowpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Lowpass aplicado a $(length(dict_EEG_Lowpass)) canales")
println()
end

# ╔═╡ 07e4fff6-77b0-48c8-93da-e354deb017aa
begin
# Ploteamos la respuesta del filtro Lowpass
p_lowpass, _, _ = plot_filter_response(flt_lowpass, fs; title = "Filtro Lowpass (150 Hz)")
p_lowpass
end

# ╔═╡ f608f7b6-8fad-4289-9909-8d57ce90a070
save_filtered_signals(dict_EEG_Lowpass, "dict_EEG_Lowpass.bin")

# ╔═╡ 2bcd6082-ca7b-4353-9d14-f20bc4aa20b7
begin
# Comparar PSD promedio del Highpass vs Lowpass filtrado (superpuestos)
p_psd_comparison_lp = compare_PSD_averages(
    dict_EEG_Highpass, 
    dict_EEG_Lowpass, 
    fs;
    label_1 = "Notch + Bandreject + Highpass",
    label_2 = "Notch + Bandreject + Highpass + Lowpass (150 Hz)",
    title = "Comparación PSD Promedio: Highpass vs Highpass+Lowpass"
)
end

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

# ╔═╡ af42ccde-b32d-4c24-aacc-7b96e158012d
md"
# Postprocesamiento EEG
"

# ╔═╡ b73dff40-9bf4-439e-b413-011cefa8e4cf
md"
## Análisis de Componentes Independientes (ICA) 

Después del preprocesamiento (**Fase 1**), el **EEG filtrado** se descompone en **componentes estadísticamente independientes** para separar la actividad neuronal de artefactos típicos (oculares, musculares y ruido de línea).

**ICA** asume que la señal multicanal observada es una **mezcla lineal de fuentes latentes** que son mutuamente independientes y no gaussianas. Bajo este modelo, las observaciones pueden escribirse como ``X = A S`` donde ``X \in \mathbb{R}^{C \times N}`` es la matriz de C **canales** y N **muestras temporales**,``A \in \mathbb{R}^{C \times C}`` es la **matriz de mezcla** (pesos espaciales), y ``S \in \mathbb{R}^{C \times N}`` contiene las **señales fuente independientes**.

El objetivo es estimar una **matriz de desmezcla** ``W`` tal que ``S = W X`` recupere las fuentes originales (salvo permutación y signo).

Debido a que la suma de variables aleatorias independientes tiende a ser gaussiana (*teorema central del límite*), **ICA** identifica las fuentes **maximizando la no gaussianidad** de las señales proyectadas.
"

# ╔═╡ 996bfc87-a4f3-494d-a809-d7df30857eb2
md"
---

### Algoritmo FastICA simétrico

La implementación utiliza el algoritmo **FastICA simétrico**. Primero, los datos se **centran** (media cero por canal) y posteriormente se **blanquean** mediante PCA para que los datos blanqueados ``Z`` tengan **covarianza identidad**.

El número de componentes se fija igual al número de canales ``k = C``. Esto conserva todas las dimensiones espaciales y permite que el paso posterior de limpieza asigne varianza a componentes neuronales o artefactuales.

La proporción de varianza explicada por cada componente puede calcularse a partir de los autovalores de **PCA** y utilizarse para control de calidad (por ejemplo, detectar sujetos con varianza anormalmente baja en los primeros componentes).

En el espacio blanqueado, la matriz de desmezcla es **ortogonal**. **FastICA** actualiza iterativamente una matriz ortogonal ``W`` maximizando una función de contraste que aproxima la **negentropía**. En esta implementación se utiliza la no linealidad ``g(u) = \tanh(a u)`` con ``a = 1`` adecuada para fuentes supergaussianas (por ejemplo, artefactos oculares o musculares).

En cada iteración se calcula ``Y = W Z`` se aplica la función \(g\) y su derivada, y se actualiza \(W\) mediante un paso tipo Newton seguido de **decorrelación simétrica**, de modo que

```math
W W^{\top} \approx I
```

La convergencia se verifica mediante

```math
\max_i \left|1 - \left| w_i^{new} \cdot w_i^{old} \right|\right| < \tau
```

por ejemplo con ``\tau = 10^{-7}``.

Finalmente, las salidas del algoritmo son:

- los **componentes independientes**

```math
S
```

- la **matriz de desmezcla total**

```math
W_{\mathrm{total}} = W V_{\mathrm{whit}}
```

- la **matriz de mezcla**

```math
A = W_{\mathrm{total}}^{-1}
```

cuyas columnas representan los **mapas topográficos** de cada componente.

El **Algoritmo FastICA** resume el procedimiento **FastICA** simétrico utilizado en `src/ICA.jl`.
"

# ╔═╡ a58fb4d7-12d7-4f82-a00a-926f6772ea5b
md"""
!!! info "Algoritmo FastICA simétrico (implementación en `ICA.jl`)"

    **Entrada:**  
    ``X ∈ ℝ^{C × N}`` (**EEG filtrado**),  
    `k = C`, `max_iter`, `tol`, `a`, `seed`

    1. Centrado de los datos

       ```julia
       for i = 1 to C
           Xc[i,:] ← X[i,:] − mean(X[i,:])
       end
       ```

    2. Blanqueamiento (**whitening**)

       ```text
       Cx ← (1/N) · Xc · Xcᵀ
       [E, Λ] ← eig(Cx)
       Z ← Λ^{-1/2} · Eᵀ · Xc      (usar k componentes)
       ```

    3. Inicialización de la matriz de separación

       ```julia
       set random seed
       W ← randn(k,k)
       W ← sym_decorrelation(W)
       ```

    4. Iteración principal (**FastICA simétrico**)

       ```julia
       for iter = 1 to max_iter

           Y  ← W · Z

           GY ← tanh(a · Y)
           G' ← a · (1 − tanh(a · Y)^2)

           D  ← diag(mean(G', dims=2))

           W_new ← (1/N) · GY · Zᵀ − D · W
           W_new ← sym_decorrelation(W_new)

           M ← W_new · Wᵀ
           max_diff ← max_i |1 − |M_ii||

           if max_diff < tol
               W ← W_new
               break
           end

           W ← W_new

       end
       ```

    5. Reconstrucción de componentes independientes

       ```julia
       S        ← W · Z
       W_total  ← W · V_whit
       A        ← inv(W_total)
       ```

    **Salida:**  
    - `S` (componentes independientes),  
    - `W_total` (matriz de separación total),  
    - `A` (matriz de mezcla).
"""

# ╔═╡ 5dcac750-2be2-4302-92ae-b01464a4ac2c
begin
"""
fastica(X; n_comp, max_iter, tol, a, seed)

Aplica FastICA simétrico a la matriz X para extraer componentes independientes.

El algoritmo busca maximizar la no-gaussianidad de las componentes extraídas
usando la función tanh como aproximación de la negentropía (medida de no-gaussianidad).

Entrada:
  X         :: Matrix{Float64}     # (n_channels × n_samples), ya filtrado
  n_comp    :: Int                 # nº de componentes ICA a extraer
  max_iter  :: Int                 # máximo nº de iteraciones
  tol       :: Float64             # tolerancia de convergencia
  a         :: Float64             # parámetro de tanh (g(u) = tanh(a*u))
  seed      :: Int                 # semilla para reproducibilidad

Salida:
  S        :: Matrix{Float64}  # (n_comp × n_samples), señales independientes
  W_total  :: Matrix{Float64}  # (n_comp × n_channels), matriz de desmezcla total
  A        :: Matrix{Float64}  # (n_channels × n_comp), matriz de mezcla
"""
function fastica(
    X::AbstractMatrix{<:Real};
    n_comp::Int,
    max_iter::Int = 512,
    tol::Float64 = 1e-7,
    a::Float64 = 1.0,
    seed::Int = 1234,
)

    # --------------------------
    # 0) Preparativos
    # --------------------------
    n_channels, n_samples = size(X)
    @assert 1 ≤ n_comp ≤ n_channels "n_comp debe estar entre 1 y n_channels"

    # --------------------------
    # 1) Centrado por canal
    #    (media cero en cada fila)
    # --------------------------
    # ICA requiere que los datos tengan media cero. Se centra cada canal
    # restando su media, lo que elimina el componente DC (corriente continua).
    Xc = copy(Matrix{Float64}(X))  # copiamos a matriz densa de Float64
    for i in 1:n_channels
        μ = mean(Xc[i, :])
        Xc[i, :] .-= μ
    end

    # --------------------------
    # 2) Blanqueo (PCA whitening)
    #    Obtendremos:
    #    - Z: datos blanqueados (n_comp × n_samples)
    #    - V_whit: matriz de blanqueo (n_comp × n_channels)
    # --------------------------
    # El blanqueo reduce la dimensionalidad y elimina correlaciones de segundo orden.
    # Esto simplifica el problema: en lugar de buscar componentes independientes
    # en el espacio original, las buscamos en el espacio blanqueado donde solo
    # necesitamos una matriz ortogonal W (más fácil de optimizar).
    Z, V_whit = whiten_pca(Xc, n_comp)

    # --------------------------
    # 3) Inicialización de W (k × k)
    #    W trabaja en el espacio blanqueado
    # --------------------------
    # W es la matriz de desmezcla en el espacio blanqueado.
    # Se inicializa aleatoriamente y luego se ortonormaliza para que sus filas
    # sean ortogonales (esto asegura que las componentes extraídas sean ortogonales).
    Random.seed!(seed)    # para reproducibilidad
    W = randn(n_comp, n_comp)  # inicialización aleatoria

    # Ortonormalizar W inicialmente
    # Esto asegura que las filas de W sean ortogonales desde el inicio
    W = sym_decorrelation(W)

    # --------------------------
    # 4) Bucle principal de FastICA
    # --------------------------
    # En cada iteración, se actualiza W para maximizar la no-gaussianidad
    # de las componentes Y = W * Z. La función tanh aproxima la negentropía,
    # que es una medida de no-gaussianidad.
    for iter in 1:max_iter
        # Proyección actual: Y = W * Z
        # Y: (n_comp × n_samples)
        # Y contiene las componentes actuales en el espacio blanqueado
        Y = W * Z

        # No linealidad supergaussiana:
        #   g(u)  = tanh(a*u)  - función de contraste (aproxima negentropía)
        #   g'(u) = a * (1 - tanh(a*u)^2)  - derivada de g
        # La función tanh es adecuada para señales supergaussianas (con picos),
        # que es típico en artefactos como parpadeos o actividad muscular.
        GY = tanh.(a .* Y)                      # g(Y)
        Gp = a .* (1 .- tanh.(a .* Y).^2)       # g'(Y)

        # Medias de g'(Y) sobre las muestras (por fila)
        # Estas medias se usan en la actualización de W
        mean_Gp = mean(Gp, dims = 2)           # (n_comp × 1)
        D = Diagonal(vec(mean_Gp))             # diag de medias

        # Actualización tipo FastICA simétrico:
        #   W_new = (1/N) * GY * Z' - D * W
        # Esta fórmula proviene del método de Newton aplicado a la maximización
        # de la negentropía. El primer término es el gradiente, el segundo
        # es un término de normalización.
        W_new = (GY * Z') / n_samples - D * W

        # Ortonormalización simétrica
        # Asegura que las filas de W_new sean ortogonales, lo cual garantiza
        # que las componentes extraídas sean ortogonales (y por tanto independientes
        # en el espacio blanqueado).
        W_new = sym_decorrelation(W_new)

        # ----------------------
        # 4.1) Criterio de convergencia
        #      max_i |1 - |w_i_new ⋅ w_i_old|| < tol
        # ----------------------
        # Se verifica si las filas de W han cambiado significativamente.
        # Si el producto interno entre w_i_new y w_i_old es cercano a 1 (o -1),
        # significa que la dirección no ha cambiado mucho → convergencia.
        M = W_new * W'                         # productos internos entre filas
        diagM = diag(M)                         # productos internos de cada fila consigo misma
        max_diff = maximum(abs.(1 .- abs.(diagM)))  # diferencia máxima

        println("Iteración $iter  ->  max_diff = $max_diff")

        # Si convergió, actualizamos W y salimos
        if max_diff < tol
            println("✔ Convergencia alcanzada en $iter iteraciones (tol = $tol)")
            W = W_new
            break
        end

        # Preparar siguiente iteración
        W = W_new

        # Si llegamos al final sin romper, informamos
        if iter == max_iter
            println("⚠ Se alcanzó max_iter = $max_iter sin converger (max_diff = $max_diff)")
        end
    end

    # --------------------------
    # 5) Construir salidas
    # --------------------------
    # Una vez convergido el algoritmo, construimos las matrices finales:
    # - S: componentes independientes (señales fuente)
    # - W_total: matriz de desmezcla total (del espacio original al de componentes)
    # - A: matriz de mezcla (del espacio de componentes al original)

    # Componentes independientes en el espacio blanqueado:
    # S = W * Z   (n_comp × n_samples)
    # Cada fila de S es una componente independiente (IC)
    S = W * Z

    # Matriz de desmezcla total:
    #   S = W * Z = W * V_whit * Xc = W_total * Xc
    #   => W_total = W * V_whit
    # W_total transforma directamente los datos centrados Xc a componentes S
    # Tamaño: (n_comp × n_channels)
    W_total = W * V_whit

    # Matriz de mezcla aproximada:
    #   Xc ≈ A * S
    # A es la inversa de W_total y representa cómo se mezclan las componentes
    # para reconstruir los datos originales. Las columnas de A son los mapas
    # topográficos de cada componente (pesos espaciales).
    # Tamaño: (n_channels × n_comp)
    A = inv(W_total)

    return S, W_total, A
end
end

# ╔═╡ 0b2fbbd6-faed-48f8-b35c-f3fe7de7526a
begin
println("=============================")
println("Carga de datos filtrados    ")
println("=============================")

# Directorio base para datos ICA (donde se guardarán los resultados)
dir_ica = joinpath(@__DIR__, "..", "data", "ICA")
path_dict_ica = joinpath(dir_ica, "dict_EEG_ICA.bin")

# Directorio de datos filtrados (entrada para ICA)
path_dict_lowpass = joinpath(dir_filtering, "dict_EEG_Lowpass.bin")

# Cargar diccionario con señales por canal
dict_EEG_filtered = if isfile(path_dict_lowpass)
	Serialization.deserialize(path_dict_lowpass)
else
	@warn "No se encontró dict_EEG_Lowpass.bin; generando dict demo desde data_raw (modo público/export)." path_dict_lowpass
	if @isdefined(data_raw) && nrow(data_raw) > 0 && (size(data_raw, 2) > 1)
		chans = Vector{String}(data_raw.Channel)
		mat = Matrix(data_raw[:, Not(:Channel)])
		Dict(chans[i] => vec(mat[i, :]) for i in 1:length(chans))
	else
		Dict("Cz" => randn(2000), "Pz" => randn(2000), "Fz" => randn(2000), "Oz" => randn(2000))
	end
end

println("✓ Archivo: $(basename(path_dict_lowpass))")
println()
display(dict_EEG_filtered)
println()
end

# ╔═╡ eeac27b1-7fff-46a3-acfa-3a540cc989ca
begin
# Orden estable y reproducible de canales (orden alfabético)
# Esto asegura que la matriz X tenga un orden consistente
channels_ICA   = sort(collect(keys(dict_EEG_filtered)))
n_channels = length(channels_ICA)
first_ch   = first(channels_ICA)
n_samples  = length(dict_EEG_filtered[first_ch])

println("✓ Dimensiones: $n_channels canales × $n_samples muestras")
println("✓ Canales (ordenados): $(join(channels, ", "))")
println()
end

# ╔═╡ 47e5c550-a820-44a4-a1bc-ac05851665f4
begin
# Construir matriz X: canales × muestras
# Cada fila corresponde a un canal, cada columna a una muestra temporal
X = Array{Float64}(undef, n_channels, n_samples)

# Rellenar la matriz X con los datos de cada canal
# Se verifica que todos los canales tengan la misma longitud
for (i, ch) in enumerate(channels_ICA)
    sig = dict_EEG_filtered[ch]
    @assert length(sig) == n_samples "Canal $ch con longitud distinta"
    X[i, :] = sig
end

# Número de componentes ICA a extraer
# En ICA, típicamente se extraen tantos componentes como canales hay
k = n_channels  

println("ICA: $n_channels canales, $n_samples muestras, $k componentes")
end

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

# ╔═╡ 230c3dbb-af27-4147-a972-ab1c45578fdb
md"
### Pipeline fase 2b (ICA cleaning)

La etapa de limpieza (`src/ICA_cleaning.jl`) toma los resultados de la Fase 2a y evalúa cada componente independiente utilizando características espaciales, espectrales y estadísticas.

Las **características espaciales** incluyen `frontal_ratio` y `temporal_ratio`, calculadas a partir de las columnas de la matriz de mezcla **A**.

Las **características espectrales** incluyen:

- `blink_ratio` (0.5–4 Hz frente a 4–40 Hz)
- `emg_ratio` (30–80 Hz frente a 1–30 Hz)
- `line_ratio` (48–52 Hz frente al espectro total)

Las **características estadísticas** incluyen:

- **kurtosis**
- `extreme_frac`, definido como la fracción de muestras que superan ±5σ.

Opcionalmente, la característica `corrEOG` utiliza canales EOG cuando están disponibles.

Todas las características se normalizan mediante **z-score** y se combinan para obtener las puntuaciones:

- `ocular_score`
- `muscle_score`
- `line_score`
- `jump_score`

La **puntuación global de artefacto**, `artifact_score`, se define como el máximo de estas puntuaciones.

Los componentes con

```math
artifact\_score > 1.5
```
"

# ╔═╡ bfb33303-de8d-42ef-9d7e-e6b5e813bee8
md"""
!!! info "Pipeline Fase 2b: Limpieza mediante ICA"

    **Objetivo.**
    - Evaluar los componentes independientes (ICs) mediante características espaciales, espectrales y estadísticas.
    - Etiquetar los componentes artefactuales.
    - Reconstruir el EEG limpio anulando las filas artefactuales en ``S`` y calculando ``X_{\text{clean}} = A S_{\text{clean}}``.
    - Guardar los datos limpios y visualizaciones comparativas opcionales.

    **Entrada.**
    - Resultados de ICA de la Fase 2a: `data/ICA/dict_EEG_ICA.bin` (``S``, ``W_{\text{total}}``, ``A``, canales).
    - Posiciones de electrodos en TSV (por ejemplo `data/electrodes/sub-M05_ses-T2_electrodes.tsv`) para visualizaciones topográficas.
    - Opcional: métricas pre-ICA desde `data/filtering/dict_EEG_Lowpass.bin` para comparación.

    **Salida.**
    - `data/ICA/dict_EEG_ICA_clean.bin`: EEG limpio indexado por canal.
    - `data/ICA/dict_EEG_ICA_full.bin`: metadatos completos (ICs descartados, tipo de artefacto, tabla de evaluación).
    - `results/figures/ICA_cleaning/`: topografía + señal temporal + PSD por IC (`IC_001.png`, ...).
    - `results/tables/ICA_components_evaluation.csv`: características, puntuaciones y etiquetas por componente.

    **Pasos de procesamiento (alineados con `src/ICA_cleaning.jl`).**

    1. **Preliminares (opcional).**  
       Cargar `dict_EEG_Lowpass.bin`; calcular curtosis y energía por canal para comparación posterior.

    2. **Carga de ICA y posiciones de electrodos.**  
       Cargar `dict_EEG_ICA.bin` y las posiciones de electrodos; calcular la malla de topografía  
       (proyección azimutal). Generar `plot_ic_summary` (topografía, señal temporal y PSD)  
       para cada IC; guardar en `results/figures/ICA_cleaning/`.

    3. **Extracción de características.**  
       Para cada IC calcular: `frontal_ratio`, `temporal_ratio`, `blink_ratio`, `emg_ratio`,  
       `line_ratio`, curtosis, `extreme_frac`, `corrEOG` (si hay EOG disponible).  
       Normalizar todas las características mediante *z-score*.

    4. **Puntuaciones y etiquetado.**  
       Calcular `ocular_score`, `muscle_score`, `line_score`, `jump_score`.  
       Definir `artifact_score` como el máximo de estas puntuaciones.  
       Etiquetar como **parpadeo / músculo / línea / salto** si `artifact_score > 1.5`;  
       en caso contrario etiquetar como **neuronal**.  
       Guardar resultados en `results/tables/ICA_components_evaluation.csv`.

    5. **Reconstrucción de la señal.**  
       Establecer ``S_{\text{clean}}[i,:] = 0`` para cada IC artefactual ``i``.  
       Calcular

       ```math
       X_{\text{clean}} = A S_{\text{clean}}
       ```

       Guardar el diccionario por canal en `dict_EEG_ICA_clean.bin`  
       y los metadatos completos en `dict_EEG_ICA_full.bin`.

    6. **Opcional.**  
       Generar gráficos comparativos antes/después de la limpieza  
       (segmentos EEG, PSD, curtosis y energía por canal).
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

# ╔═╡ 19a03e63-b353-4206-9deb-823d52bcd2f3
md"
## Corrección de Línea Base (antes de artefactos)

Tras la segmentación en épocas de longitud fija, cada segmento puede conservar un desplazamiento DC (offset) o una deriva lenta que difiere entre segmentos. La corrección de baseline elimina este desplazamiento restando, de cada segmento, la media de un intervalo de referencia predefinido (normalmente el inicio del segmento). El resultado es que la ventana de baseline tiene media aproximadamente cero en cada segmento, de modo que las estimaciones de amplitud y espectrales no se vean sesgadas por offsets específicos de cada segmento.

En EEG en reposo sin estímulos discretos, este paso no corresponde a un baseline fisiológico de tipo pre-estímulo en el sentido de potenciales evocados; actúa como una normalización por segmento que elimina el nivel DC de la porción inicial de cada época y homogeneiza los segmentos para análisis posteriores de conectividad o espectro. Una alternativa sería restar la media de todo el segmento (demean de toda la época), lo cual también es válido para análisis puramente espectrales; sin embargo, aquí se mantiene el enfoque basado en un intervalo para conservar la consistencia con pipelines previos y permitir una ventana de baseline bien definida (por ejemplo, los primeros 100 ms) para informes y comprobaciones de calidad.

La implementación (`src/baseline.jl`) asume datos segmentados almacenados como un arreglo tridimensional

```math
\mathbf{E} \in \mathbb{R}^{C \times L \times K}
```

(canales × muestras por segmento × segmentos).

La ventana de baseline se define como

```math
[t_{\mathrm{start}}, t_{\mathrm{end}}]
```

en segundos. A una frecuencia de muestreo

```math
F_s
```

esto corresponde a los índices de muestra

```math
n_{\mathrm{start}} \text{ -- } n_{\mathrm{end}}
```

Para cada canal

```math
c
```

y segmento

```math
k
```

la media del baseline es

```math
\mu_{c,k} = \frac{1}{N_{\mathrm{bl}}} \sum_{n=n_{\mathrm{start}}}^{n_{\mathrm{end}}} E_{c,n,k}
```

donde

```math
N_{\mathrm{bl}} = n_{\mathrm{end}} - n_{\mathrm{start}} + 1
```

El segmento corregido se define como

```math
E_{c,:,k}^{\mathrm{corr}} = E_{c,:,k} - \mu_{c,k}
```

De este modo, el intervalo de baseline tiene media cero y el resto del segmento se desplaza por la misma constante. La desviación estándar del segmento permanece inalterada.

La configuración por defecto utiliza

```math
t_{\mathrm{start}} = 0
```

y

```math
t_{\mathrm{end}} = 0.1\,\text{s}
```

(primeros 100 ms). Con

```math
F_s = 500\,\text{Hz}
```

esto corresponde a las muestras 1–50.

La ventana de 100 ms es lo suficientemente corta para evitar la contaminación por deriva lenta dentro del segmento y lo suficientemente larga para proporcionar una estimación estable de la media del segmento utilizada en la sustracción del baseline.

El Algoritmo 2 resume el procedimiento; la Tabla 2 y el cuadro inferior formalizan la Fase 4 del pipeline.
"

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

# ╔═╡ 3713b98b-ae0a-440f-b3a3-3bfebe68a937
md"
## Rechazo por Artefactos

Después de la corrección de baseline, las épocas segmentadas pueden seguir conteniendo artefactos residuales como breves ráfagas musculares, saltos de electrodo o saturaciones que sobreviven a la limpieza basada en ICA. Para evitar contaminar las estimaciones de conectividad y espectro, se aplica un paso final de rechazo de artefactos que descarta segmentos completos cuya amplitud supera umbrales conservadores.

La rutina `src/artifact_rejection.jl` implementa un criterio simple e interpretable basado en amplitud: solo se inspeccionan los primeros

```math
N_{\mathrm{used}}
```

canales (por defecto: 30), y cualquier segmento en el que al menos una muestra en estos canales quede fuera del rango

```math
[A_{\min}, A_{\max}]
```

(por defecto

```math
A_{\min} = -70\,\mu\text{V}, \quad A_{\max} = +70\,\mu\text{V}
```

) se marca como inválido y se elimina.

Esto produce una hipermatriz con menos segmentos pero con las mismas dimensiones de canales y tiempo. Los umbrales de

```math
\pm 70\,\mu\text{V}
```

se eligen para detectar artefactos evidentes (saturaciones, ráfagas EMG extremas) manteniéndose por encima de las amplitudes típicas en EEG en reposo; representan un compromiso conservador entre sensibilidad a valores atípicos y conservación de datos utilizables.

Los parámetros `before_event_ms` y `after_event_ms` se almacenan para un posible rechazo relacionado con eventos, aunque la implementación actual utiliza únicamente umbrales globales de amplitud sobre todo el segmento.

Formalmente, sea

```math
\mathbf{E}^{\mathrm{bl}} \in \mathbb{R}^{C \times L \times K}
```

el EEG segmentado y corregido por baseline (canales × muestras × segmentos), y sea

```math
\mathcal{C}_{\mathrm{used}} = \{1,\dots,N_{\mathrm{used}}\}
```

el conjunto de índices de canales que se inspeccionan.

Para el segmento

```math
k
```

los datos inspeccionados son

```math
\mathbf{E}^{\mathrm{bl}}[\mathcal{C}_{\mathrm{used}}, :, k]
```

Se definen

```math
m_k = \min_{c \in \mathcal{C}_{\mathrm{used}},\, 1 \le n \le L} E^{\mathrm{bl}}_{c,n,k}
```

y

```math
M_k = \max_{c \in \mathcal{C}_{\mathrm{used}},\, 1 \le n \le L} E^{\mathrm{bl}}_{c,n,k}
```

El segmento se marca como inválido si

```math
m_k < A_{\min} \quad \text{o} \quad M_k > A_{\max}
```

En caso contrario, se conserva.

La tasa de rechazo se define como

```math
\mathrm{rate} = \frac{K - K_{\mathrm{AR}}}{K}
```

donde

```math
K_{\mathrm{AR}}
```

es el número de segmentos retenidos. Esta magnitud sirve como métrica de calidad (por ejemplo, reportada en `artifact_rejection_statistics.csv`).

Con fines diagnósticos, también se registra el conjunto de canales que violan los umbrales:

```math
\mathcal{V}_k =
\{c \in \mathcal{C}_{\mathrm{used}} :
\min_n E^{\mathrm{bl}}_{c,n,k} < A_{\min}
\ \text{o} \
\max_n E^{\mathrm{bl}}_{c,n,k} > A_{\max}\}
```

La hipermatriz filtrada

```math
\mathbf{E}^{\mathrm{AR}}
```

se construye concatenando únicamente los segmentos válidos.

El Algoritmo 3 resume el procedimiento; la Tabla 3 y el cuadro inferior formalizan la Fase 5 del pipeline.
"

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

# ╔═╡ 1c9a3df0-545c-4695-aef4-face8703449d
md"
## Corrección de baseline (después de artefactos)

Una segunda corrección de baseline opcional puede aplicarse a los segmentos que han pasado el rechazo de artefactos (salida de la Fase 5) si el análisis lo requiere; este paso se formaliza como la Fase 6 del pipeline.

Si se utiliza, se aplica el mismo procedimiento basado en segmentos que en la Fase 4 (restar, para cada segmento, la media de una ventana de baseline) al conjunto reducido de segmentos.

El resultado se utiliza posteriormente como entrada para la estimación de conectividad (Fase 7).
"

# ╔═╡ 20569852-d89f-49f0-8d73-30ad7d2261d9
md"""
### Pipeline fase 6 (Baseline after artifacts)

Idem fase 4
"""

# ╔═╡ 3cac8ebb-b5ab-4710-8a92-7627813b579e
md"
## Análisis Espectral

Tras la segunda corrección de baseline opcional (Fase 6), el EEG segmentado se utiliza para estimar espectros de potencia en el dominio de la frecuencia.

La rutina `src/FFT.jl` implementa un pipeline espectral compatible con BrainVision Analyzer (BVA): carga los datos desde `data/baseline/dict_2nd_baseline_correction.bin` (canales × muestras × segmentos), elimina el desplazamiento DC por segmento y canal (estilo BVA, antes del ventaneo), aplica una atenuación basada en ventana de Hamming (Tukey con longitud de ventana del 10 %), y posteriormente aplica zero-padding hasta

```math
N_{\mathrm{fft}} = 512
```

puntos para obtener una resolución en frecuencia

```math
\Delta f = \frac{F_s}{N_{\mathrm{fft}}} \approx 0.976\,\text{Hz}
```

Se calcula la FFT real (rFFT); el componente DC se reintroduce en el bin de

```math
0\,\text{Hz}
```

La corrección de varianza (división por

```math
\overline{w^2}
```

, la media del cuadrado de la ventana) y el plegado del espectro completo (“Use Full Spectrum”: duplicar los bins interiores y normalizar por

```math
N_{\mathrm{fft}}^2
```

) producen la potencia por bin en

```math
\mu\text{V}^2
```

La potencia específica por banda se obtiene promediando los bins cuya frecuencia central cae dentro de cada banda (Delta, Theta, Alpha, Beta Low/Mid/High, Gamma).

Los resultados se promedian entre segmentos para cada canal y se guardan para su uso en gráficos y en los flujos de trabajo de conectividad.

El Algoritmo 4, la Tabla 4, la Tabla 5 y el cuadro inferior resumen este paso del pipeline.
"

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

# ╔═╡ f8bd3fd9-93c3-4b16-8b68-5c780a013e9b
md"""
### Pipeline fase 7 (Spectral Analysis)
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
md"
# Conectividad cerebral

La **conectividad cerebral** se refiere a los patrones de interacción o acoplamiento entre distintas poblaciones neuronales. 

Generalmente se clasifica en tres niveles conceptuales:

**Conectividad estructural (CE)**, describe el cableado físico del cerebro —vías axonales y sinapsis— revelado mediante técnicas como la imagen por difusión o la histología. Esta conectividad determina qué regiones pueden interactuar, pero no especifica el acoplamiento dinámico entre ellas.

**Conectividad funcional (CF)**, cuantifica dependencias estadísticas entre señales registradas (por ejemplo, correlación, coherencia o sincronización de fase) sin asumir una dirección causal. Permite determinar si dos regiones covarían o oscilan de forma coordinada.

**Conectividad efectiva (CE)**, busca inferir influencias dirigidas o causales (por ejemplo, qué región impulsa a otra), normalmente mediante enfoques basados en modelos como la causalidad de Granger o el modelado causal dinámico. Este tipo de análisis requiere más datos y supuestos que la conectividad funcional.

En **EEG** y **MEG**, la conectividad se evalúa habitualmente a nivel **funcional**, estimando cuán fuerte y consistente es la relación entre señales de sensores o fuentes en el tiempo o en la frecuencia. 

Los registros en **estado de reposo (*resting-state*)** son especialmente adecuados para caracterizar este acoplamiento en ausencia de estructura de tarea. La conectividad funcional se ha utilizado ampliamente para estudiar la organización de redes cerebrales tanto en condiciones saludables como patológicas.

En **Esclerosis Múltiple**, las alteraciones de la sustancia blanca y gris pueden modificar tanto la conectividad estructural como la funcional. Las medidas de conectividad funcional basadas en EEG pueden capturar reorganizaciones o debilitamientos en la coordinación oscilatoria de larga distancia, incluso cuando el daño estructural aún no es evidente o se encuentra distribuido.

Debido a que el **EEG** del cuero cabelludo está afectado por la **conducción de volumen**, las métricas de conectividad sensibles a correlaciones de desfase cero o cercano a cero (por ejemplo, correlación o coherencia) pueden verse confundidas por fuentes compartidas. Por esta razón, suelen preferirse métricas basadas en fase que enfatizan relaciones de fase consistentes y reducen la contribución de interacciones de desfase cero en el espacio de sensores.

---

## Conectividad funcional vs. conectividad efectiva

La **conectividad funcional (FC)** describe dependencias estadísticas entre eventos neurofisiológicos remotos sin inferir direccionalidad. No se basa en un modelo generativo explícito del acoplamiento y, normalmente, implica la comparación con una hipótesis nula de independencia.

Por el contrario, la **conectividad efectiva (EC)** modela influencias causales dirigidas entre regiones cerebrales. Este enfoque evalúa hipótesis sobre la arquitectura de acoplamiento y depende de un modelo explícito de las interacciones entre regiones.

---

## Estimación de conectividad funcional estática

Esta sección describe cómo estimamos **conectividad funcional estática** a partir de EEG en estado de reposo. Para cada sujeto y condición se obtiene una matriz de conectividad por banda de frecuencia, que resume la intensidad del acoplamiento entre todos los pares de sensores a lo largo de la grabación.

El pipeline consta de **cuatro etapas**, cada una diseñada para abordar un problema metodológico específico.

### Justificación del orden del pipeline

Los segmentos preprocesados (después de la corrección de baseline y del rechazo de artefactos) aún reflejan la **mezcla de fuentes causada por la conducción de volumen**. Cada canal EEG registra una combinación de múltiples generadores neuronales, por lo que medidas directas como la correlación o la coherencia entre canales pueden verse fuertemente influenciadas por la propagación espacial del campo eléctrico y sobreestimar las verdaderas interacciones neuronales.

Por esta razón, el pipeline sigue el siguiente orden:

**Primero**, reducimos los efectos de la conducción de volumen transformando los potenciales en **Current Source Density (CSD)**. Esta transformación enfatiza la actividad local y produce señales independientes de la referencia.

**Segundo**, para cada banda de frecuencia de interés se aplican filtros pasa-banda a las señales CSD. A continuación se obtiene su representación analítica mediante la **transformada de Hilbert** y se calcula el **weighted Phase Lag Index (wPLI)** entre todos los pares de canales. El wPLI se centra en relaciones de fase consistentes y reduce el peso de las contribuciones de desfase cero asociadas a conducción de volumen.

Este procedimiento produce **una matriz de conectividad por época**. Posteriormente, estas matrices se agregan entre épocas (por ejemplo mediante la mediana) para obtener una única matriz **estática** por sujeto, condición y banda de frecuencia.

**Tercero**, se construye una **distribución nula** generando datos sustitutos (*surrogates*) mediante aleatorización de fase y recalculando el wPLI para cada conjunto de datos generado. Esto permite asignar un **valor $p$ empírico** a cada conexión y controlar falsos positivos (por ejemplo mediante corrección FDR).

**Cuarto**, las matrices se **umbralizan** y se derivan métricas de red a nivel global, como:

- fuerza global de conectividad  
- fuerza nodal  
- densidad de conexiones significativas  

Finalmente, estas métricas se comparan entre grupos (por ejemplo **pacientes con Esclerosis Múltiple frente a controles**) utilizando **tests estadísticos no paramétricos**.

### Organización de las siguientes secciones

Las subsecciones siguientes describen cada etapa del pipeline:

1. **CSD** para reducir la conducción de volumen  
2. **Estimación de conectividad mediante wPLI** y agregación de matrices estáticas  
3. **Generación de datos surrogate** e inferencia estadística  

En conjunto, este pipeline permite realizar una **comparación técnicamente sólida e interpretable de la conectividad funcional en estado de reposo entre grupos y condiciones experimentales**.
"

# ╔═╡ 2c8d9444-4e8f-4234-b31d-c7e3425f8d13
md"
## Reducción Conducción volumen (CSD)

La **conducción de volumen**, del inglés Current Source Density, **(CSD)** en la cabeza dispersa las corrientes de modo que cada electrodo del cuero cabelludo registra una mezcla ponderada de múltiples fuentes. Esto reduce la resolución espacial e incrementa artificialmente la coherencia entre canales cercanos. Para enfatizar la actividad local y mejorar la interpretabilidad de la conectividad, aplicamos una transformación de **CSD** antes de calcular las medidas de conectividad.

La **CSD** es un **Laplaciano espacial** que aproxima la segunda derivada del potencial con respecto a la superficie del cuero cabelludo, atenuando las contribuciones debidas a la conducción de volumen y haciendo que la señal sea efectivamente independiente de la referencia.

La implementación sigue el método de **Perrin et al. (1989)**: interpolación mediante *spherical splines* sobre la esfera unitaria, expansión en polinomios de Legendre para el núcleo del spline y un operador Laplaciano regularizado equivalente al utilizado en **BrainVision Analyzer**.

La rutina `src/Connectivity/CSD.jl`:

- carga los segmentos ya corregidos por línea base
- carga las posiciones de los electrodos (BIDS `.tsv`)
- verifica la consistencia entre canales y electrodos
- construye el operador de Perrin **L**
- aplica la transformación **CSD** por cada instante temporal y para cada segmento

La transformación se define como

```math
E_{CSD} = L E
```
"

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

# ╔═╡ 04e8aa59-fecb-401f-8adf-95374473dbb7
md"
## Conectividad funcional basada en fase

Las medidas de **conectividad funcional basadas en fase** cuantifican interdependencias estadísticas entre señales oscilatorias centrándose en sus **relaciones de fase instantáneas**.

En **EEG**, este enfoque es especialmente relevante porque la **conducción de volumen** puede inducir correlaciones de **desfase cero** que no reflejan un acoplamiento fisiológico genuino entre regiones cerebrales.

Las medidas derivadas del **espectro cruzado** permiten separar las contribuciones de **magnitud** y **fase**, facilitando el análisis de la sincronización entre señales neuronales.

---

### Formulación en dominio temporal y frecuencial

Sean `x_k(t)` e `y_k(t)` dos señales EEG (a nivel de sensor o de fuente).  
En el dominio temporal, su **correlación cruzada** se define como


```math
R_{xy}(\tau) = \mathbb{E}[x(t)y(t+\tau)]
```

En forma muestral,

```math
\hat{R}_{xy}(\tau) =
\frac{1}{N}\sum_{k=1}^{N} x_k y_{k+\tau}
```

Las correlaciones con retardo temporal (`\tau`) enfatizan interacciones retardadas, mientras que promediar sobre `\tau` elimina información direccional.

Aplicando la **transformada de Fourier** se obtiene la representación en frecuencia:

```math
X(f) = \mathcal{F}\{x(t)\}, \qquad
Y(f) = \mathcal{F}\{y(t)\}
```

El **espectro cruzado** se define como

```math
S_{xy}(f) = X(f)Y^*(f)
```

donde `(\cdot)^*` denota conjugación compleja.  
Los **autoespectros** correspondientes son `S_{xx}(f)` y `S_{yy}(f)`.

---

### Coherencia

La **coherencia** cuantifica el acoplamiento lineal normalizado entre dos señales en el dominio de la frecuencia:

```math
\mathrm{coh}_{xy}(f) =
\frac{|S_{xy}(f)|}{\sqrt{S_{xx}(f)S_{yy}(f)}}
```

Esta normalización restringe la coherencia al intervalo `0 \leq coh \leq 1`, donde:

- **0** indica ausencia de sincronización  
- **1** indica sincronización perfecta  

Sin embargo, debido a que incluye componentes de **fase cero**, la coherencia es sensible a efectos de **conducción de volumen**.

La **fase del espectro cruzado** viene dada por

```math
\theta(f) =
\arctan
\left(
\frac{\Im\{S_{xy}(f)\}}
{\Re\{S_{xy}(f)\}}
\right)
```

lo que refleja la **diferencia de fase** entre las señales a la frecuencia `f`.

---

### Phase Locking Value (PLV)

Para aislar la consistencia de fase entre ensayos o segmentos temporales se utiliza el **Phase Locking Value (PLV)**:

```math
PLV =
\left|
\frac{1}{N}
\sum_{k=1}^{N}
e^{i\Delta\phi_k}
\right|
```

donde `\Delta\phi_k` representa la **diferencia de fase instantánea** entre las señales en el ensayo `k`.

El **PLV** enfatiza relaciones de fase consistentes **independientemente de la amplitud de las señales**.

---

### Amplitude Envelope Correlation (AEC)

La conectividad basada en amplitud puede estimarse mediante la **correlación entre envolventes de amplitud** de las señales analíticas:

```math
AEC = \mathrm{corr}(A_x(t), A_y(t))
```

donde la envolvente de amplitud se obtiene mediante la **transformada de Hilbert**:

```math
A_x(t) = |\mathcal{H}\{x(t)\}|
```
"

# ╔═╡ 47af6c89-b4c4-43e6-a2cf-1771cb45089e
md"
## Weighted Phase Lag Index (wPLI)

La conectividad funcional basada en fase se estima utilizando el **Weighted Phase Lag Index (wPLI)**, que enfatiza relaciones de fase consistentes distintas de cero entre señales y reduce el peso de las contribuciones de **fase cero**, que probablemente se deben a **conducción de volumen**.  

Para cada banda de frecuencia, las señales se filtran mediante **band-pass** y se transforman en representaciones analíticas mediante la **transformada de Hilbert**. Posteriormente, el wPLI se calcula agregando todos los samples de los segmentos (estilo BrainVision Analyzer, *across segments*).

---

**Definición formal.**  
Sea ``S_{xy}(t) = z_i(t)\,\overline{z_j(t)}`` el **cross-spectrum** en el instante ``t`` entre las señales analíticas ``z_i`` y ``z_j`` de los canales ``i`` y ``j``. El wPLI se define como

```math
\mathrm{wPLI}_{ij} =
\frac{\left|\mathbb{E}\left[\mathrm{Im}(S_{xy})\right]\right|}
{\mathbb{E}\left[|\mathrm{Im}(S_{xy})|\right] + \epsilon}
=
\frac{\left|\mathbb{E}\left[\mathrm{sign}(\mathrm{Im}_{ij}(t))\cdot|\mathrm{Im}_{ij}(t)|\right]\right|}
{\mathbb{E}\left[|\mathrm{Im}_{ij}(t)|\right] + \epsilon}
```

donde ``\mathrm{Im}_{ij}(t) = \mathrm{Im}(S_{xy}(t)) = \mathrm{Im}(z_i(t)\overline{z_j(t)})`` y ``\epsilon`` es una pequeña constante para evitar divisiones por cero.  

La **esperanza matemática** se calcula sobre todos los samples temporales y segmentos.

---

**¿Por qué solo la parte imaginaria?**  

La **conducción de volumen** y la **referencia común** tienden a introducir acoplamientos de fase cero (*in-phase*), que contribuyen principalmente a la **parte real** de ``S_{xy}``. En cambio, la **parte imaginaria** refleja desfases consistentes entre señales y, por tanto, está menos sesgada por la propagación espacial de corrientes.

En comparación con el **Phase Lag Index (PLI)** clásico, que utiliza únicamente ``\mathrm{sign}(\mathrm{Im}_{ij})``, el **wPLI** pondera cada muestra por ``|\mathrm{Im}_{ij}|``. Esto reduce la sensibilidad a pequeñas diferencias de fase ruidosas y mejora la robustez en muestras pequeñas. Además, resulta más robusto frente a conducción de volumen que la **coherencia** tradicional, ya que atenúa contribuciones de fase cero.

---

La matriz de conectividad ``\mathbf{W} \in \mathbb{R}^{C \times C}`` es **simétrica** y tiene **diagonal nula**, con valores en el intervalo ``[0,1]``.

Las matrices de conectividad ``(C \times C)``, con ``C = 32`` en este dataset se obtienen para cada banda agregando todos los samples de los segmentos. Las entradas diagonales son cero y la matriz resultante es simétrica.
"

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

# ╔═╡ 06031333-5d1f-4cef-b1c4-c7f8f4093e66
md"
# Resultados
"

# ╔═╡ c586474d-9b33-4f1d-96fb-7f948724f4fe
md"
# Implementación en Julia

## ¿Por qué Julia?

**Julia** es un lenguaje de programación de alto nivel diseñado específicamente para **computación científica y numérica de alto rendimiento**. Fue desarrollado con el objetivo de combinar la facilidad de uso de lenguajes dinámicos con la velocidad de los lenguajes compilados tradicionales.

En proyectos de análisis de datos neurofisiológicos —como el procesamiento de EEG y el análisis de conectividad— es habitual enfrentarse a grandes volúmenes de datos, transformaciones numéricas intensivas y pipelines complejos. En este contexto, Julia ofrece una combinación especialmente adecuada de **legibilidad, rendimiento y capacidad de desarrollo científico**.

### Rendimiento y eficiencia

Una de las características más destacadas de Julia es que su código se compila mediante **Just-In-Time (JIT compilation)** usando LLVM. Esto permite que muchas operaciones numéricas se ejecuten a velocidades comparables a lenguajes compilados como **C o Fortran**, sin sacrificar la simplicidad de un lenguaje de alto nivel.

En lenguajes como Python, alcanzar este nivel de rendimiento suele requerir escribir partes críticas del código en C, Cython o utilizar librerías externas altamente optimizadas. En Julia, por el contrario, es posible escribir directamente código científico de alto rendimiento en el propio lenguaje.

Esto resulta particularmente útil en tareas comunes en análisis EEG como:

- transformadas de Fourier  
- filtrado digital  
- operaciones matriciales de gran tamaño  
- cálculo de métricas de conectividad entre múltiples canales  
- simulaciones o generación de datos surrogate

### Comparación con Python

Python es actualmente uno de los lenguajes más utilizados en ciencia de datos y neurociencia computacional gracias a su gran ecosistema de librerías. Sin embargo, Julia ofrece varias ventajas relevantes:

- **Mayor rendimiento en código numérico puro**, evitando la necesidad de extensiones en C.
- **Tipado múltiple y despacho múltiple**, lo que permite escribir código genérico y altamente optimizado.
- **Sintaxis matemática clara**, especialmente adecuada para expresar algoritmos científicos.

En términos de experiencia de usuario, Julia mantiene muchas de las ventajas de Python:

- sintaxis legible  
- uso interactivo en notebooks  
- ecosistema creciente de paquetes científicos  

pero con un rendimiento significativamente mayor en muchos casos.

### Comparación con MATLAB

MATLAB ha sido tradicionalmente una herramienta muy utilizada en procesamiento de señales y neurociencia. Julia ofrece un entorno muy similar en cuanto a estilo de programación:

- operaciones vectorizadas y matriciales nativas  
- álgebra lineal integrada en el lenguaje  
- sintaxis matemática compacta  

Sin embargo, Julia presenta varias ventajas importantes:

- **software completamente libre y open source**
- **rendimiento comparable a lenguajes compilados**
- **gestión moderna de paquetes y entornos reproducibles**
- **mejor integración con herramientas modernas de computación científica**

En resumen, Julia combina características que tradicionalmente estaban separadas entre distintos lenguajes:  
es **tan fácil de usar como Python**, **tan natural para cálculo matricial como MATLAB**, y **tan rápido como C**.

---

## Pluto.jl como entorno 

El desarrollo, documentación y exploración interactiva del pipeline de análisis EEG se realizaron utilizando **Pluto.jl**, un entorno de notebooks reactivo diseñado específicamente para el lenguaje **Julia**. Pluto combina características de programación interactiva, documentación científica y reproducibilidad computacional en una única interfaz ligera.

A diferencia de los notebooks tradicionales, Pluto está diseñado desde el principio para garantizar **consistencia, trazabilidad y ejecución reproducible**, lo que lo convierte en una herramienta especialmente adecuada para proyectos de investigación científica.

---

### ¿Qué es Pluto.jl?

**Pluto.jl** es un sistema de notebooks interactivos basado en Julia que permite combinar:

- código ejecutable  
- texto explicativo (Markdown)  
- ecuaciones matemáticas  
- visualizaciones y resultados  

en un único documento dinámico.

Cada celda del notebook contiene una única expresión de Julia, y Pluto gestiona automáticamente las dependencias entre celdas. Cuando una variable cambia, **todas las celdas que dependen de ella se recalculan automáticamente**, manteniendo el notebook siempre en un estado consistente.

Este modelo reactivo elimina problemas comunes en otros entornos de notebooks, como el uso de variables obsoletas o ejecuciones fuera de orden.

---

### Ventajas frente a Jupyter Notebooks

Aunque **Jupyter Notebook** es una herramienta ampliamente utilizada en ciencia de datos, Pluto introduce varias mejoras importantes para el trabajo científico reproducible.

#### Ejecución reactiva y coherencia del estado

En Jupyter, las celdas pueden ejecutarse en cualquier orden, lo que puede generar inconsistencias entre el estado real del kernel y el código visible en el notebook. Esto puede dificultar la reproducción exacta de resultados.

Pluto evita este problema mediante un **modelo reactivo basado en dependencias**:

- cada variable solo puede definirse una vez  
- las dependencias entre celdas se detectan automáticamente  
- los cambios en una celda actualizan todas las celdas dependientes  

Esto garantiza que el notebook siempre represente un **estado válido y reproducible del análisis**.

---

## Cómo iniciar Pluto.jl

Para trabajar con Pluto es necesario tener **Julia instalado** y el paquete `Pluto` añadido al entorno.  
Si es la primera vez que se utiliza, puede instalarse desde el gestor de paquetes de Julia:

```
using Pkg
Pkg.add('Pluto')
```

Una vez instalado, Pluto se inicia desde el intérprete de Julia ejecutando:

```
using Pluto
Pluto.run()
```

Al ejecutar este comando, Pluto lanza automáticamente un **servidor local** y abre una pestaña en el navegador web con la interfaz de notebooks.

Desde esta interfaz se pueden:

- crear notebooks nuevos
- abrir notebooks existentes
- gestionar el historial de notebooks recientes

---

## Abrir un notebook existente

Para abrir un notebook previamente creado (por ejemplo uno del pipeline de análisis EEG), basta con seleccionar el archivo `.jl` correspondiente desde la interfaz de Pluto.

También es posible abrirlo directamente desde Julia:

```
using Pluto
Pluto.run(notebook='ruta/al/notebook.jl')
```

Pluto cargará el notebook y ejecutará automáticamente las celdas necesarias respetando las dependencias entre ellas.

---

## Cerrar Pluto.jl

Para cerrar Pluto existen varias opciones:

- **Desde la interfaz web:** utilizando el botón *Shutdown* del servidor.
- **Desde el terminal o REPL de Julia:** presionando `Ctrl + C`.

Al detener el servidor, todas las sesiones activas de Pluto se cierran y el navegador deja de estar conectado al entorno de ejecución.

---

## Ventajas para el desarrollo del pipeline

En este proyecto, Pluto se utiliza no solo como herramienta de ejecución sino también como **documentación interactiva del pipeline**. Esto permite:

- describir cada fase del procesamiento EEG junto con el código correspondiente
- visualizar resultados intermedios y finales
- mantener una estructura reproducible del análisis

De este modo, el notebook actúa simultáneamente como **documentación, entorno de experimentación y registro del análisis computacional**.

---

#### Reproducibilidad automática del entorno

Pluto integra directamente el sistema de gestión de paquetes de Julia. Cada notebook puede almacenar internamente la información sobre las dependencias utilizadas en el análisis.

Esto permite que otro usuario pueda abrir el notebook y reproducir automáticamente el entorno de trabajo sin necesidad de instalar manualmente librerías específicas.

Este enfoque reduce significativamente los problemas habituales de reproducibilidad asociados con diferencias de versiones de paquetes o configuraciones del sistema.

---

#### Integración natural con Julia

Pluto está diseñado específicamente para el ecosistema Julia, lo que permite aprovechar de forma directa:

- bibliotecas de computación científica  
- procesamiento de señales  
- álgebra lineal de alto rendimiento  
- visualización avanzada  

Esto resulta especialmente útil en pipelines de análisis de EEG que incluyen:

- procesamiento de grandes matrices de datos  
- transformadas de Fourier  
- filtrado de señales  
- cálculo de métricas de conectividad entre múltiples canales  

El alto rendimiento de Julia permite realizar estas operaciones directamente en el notebook sin necesidad de optimizaciones externas.

---

#### Documentación científica integrada

Pluto facilita la creación de documentos científicos interactivos que combinan:

- explicación teórica  
- código reproducible  
- resultados visuales  

Esto permite utilizar el notebook no solo como herramienta de análisis, sino también como **documentación ejecutable del pipeline completo**.

En el contexto de este proyecto, Pluto se utiliza para:

- documentar cada etapa del pipeline EEG  
- visualizar resultados intermedios y finales  
- explicar los fundamentos metodológicos del análisis  
- generar figuras utilizadas en el informe

---

### Beneficios para la investigación reproducible

El uso de Pluto.jl aporta varias ventajas importantes para proyectos de investigación computacional:

- **reproducibilidad garantizada** del entorno y del flujo de ejecución  
- **transparencia del análisis**, al integrar código y explicación en el mismo documento  
- **facilidad para compartir notebooks completos** con otros investigadores  
- **exploración interactiva de datos y resultados**

En conjunto, Pluto permite que el pipeline de análisis EEG desarrollado en este proyecto sea **fácil de entender, reproducir y extender por otros investigadores**.

---

## Estructura y organización del proyecto

El código del proyecto está organizado como un único entorno de Julia denominado `EEG_JULIA`. Este entorno contiene todas las dependencias necesarias para ejecutar el pipeline completo de análisis.

La estructura general del proyecto separa claramente los distintos componentes del flujo de trabajo:

- **configuración**  
- **datos de entrada**  
- **resultados generados**  
- **módulos de código fuente**

Esta organización facilita la **reproducibilidad del análisis**, permitiendo ejecutar el pipeline completo a partir del mismo entorno de Julia y garantizando que todas las dependencias se encuentren en versiones controladas.

| Directorio | Contenido y propósito |
|---|---|
| `config/` | Archivos de configuración como `default_config.jl`: rutas, listas de sujetos y parámetros que controlan el pipeline sin necesidad de modificar el código. |
| `data/` | Datos de entrada e intermedios para cada etapa del procesamiento (por ejemplo `raw/`, `filtering/`, `ICA/`, `segmentation/`, `FFT/`, `CSD/`, `wPLI/`). Mantiene la organización del disco alineada con las etapas del pipeline. |
| `results/` | Resultados generados agrupados por tipo: `figures/` (CSD, FFT, ICA_cleaning, wPLI), `logs/` y `tables/`. Incluye todo aquello que puede regenerarse a partir del código y los datos. |
| `src/` | Módulos principales de Julia, generalmente un archivo `.jl` por cada etapa del pipeline (IO, filtrado, ICA, segmentación, FFT, baseline, rechazo de artefactos y módulos de conectividad como `CSD.jl` o `wPLI.jl`). |
| `script/` | Scripts ejecutables (por ejemplo drivers del pipeline o scripts de análisis) que integran los módulos definidos en `src/`. |

---

## Generación de figuras (Makie)

Todas las figuras de resultados (matrices de conectividad, espectros, comparaciones entre grupos) se generan directamente en Julia utilizando la librería **Makie** y se exportan como archivos **PDF vectoriales** para su inclusión directa en este informe.

Makie fue elegida en lugar de librerías como **Plots.jl** por varias razones técnicas.

Plots.jl proporciona una API unificada sobre múltiples backends gráficos (GR, PyPlot, entre otros) y es muy conveniente para visualización rápida o exploratoria. Sin embargo, cuando se trabaja con:

- matrices grandes (por ejemplo matrices canal × canal de conectividad),
- series temporales largas,
- o figuras complejas con múltiples paneles,

su rendimiento puede ser limitado y el comportamiento puede variar entre backends.

Makie, en cambio, es una **arquitectura gráfica moderna basada en un único stack**, lo que garantiza mayor consistencia y rendimiento. Sus principales ventajas son:

- **renderizado acelerado por GPU**
- **manejo eficiente de grandes matrices**
- **salida vectorial de alta calidad (PDF o SVG)**
- **control fino del diseño de figuras multi-panel**

Además, Makie utiliza un sistema de **scene graph y layouts jerárquicos**, que facilita la construcción de figuras complejas, como por ejemplo:

- una matriz de conectividad por banda de frecuencia  
- comparaciones entre grupos experimentales  
- paneles múltiples con escalas compartidas  

Para este proyecto se utiliza el backend **CairoMakie**, optimizado para la generación de figuras estáticas en formato vectorial.

---

## Ejemplo mínimo de generación de figura

El siguiente ejemplo muestra el flujo básico para crear una figura y exportarla a PDF.

```
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel = 'Channel',
    ylabel = 'Channel'
)

heatmap!(ax, rand(8,8))

save('connectivity-schematic.pdf', fig)
```

Este mismo patrón se utiliza en el pipeline para visualizar matrices reales de conectividad wPLI y figuras comparativas entre grupos, utilizando escalas de color y etiquetas apropiadas para cada análisis.
"

# ╔═╡ 0aa3b638-ab24-4fbb-8dfc-f5bd5e931f0c
md"
# Control de versiones con Git y GitHub

El **control de versiones** se gestiona mediante **Git**, lo que permite registrar los cambios realizados en el código del proyecto, mantener un historial completo del desarrollo y facilitar la colaboración o reutilización del pipeline en diferentes equipos.

El repositorio del proyecto está alojado en **GitHub**, lo que proporciona varias ventajas importantes:

- **Respaldo automático del código**
- **Historial completo de cambios**
- **Facilidad para compartir el proyecto**
- **Posibilidad de clonar el repositorio en otros ordenadores**
- **Gestión de colaboración entre múltiples desarrolladores**

Gracias a este sistema, cualquier versión anterior del código puede recuperarse fácilmente y es posible rastrear cuándo y por qué se realizaron cambios específicos en el pipeline de análisis.

| Escenario | Comandos | Explicación |
|---|---|---|
| Crear un nuevo repositorio en GitHub | `git init` | Inicializa un nuevo repositorio Git en la carpeta actual. |
|  | `git remote add origin <url>` | Conecta el repositorio local con el repositorio remoto en GitHub. |
|  | `git add .` | Añade todos los archivos al área de preparación (*staging*). |
|  | `git commit -m 'mensaje'` | Crea un snapshot de los cambios preparados con un mensaje descriptivo. |
|  | `git push -u origin main` | Sube la rama `main` a GitHub y la establece como rama remota por defecto. |
| Clonar un repositorio existente de GitHub | `git clone <url>` | Descarga el repositorio remoto en una nueva carpeta local. |
|  | `cd <folder>` | Cambia al directorio del proyecto clonado. |
|  | `julia -e 'using Pkg; Pkg.instantiate()'` | Instala el entorno de Julia definido por `Project.toml` y `Manifest.toml`. |
|  | *(posteriormente: `git push`, `git pull`)* | Permite subir commits locales o descargar cambios remotos. |

---

## Organización del repositorio

El repositorio contiene únicamente los elementos necesarios para reproducir el análisis:

- código fuente en Julia  
- archivos de configuración  
- scripts del pipeline  
- documentación y notebooks  
- pequeños archivos de metadatos  

Los **datos grandes o generados automáticamente** no se incluyen en el control de versiones para evitar repositorios excesivamente pesados y mejorar la eficiencia del trabajo colaborativo.

Por este motivo, directorios como:

```
data/
results/
```

se excluyen del repositorio mediante el archivo:

```
.gitignore
```

Este archivo indica a Git qué archivos o carpetas deben ignorarse durante el seguimiento de cambios.

| Qué incluir en el repositorio | Por qué |
|---|---|
| `src/` | Código fuente principal (módulos y funciones del pipeline). |
| `script/` | Scripts ejecutables utilizados para lanzar o coordinar el pipeline. |
| `config/` | Archivos de configuración que definen parámetros y rutas de forma portable. |
| `Project.toml` | Declara las dependencias del proyecto para que otros puedan reproducir el entorno. |
| `Manifest.toml` *(opcional)* | Fija las versiones exactas de las dependencias para garantizar reproducibilidad completa entre máquinas. |

### `.gitignore` recomendado

Un archivo `.gitignore` mínimo para un proyecto en **Julia** debería excluir artefactos generados automáticamente, configuraciones locales del editor y salidas de datos grandes.

| Patrón / ruta | Por qué ignorarlo |
|---|---|
| `*.jl.cov` | Archivos de cobertura generados durante las pruebas. |
| `Manifest.toml` | Opcional: se puede ignorar si se prefiere un entorno flexible sin versiones bloqueadas (inclúyelo si quieres reproducibilidad completa). |
| `data/` | Datos brutos o intermedios, que suelen ser grandes y no adecuados para control de versiones. |
| `results/` | Resultados generados automáticamente (figuras, exportaciones o cálculos cacheados). |
| `.DS_Store` | Archivos de metadatos generados por macOS. |
| `.vscode/`, `.cursor/` | Configuración local de editores o IDE (ignorar salvo que quieras compartir configuración de equipo). |

Ignorar estos archivos ayuda a mantener el repositorio pequeño y evita commits accidentales de archivos grandes o temporales.

---

## Manejo de archivos grandes

En los casos en que ciertos archivos de datos deban mantenerse versionados (por ejemplo, modelos preentrenados o datasets pequeños pero importantes), puede utilizarse **Git Large File Storage (Git LFS)**.

Git LFS permite gestionar archivos grandes almacenando en el repositorio únicamente referencias ligeras, mientras que el contenido real se mantiene en un almacenamiento externo optimizado para este tipo de datos.

---

## Beneficios para la reproducibilidad científica

El uso de Git y GitHub contribuye de forma significativa a la **reproducibilidad del proyecto**:

- permite reconstruir exactamente qué versión del código produjo determinados resultados  
- facilita la revisión y auditoría del pipeline de análisis  
- simplifica la distribución del proyecto a otros investigadores  

En conjunto, el control de versiones garantiza que el desarrollo del pipeline de análisis EEG sea **transparente, trazable y reproducible**.
"

# ╔═╡ 749a5ee5-b336-43e7-9cea-fbe323eb478b
md"""
# Uso responsable de herramientas IA

Durante la elaboración de este informe y el desarrollo del código en **Julia**, se utilizaron herramientas de **inteligencia artificial (IA)** como apoyo al proceso de trabajo. Estas herramientas se emplearon de forma **transparente, responsable y supervisada**, con el objetivo de mejorar la productividad, facilitar la exploración de ideas y apoyar la redacción técnica.

Es importante destacar que las herramientas de IA se utilizaron **como asistencia**, no como sustituto del razonamiento científico ni del trabajo personal. Todas las decisiones metodológicas, la implementación del pipeline y la interpretación de los resultados corresponden al autor del proyecto.

---

## Cursor

El editor **Cursor** (<https://cursor.com>) fue utilizado como entorno de desarrollo para el código en Julia. Cursor integra un asistente basado en inteligencia artificial que permite interactuar con el código del proyecto de forma contextual.

Este entorno facilitó varias tareas durante el desarrollo del pipeline:

- escritura inicial de fragmentos de código  
- refactorización y mejora de funciones existentes  
- ayuda en la depuración (*debugging*)  
- navegación y comprensión de la estructura del proyecto  
- generación o adaptación de pequeños fragmentos de documentación  

El asistente integrado permitió acelerar tareas repetitivas o mecánicas, manteniendo siempre la revisión manual del código generado.

---

## ChatGPT

**ChatGPT (OpenAI)** se utilizó como herramienta de apoyo conceptual y de redacción durante el desarrollo del proyecto.

Entre los usos principales se incluyen:

- clarificación de conceptos metodológicos (por ejemplo **BIDS**, **wPLI**, métodos **surrogate**, o aspectos de procesamiento de EEG)  
- generación de ideas preliminares para estructurar secciones del informe  
- redacción inicial de resúmenes o explicaciones técnicas  
- sugerencias de organización para la documentación del pipeline  

Las respuestas generadas por ChatGPT se utilizaron únicamente como **punto de partida para el trabajo personal**. Todo el contenido final fue revisado, editado y adaptado manualmente para garantizar su exactitud científica y coherencia con el proyecto.

---

## Consideraciones éticas y de reproducibilidad

El uso de herramientas de inteligencia artificial se realizó siguiendo principios de **transparencia y responsabilidad académica**. En particular:

- las herramientas de IA se emplearon como apoyo, no como sustituto del análisis científico  
- el diseño metodológico y la implementación final del código fueron realizados y verificados manualmente  
- las descripciones del pipeline se basan en el código y en la comprensión directa del proceso de análisis  

De este modo, el uso de IA contribuyó principalmente a **mejorar la claridad y eficiencia del proceso de desarrollo**, sin comprometer la integridad científica del trabajo.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
MakieThemes = "e296ed71-da82-5faf-88ab-0034a9761098"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyMNE = "6c5003b2-cbe8-491c-a0d1-70088e6a0fd6"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
TopoPlots = "2bdbdf9c-dbd8-403f-947b-1a4e0dd41a7a"
Unfold = "181c99d8-e21b-4ff3-b70b-c233eddec679"
UnfoldMakie = "69a5ce3b-64fb-4f22-ae69-36dd4416af2a"
UnfoldSim = "ed8ae6d2-84d3-44c6-ab46-0baf21700804"

[compat]
CSV = "~0.10.16"
CairoMakie = "~0.15.9"
DSP = "~0.8.4"
DataFrames = "~1.8.1"
DataFramesMeta = "~0.15.6"
FFTW = "~1.10.0"
MakieThemes = "~0.1.5"
Plots = "~1.41.6"
PlutoUI = "~0.7.79"
PyMNE = "~0.2.4"
StatsBase = "~0.34.10"
StatsPlots = "~0.15.8"
TopoPlots = "~0.3.0"
Unfold = "~0.8.9"
UnfoldMakie = "~0.5.22"
UnfoldSim = "~0.5.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "6bc127e4220ea8a9f9be8bff7c68427e3bc56f19"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "35ea197a51ce46fcd01c4a44befce0578a1aaeca"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AlgebraOfGraphics]]
deps = ["Accessors", "Colors", "DataAPI", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "Isoband", "KernelDensity", "Loess", "Makie", "NaturalSort", "PlotUtils", "PolygonOps", "PooledArrays", "PrecompileTools", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "f97855293b04707e1cc848e5ade97ffa473e8795"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.12.4"

    [deps.AlgebraOfGraphics.extensions]
    AlgebraOfGraphicsDynamicQuantitiesExt = "DynamicQuantities"
    AlgebraOfGraphicsUnitfulExt = "Unitful"

    [deps.AlgebraOfGraphics.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "libblastrampoline_jll"]
git-tree-sha1 = "7f54761502ff149a9d492e4acefe9805898e29b3"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.2+0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "78b3a7a536b4b0a747a0f296ea77091ca0a9f9a3"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.23.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "e0b47732a192dd59b9d079a06d04235e2f833963"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.12.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "ConcurrentUtilities", "DataAPI", "Dates", "EnumX", "Mmap", "PooledArrays", "SentinelArrays", "StringViews", "Tables", "TimeZones", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "4a69a3eadc1f7da78d950d1ef270c3a62c1f7e01"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.8.1"

[[deps.ArrowTypes]]
deps = ["Sockets", "UUIDs"]
git-tree-sha1 = "404265cd8128a2515a81d5eae16de90fdef05101"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "2.3.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.BSplineKit]]
deps = ["ArrayLayouts", "BandedMatrices", "FastGaussQuadrature", "ForwardDiff", "LinearAlgebra", "PrecompileTools", "Random", "Reexport", "SparseArrays", "Static", "StaticArrays", "StaticArraysCore", "StatsAPI"]
git-tree-sha1 = "02d491054afeb89b7f34331701e4474eb0b904f7"
uuid = "093aae92-e908-43d7-9660-e50ee39d5a0a"
version = "0.19.2"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "02fa77c70ba84361b9bc9ff28523bd9d78519265"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.11.0"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "091d591a060e43df1dd35faab3ca284925c48e46"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.3.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8d8e0b0f350b8e1c91420b5e64e5de774c2f0f4d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.16"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "fa072933899aae6dc61dde934febed8254e66c6a"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.15.9"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a21c5464519504e41e0cbc91f0188e8ca23d7440"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.CategoricalArrays]]
deps = ["Compat", "DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "73acb4ed51b1855e1b5ce5c610334363a98d13f1"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "1.0.2"
weakdeps = ["Arrow", "JSON", "RecipesBase", "SentinelArrays", "StatsBase", "StructTypes"]

    [deps.CategoricalArrays.extensions]
    CategoricalArraysArrowExt = "Arrow"
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStatsBaseExt = "StatsBase"
    CategoricalArraysStructTypesExt = "StructTypes"

[[deps.Chain]]
git-tree-sha1 = "765487f32aeece2cf28aa7038e29c31060cb5a69"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "1.0.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChunkSplitters]]
git-tree-sha1 = "1c52c8e2673edc030191177ff1aee42d25149acb"
uuid = "ae650224-84b6-46f8-82ea-d812ca08434e"
version = "3.2.0"

[[deps.CloughTocher2DInterpolation]]
deps = ["LinearAlgebra", "MiniQhull", "PrecompileTools", "StaticArrays"]
git-tree-sha1 = "7b250c72cb6ec97bd77fa025e350ba848f623ce9"
uuid = "b70b374f-000b-463f-88dc-37030f004bd0"
version = "0.1.1"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "b7231a755812695b8046e8471ddc34c8268cbad5"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "3.0.0"

[[deps.CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "d58afcd2833601636b48ee8cbeb2edcb086522c2"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "c761b00e7755700f9cdf5b02039939d1359330e1"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.1.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "3b4be73db165146d8a88e47924f464e55ab053cd"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.7"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.CondaPkg]]
deps = ["JSON", "Markdown", "MicroMamba", "Pidfile", "Pkg", "Preferences", "Scratch", "TOML", "pixi_jll"]
git-tree-sha1 = "0300af904a8c8d41ff715a60a6959136d22b8572"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.34"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "a692f5e257d332de1e554e4566a4e5a8a72de2b2"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.4"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "PrettyTables", "Reexport", "TableMetadataTools"]
git-tree-sha1 = "b0652fb7f3c094cf453bf22e699712a0bed9fc83"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.15.6"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.Delaunator]]
git-tree-sha1 = "93053ec77347697441fa764dca35ebe0b41be6c3"
uuid = "466f8f70-d5e3-4806-ac0b-a54b75a91218"
version = "0.1.3"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "a55766a9c8f66cf19ffcdbdb1444e249bb4ace33"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.4.6"

[[deps.Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "7da4b14cc4c3443a1afc64abee17f4fcb45ad837"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.4"

[[deps.Dierckx_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3251f44b3cac6fec4cec8db45d3ab0bfed51c4d8"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.2.0+0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "7ae99144ea44715402c6c882bfef2adbeadbc4ce"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.16"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Effects]]
deps = ["Combinatorics", "DataFrames", "Distributions", "ForwardDiff", "LinearAlgebra", "Statistics", "StatsAPI", "StatsBase", "StatsModels", "Tables"]
git-tree-sha1 = "7fa33d3a2b536008e45c106ba5e776156d5fac79"
uuid = "8f03c58b-bd97-4933-a826-f71b64d2cca2"
version = "1.7.0"
weakdeps = ["GLM", "MixedModels"]

    [deps.Effects.extensions]
    EffectsGLMExt = "GLM"
    EffectsMixedModelsExt = ["GLM", "MixedModels"]

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "75e5697f521c9ab89816d3abeea806dfc5afb967"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.12"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "0044e9f5e49a57e88205e8f30ab73928b05fe5b6"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.1.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "eef4c86803f47dcb61e9b8790ecaa96956fdd8ae"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.2"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "3bcb30438ee1655e3b9c42d97544de7addc9c589"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.3"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "7528a7956248c723d01a0a9b0447bf254bf4da52"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.5"

[[deps.GeoInterface]]
deps = ["DataAPI", "Extents", "GeoFormatTypes"]
git-tree-sha1 = "2b0312a0c06b4408773c6dc1829b472ea706f058"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.6.1"
weakdeps = ["GeometryBasics", "Makie", "RecipesBase"]

    [deps.GeoInterface.extensions]
    GeoInterfaceMakieExt = ["Makie", "GeometryBasics"]
    GeoInterfaceRecipesBaseExt = "RecipesBase"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"
weakdeps = ["GeoInterface"]

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "MPIPreferences", "Mmap", "Preferences", "Printf", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "e856eef26cf5bf2b0f95f8f4fc37553c72c8641c"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.17.2"

    [deps.HDF5.extensions]
    MPIExt = "MPI"

    [deps.HDF5.weakdeps]
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "e94f84da9af7ce9c6be049e9067e511e17ff89ec"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.6+0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.Highlights]]
deps = ["DocStringExtensions", "InteractiveUtils", "REPL"]
git-tree-sha1 = "9e13b8d8b1367d9692a90ea4711b4278e4755c32"
uuid = "eafb193a-b7ab-5a9e-9068-77385905fa72"
version = "0.5.3"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XML2_jll", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "157e2e5838984449e44af851a52fe374d56b9ada"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.13.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "52116260a234af5f69969c5286e6a5f8dc3feab8"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.12"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "dfde81fafbe5d6516fb864dc79362c5c6b973c82"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.10.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"
weakdeps = ["ArrowTypes", "Parsers"]

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "02b61501dbe6da3b927cc25dacd7ce32390ee970"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.2"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "ScopedValues", "TranscodingStreams"]
git-tree-sha1 = "d97791feefda45729613fafeccc4fbef3f539151"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.5.15"
weakdeps = ["UnPack"]

    [deps.JLD2.extensions]
    UnPackExt = "UnPack"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"
weakdeps = ["ArrowTypes"]

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"
weakdeps = ["ArrowTypes"]

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "9ea3422d03222c6de679934d1c08f0a99405aa03"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.5.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI", "StatsFuns"]
git-tree-sha1 = "b1ad83b367b915e2dc485dee3d62a6a6317d7ad4"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.5"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "191686b1ac1ea9c89fc52e996ad15d1d241d1e33"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.10.1+0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MLBase]]
deps = ["IterTools", "Random", "Reexport", "StatsBase"]
git-tree-sha1 = "ac79beff4257e6e80004d5aee25ffeee79d91263"
uuid = "f0e99cf1-93fa-52ec-9ecc-5026115318e0"
version = "0.9.2"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "9341048b9f723f2ae2a72a5269ac2f15f80534dc"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.3.2+0"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "8e98d5d80b87403c311fd51e8455d4546ba7a5f8"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.12"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "36c2d142e7d45fb98b5f83925213feb3292ca348"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.5.5+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "68af66ec16af8b152309310251ecb4fbfe39869f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.9"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MakieThemes]]
deps = ["Colors", "Makie", "Random"]
git-tree-sha1 = "7df2be39e8f968dce9a31e8f97ed87981d28d472"
uuid = "e296ed71-da82-5faf-88ab-0034a9761098"
version = "0.1.5"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.MetaArrays]]
deps = ["Requires"]
git-tree-sha1 = "6647f7d45a9153162d6561957405c12088caf537"
uuid = "36b8f3f0-b776-11e8-061f-1f20094e1fc8"
version = "0.2.10"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "535656ce55266bfed0575cd051acc4f36dc869a0"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.15"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bc95bf4149bf535c09602e3acdf950d9b4376227"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+3"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MixedModels]]
deps = ["Arrow", "BSplineKit", "Compat", "DataAPI", "Distributions", "GLM", "JSON3", "LinearAlgebra", "Markdown", "MixedModelsDatasets", "NLopt", "PooledArrays", "PrecompileTools", "ProgressMeter", "Random", "SparseArrays", "StaticArrays", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels", "StructTypes", "Tables", "TypedTables"]
git-tree-sha1 = "9292e6f1a91ff530ecb3387ccf68ade3a57ad07c"
uuid = "ff71e718-51f3-5ec2-a782-8ffcbfa3c316"
version = "4.38.1"

    [deps.MixedModels.extensions]
    MixedModelsFiniteDiffExt = ["FiniteDiff"]
    MixedModelsForwardDiffExt = ["ForwardDiff"]
    MixedModelsPRIMAExt = ["PRIMA"]

    [deps.MixedModels.weakdeps]
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    PRIMA = "0a7d04aa-8ac2-47b3-b7a7-9dbd6ad661ed"

[[deps.MixedModelsDatasets]]
deps = ["Arrow", "Artifacts", "LazyArtifacts"]
git-tree-sha1 = "ac0036e4f1829db000db46aad4cd5a207bba8465"
uuid = "7e9fb7ac-9f67-43bf-b2c8-96ba0796cbb6"
version = "0.1.2"

[[deps.MixedModelsSim]]
deps = ["LinearAlgebra", "MixedModels", "PooledArrays", "PrettyTables", "Random", "Statistics", "Tables"]
git-tree-sha1 = "9bbcfdc15ed571358e9dafe160dcfd44f0cac353"
uuid = "d5ae56c5-23ca-4a1f-b505-9fc4796fc1fe"
version = "0.2.13"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "2c140d60d7cb82badf06d8783800d0bcd1a7daa2"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.8.1"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.MyterialColors]]
git-tree-sha1 = "01d8466fb449436348999d7c6ad740f8f853a579"
uuid = "1c23619d-4212-4747-83aa-717207fae70f"
version = "0.3.0"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NLopt]]
deps = ["CEnum", "NLopt_jll"]
git-tree-sha1 = "624785b15005a0e0f4e462b27ee745dbe5941863"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "1.2.1"

    [deps.NLopt.extensions]
    NLoptMathOptInterfaceExt = ["MathOptInterface"]

    [deps.NLopt.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b0154a615d5b2b6cf7a2501123b793577d0b9950"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.10.0+0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NaturalNeighbours]]
deps = ["ChunkSplitters", "DelaunayTriangulation", "ElasticArrays", "LinearAlgebra", "Random"]
git-tree-sha1 = "c07721e921b9973c550e8560725bef98580e6555"
uuid = "f16ad982-4edb-46b1-8125-78e5a8b5a9e6"
version = "1.3.6"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "e2c3bba08dd6dedfe17a17889131b885b8c082f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.27"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f2b3b9e52a5eb6a3434c8cca67ad2dde011194f4"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.30+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "df9b7c88c2e7a2e77146223c526bf9e236d5f450"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML", "Zlib_jll"]
git-tree-sha1 = "2f3d05e419b6125ffe06e55784102e99325bdbe2"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "5.0.10+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "48968edaf014f67e58fe4c8a4ce72d392aed3294"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.3"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Peaks]]
deps = ["RecipesBase", "SIMD"]
git-tree-sha1 = "75d0ce1c30696d77bc60840222d7fc5d549ebf5f"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.5.3"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "211530a7dc76ab59087f4d4d1fc3f086fbe87594"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.2.3"

    [deps.PrettyTables.extensions]
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "f0803bc1171e455a04124affa9c21bba5ac4db32"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.6"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.PyMNE]]
deps = ["CondaPkg", "PythonCall", "Reexport"]
git-tree-sha1 = "2135adc42522d88293c797a7b15bfc7d0c856831"
uuid = "6c5003b2-cbe8-491c-a0d1-70088e6a0fd6"
version = "0.2.4"

    [deps.PyMNE.extensions]
    PyMNEOndaExt = ["Onda", "TimeSpans"]

    [deps.PyMNE.weakdeps]
    Onda = "e853f5be-6863-11e9-128d-476edb89bfb5"
    TimeSpans = "bb34ddd2-327f-4c4a-bfb0-c98fc494ece1"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "982f3f017f08d31202574ef6bdcf8b3466430bea"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.31"

    [deps.PythonCall.extensions]
    CategoricalArraysExt = "CategoricalArrays"
    PyCallExt = "PyCall"

    [deps.PythonCall.weakdeps]
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6f3ac0623e1173c006cc7377798ec3fb33fa504"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1004+0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "4d8c1b7c3329c1885b857abb50d08fa3f4d9e3c8"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.7"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.ScatteredInterpolation]]
deps = ["Combinatorics", "Distances", "LinearAlgebra", "NearestNeighbors"]
git-tree-sha1 = "0d642a08199bbeccd874b33fe3a1b699d345ca79"
uuid = "3f865c0f-6dca-5f4d-999b-29fe1e7e3c92"
version = "0.3.6"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "ac4b837d89a58c848e85e698e2a2514e9d59d8f6"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.6.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignalAnalysis]]
deps = ["DSP", "Distributions", "DocStringExtensions", "FFTW", "LinearAlgebra", "MetaArrays", "Optim", "PaddedViews", "Peaks", "PrecompileTools", "Random", "Reexport", "Requires", "SignalBase", "Statistics", "Test", "WAV"]
git-tree-sha1 = "a9f9b846649b264f2eca91049ab35e44c533f62a"
uuid = "df1fea92-c066-49dd-8b36-eace3378ea47"
version = "0.10.4"

[[deps.SignalBase]]
deps = ["Unitful"]
git-tree-sha1 = "14cb05cba5cc89d15e6098e7bb41dcef2606a10a"
uuid = "00c44e92-20f5-44bc-8f45-a1dcef76ba38"
version = "0.1.2"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "c06d695d51cfb2187e6848e98d6252df9101c588"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.3"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "SciMLPublic", "Static"]
git-tree-sha1 = "aa1ea41b3d45ac449d10477f65e2b40e3197a0d2"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.9.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "08786db4a1346d17d0a8d952d2e66fd00fa18192"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.9"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "88cf3587711d9ad0a55722d339a013c4c56c5bbc"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.8"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d05693d339e37d6ab134c5ab53c29fce5ee5d7d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.4"

[[deps.StringViews]]
git-tree-sha1 = "f2dcb92855b31ad92fe8f079d4f75ac57c93e4b8"
uuid = "354b36f9-a18e-4713-926e-db85100087ba"
version = "1.3.7"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "6dbb5b635c5437c68c28c2ac9e39b87138f37c0a"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.8"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "72df96b3a595b7aab1e101eb07d2a435963a97e2"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.5.0+2025b"

[[deps.TableMetadataTools]]
deps = ["DataAPI", "Dates", "TOML", "Tables", "Unitful"]
git-tree-sha1 = "c0405d3f8189bb9a9755e429c6ea2138fca7e31f"
uuid = "9ce81f87-eacc-4366-bf80-b621a3098ee2"
version = "0.1.0"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Term]]
deps = ["AbstractTrees", "CodeTracking", "Dates", "Highlights", "InteractiveUtils", "Logging", "Markdown", "MyterialColors", "OrderedCollections", "Parameters", "PrecompileTools", "ProgressLogging", "REPL", "Tables", "UUIDs", "Unicode", "UnicodeFun"]
git-tree-sha1 = "dfb5cd0af7f67f7487287fa7f83860fe0cb45899"
uuid = "22787eb5-b846-44ae-b979-8e399b8463ab"
version = "2.0.8"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "08c10bc34f4e7743f530793d0985bf3c254e193d"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.8"

[[deps.TiledIteration]]
deps = ["OffsetArrays", "StaticArrayInterface"]
git-tree-sha1 = "1176cc31e867217b06928e2f140c90bd1bc88283"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.5.0"

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "d422301b2a1e294e3e4214061e44f338cafe18a2"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.22.2"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.ToeplitzMatrices]]
deps = ["AbstractFFTs", "DSP", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "338d725bd62115be4ba7ffa891d85654e0bfb1a1"
uuid = "c751599d-da0a-543b-9d20-d0a503d91d24"
version = "0.8.5"
weakdeps = ["StatsBase"]

    [deps.ToeplitzMatrices.extensions]
    ToeplitzMatricesStatsBaseExt = "StatsBase"

[[deps.TopoPlots]]
deps = ["CloughTocher2DInterpolation", "Delaunator", "Dierckx", "GeometryBasics", "InteractiveUtils", "LinearAlgebra", "Makie", "NaturalNeighbours", "Parameters", "PrecompileTools", "ScatteredInterpolation", "Statistics"]
git-tree-sha1 = "2c9ac50c8a36aa4abc44fec11c8f3af000f976db"
uuid = "2bdbdf9c-dbd8-403f-947b-1a4e0dd41a7a"
version = "0.3.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TypedTables]]
deps = ["Adapt", "Dictionaries", "Indexing", "SplitApplyCombine", "Tables", "Unicode"]
git-tree-sha1 = "84fd7dadde577e01eb4323b7e7b9cb51c62c60d4"
uuid = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"
version = "1.4.6"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unfold]]
deps = ["DSP", "DataFrames", "Distributions", "DocStringExtensions", "Effects", "FileIO", "GLM", "ImageTransformations", "Interpolations", "IterativeSolvers", "JLD2", "LinearAlgebra", "Logging", "MLBase", "Missings", "OrderedCollections", "PooledArrays", "PrecompileTools", "ProgressMeter", "Random", "SimpleTraits", "SparseArrays", "StaticArrays", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels", "Suppressor", "Tables", "Term", "TimerOutputs", "TypedTables"]
git-tree-sha1 = "316cd21b661fa8958a18b4abfcb49104373fe688"
uuid = "181c99d8-e21b-4ff3-b70b-c233eddec679"
version = "0.8.9"

    [deps.Unfold.extensions]
    UnfoldBSplineKitExt = "BSplineKit"
    UnfoldCUDAExt = "CUDA"
    UnfoldKrylovExt = ["Krylov", "CUDA"]
    UnfoldRobustModelsExt = "RobustModels"

    [deps.Unfold.weakdeps]
    BSplineKit = "093aae92-e908-43d7-9660-e50ee39d5a0a"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Krylov = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
    RobustModels = "d6ea1423-9682-4bbd-952f-b1577cbf8c98"

[[deps.UnfoldMakie]]
deps = ["AlgebraOfGraphics", "BSplineKit", "CategoricalArrays", "ColorSchemes", "ColorTypes", "Colors", "CoordinateTransformations", "DataFrames", "DataStructures", "DocStringExtensions", "GeometryBasics", "GridLayoutBase", "ImageFiltering", "Interpolations", "LinearAlgebra", "Makie", "MakieThemes", "Random", "SparseArrays", "StaticArrays", "Statistics", "TopoPlots", "Unfold"]
git-tree-sha1 = "287269f4e0bce4cf81a4e114b75249ef210d717a"
uuid = "69a5ce3b-64fb-4f22-ae69-36dd4416af2a"
version = "0.5.22"
weakdeps = ["PyMNE", "UnfoldSim"]

    [deps.UnfoldMakie.extensions]
    UnfoldMakiePyMNEExt = "PyMNE"
    UnfoldMakieUnfoldSimExt = "UnfoldSim"

[[deps.UnfoldSim]]
deps = ["Artifacts", "Automa", "DSP", "DataFrames", "Distributions", "FileIO", "HDF5", "ImageFiltering", "LinearAlgebra", "MixedModels", "MixedModelsSim", "Parameters", "Random", "SignalAnalysis", "Statistics", "StatsModels", "ToeplitzMatrices"]
git-tree-sha1 = "8a4eeeec174d05b8488e5ee1f166309a06bbbe69"
uuid = "ed8ae6d2-84d3-44c6-ab46-0baf21700804"
version = "0.5.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "NaNMath", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.WAV]]
deps = ["Base64", "FileIO", "Libdl", "Logging"]
git-tree-sha1 = "7e7e1b4686995aaf4ecaaf52f6cd824fa6bd6aa5"
uuid = "8149f6b0-98f6-5db9-b78f-408fbbb8ef88"
version = "1.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "e9aeb174f95385de31e70bd15fa066a505ea82b9"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.7"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "80d3930c6347cfce7ccf96bd3bafdf079d9c0390"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.9+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "13b760f97c6e753b47df30cb438d4dc3b50df282"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.5+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "717df6f6892af4ee13279a73aa58474e58a88667"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "2.3.1+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.pixi_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "f349584316617063160a947a82638f7611a8ef0f"
uuid = "4d7b5844-a134-5dcd-ac86-c8f19cd51bed"
version = "0.41.3+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╠═1c2137c1-9117-405f-97dd-d2392fcc5911
# ╠═60620496-6fc4-45df-a2fc-eca64e76a52f
# ╟─f36b8ac2-d883-4d32-a7dc-2072e987165d
# ╟─a73289ad-5d2c-4452-b19e-8a0ef3e42836
# ╟─fb6880f3-5181-4299-af69-2675820284ca
# ╟─b8e9ae7c-942f-4f7c-9ed8-02dbbd0cb657
# ╟─3b1ca46d-580a-4a4c-83fc-d5a8249532a0
# ╟─7b611338-d3c2-4e05-bf4e-bb87b1b56e5f
# ╟─dadc3838-9f52-4397-a36d-e74626f873cb
# ╟─aa3253d3-8a64-4e0f-9e68-dd9b1c07db2d
# ╟─43e0934e-d7e9-421d-a678-ee7766bc90bc
# ╟─c62af5d5-c30d-44b7-a067-e01bbdc5e835
# ╟─d4f94266-9f5f-469c-abc9-90ed2995b95a
# ╠═b15557a3-4733-4bc7-85c7-313b08aa1039
# ╟─15966fd1-e39d-4c1c-ba38-f10a09e598dd
# ╟─8d4caa60-3407-4510-9db4-8132946572fb
# ╟─4001ae18-0212-4445-8950-75c3c21f40d3
# ╟─3ca6e623-ba24-413d-b383-4f82bbe1de87
# ╟─f6f00a09-de36-4617-b669-a633ed9c9664
# ╟─42756cc1-b01a-43bc-847a-9263233984fd
# ╟─44d68fb3-8b55-45ae-a1dc-834acde8556e
# ╟─937f5801-e2d0-411c-a7e9-2d8358e8ac95
# ╟─64f4f3c1-8b34-48fb-8a5f-3c414f09ad17
# ╟─80bfe9b0-af73-4cb8-be31-65455652ad59
# ╟─96c6f9db-a65b-4769-8ea7-dcd4b26a4378
# ╟─aa0047fa-a8e4-4e59-ae92-4e9ee1344811
# ╟─a1e97a97-3e0b-4352-8211-6ca4a1daaa7a
# ╟─1f9ebe5a-ad63-4be0-8e92-7aa1ee2b5bea
# ╟─eef886e7-89e9-4fa9-8e5a-819a5396aab2
# ╟─1ff72c58-70fb-4f55-aa0f-a883c79f72bc
# ╟─86335a85-f9da-4ebf-91b0-84e6212048d3
# ╟─009b77aa-a435-47d6-8c01-c547427f7155
# ╟─6343e583-9f06-42bf-974a-ee722805c7f5
# ╟─51ea38ab-dc1d-4cad-bacd-1922f7dddbd2
# ╟─fb1e4bb0-51db-41d6-8086-dc625b2e45d6
# ╟─fb669a4f-0320-4f2c-9589-d91484eac516
# ╟─7a6de37e-6370-4b16-8da9-f857fc79c01f
# ╟─424c9b90-8d93-407a-9fc7-840e2c1eec52
# ╟─94f78fd4-f8f8-4dfb-a035-4e2b9b227602
# ╟─042b4a34-d0ed-42f2-9595-f533807b7218
# ╟─8a99df75-275d-48ab-be48-1aee6194b9ad
# ╟─a31e013d-243e-430b-9d5b-ceaeb84eb378
# ╟─1f95dd62-aa74-4cf4-8096-9ad6faa9dcdb
# ╟─8ffd90e5-8aa9-4bca-a68a-c7bb88c3b349
# ╟─05e66c38-32c9-4e5f-81bd-f066e0e9fcd9
# ╟─9f7c42a7-dc71-4db5-9b51-07c38cea79ba
# ╟─800c5880-f5de-45bd-8458-3666931ec0ff
# ╟─ac015a68-1c14-46e8-9ef2-b00962add288
# ╟─188a16df-19e8-4455-9bf1-15bedb31efea
# ╟─d38f84c5-f99d-4d21-9113-7b7697b9ac59
# ╟─2eb0aada-37e4-4b30-87d1-76822d782adc
# ╟─611b3e3b-9931-4d8b-b9d3-3ad5d048f60f
# ╟─830fd9d7-7916-4f56-bb3c-75102e4a26fd
# ╟─fd19d3b2-e552-457a-aa8c-dd170cb05d7c
# ╟─016a1c47-ed53-4dee-9322-274e599d8e88
# ╠═668144ad-6b12-4c89-8280-0a1ba5ed3609
# ╟─b1337098-3e65-4b6b-8c88-f5becac7fff3
# ╟─8db31d27-4aa6-40bc-8c83-fb70727981fa
# ╟─359e26ad-2713-4aff-873f-bd2e765be8d2
# ╟─b96bf1b3-2e15-49c6-b7d0-877f41569e8a
# ╟─a55e367e-be5a-43da-8bce-279ebaca8ca6
# ╟─b59da6b5-0b1a-4448-86fa-f5fd8196a5b6
# ╟─91d53d89-df18-42d5-8943-af8a7d6cb74b
# ╟─bc98c842-7ce3-4a31-bfa0-bec06f3997c3
# ╟─94d21c47-149b-4168-b090-410c45e16f4e
# ╟─0bbe4118-3e47-4854-98e8-73b79aecf55f
# ╟─c780bd60-ed60-4fc0-877a-36a85c328a1f
# ╟─a753bee5-a981-4a42-8dda-2b1e10d0f003
# ╟─47511ab1-2099-4385-93b0-32f1772713b8
# ╟─7b988af1-3aab-464e-b3e6-dc9a07128877
# ╟─bd16c06c-16bf-4368-8389-5bb413325649
# ╟─add5d704-abbd-4f5c-a626-05600a057bd4
# ╟─06f567b9-feda-4381-a032-ffe2176926e6
# ╟─64e75122-d2b9-425d-ad99-0d20c891d910
# ╟─2a43830c-9ffe-4dae-9cd2-e23decfdf230
# ╟─d6a7a4c0-5680-4fb4-bc91-4d57d5d3adc5
# ╟─4f6f1059-01e1-4213-af37-5fd304474212
# ╟─712f89de-087f-4e88-83de-40d58c886896
# ╟─c4a28722-e0b9-4726-ba15-2aa8bb28ebea
# ╟─07e4fff6-77b0-48c8-93da-e354deb017aa
# ╟─f608f7b6-8fad-4289-9909-8d57ce90a070
# ╟─2bcd6082-ca7b-4353-9d14-f20bc4aa20b7
# ╟─078ecfa6-0ca4-42d3-8ffb-867e45c65c1d
# ╟─1e3c7d2a-3f68-4107-81f9-8ef9f2245802
# ╟─af42ccde-b32d-4c24-aacc-7b96e158012d
# ╟─b73dff40-9bf4-439e-b413-011cefa8e4cf
# ╟─996bfc87-a4f3-494d-a809-d7df30857eb2
# ╟─5dcac750-2be2-4302-92ae-b01464a4ac2c
# ╟─a58fb4d7-12d7-4f82-a00a-926f6772ea5b
# ╟─0b2fbbd6-faed-48f8-b35c-f3fe7de7526a
# ╟─eeac27b1-7fff-46a3-acfa-3a540cc989ca
# ╟─47e5c550-a820-44a4-a1bc-ac05851665f4
# ╟─bc054369-e61b-4a71-92a1-a5d08ff69224
# ╟─3cacaf2c-57df-4eba-a29e-f5b82feeb79c
# ╟─230c3dbb-af27-4147-a972-ab1c45578fdb
# ╟─bfb33303-de8d-42ef-9d7e-e6b5e813bee8
# ╟─89849536-05c6-41a9-bb94-7eabe4057549
# ╟─4c757bb9-70a7-4fc2-b9c3-d17a019ba75e
# ╟─5839954c-87c2-4333-81a3-9a5f64021000
# ╟─9ea266c8-0ffb-49c9-94f8-b77d206b087f
# ╟─19a03e63-b353-4206-9deb-823d52bcd2f3
# ╟─2ba2934f-e388-462e-a632-8bc906cc6a4c
# ╟─4450bc71-a1f1-4f8b-b390-99a8221b9b7d
# ╟─4f237b4c-6598-4fcc-8dd3-2b9b07e141bd
# ╟─3713b98b-ae0a-440f-b3a3-3bfebe68a937
# ╟─97c5d3c1-41f1-42cb-ad01-192e7507c528
# ╟─95d47537-1faf-47d0-ac9a-4be72a90c832
# ╟─ece3c6d2-723c-4489-83eb-ff36cf9bc10e
# ╟─1c9a3df0-545c-4695-aef4-face8703449d
# ╟─20569852-d89f-49f0-8d73-30ad7d2261d9
# ╟─3cac8ebb-b5ab-4710-8a92-7627813b579e
# ╟─e79ba77a-2e08-47db-a1f9-91efe34c2510
# ╟─f8bd3fd9-93c3-4b16-8b68-5c780a013e9b
# ╟─c9854ffc-344b-4f87-bb93-2e71e71aa504
# ╟─e28cf442-e820-47ef-947c-85398eadadd7
# ╟─2c8d9444-4e8f-4234-b31d-c7e3425f8d13
# ╟─1c2bd8c6-3309-4a92-8a10-5437f5036c7a
# ╟─04e8aa59-fecb-401f-8adf-95374473dbb7
# ╟─47af6c89-b4c4-43e6-a2cf-1771cb45089e
# ╟─5b99601e-2fb5-4a21-9ab9-ead4df41c9e6
# ╟─35f3cbd4-8a04-4f72-9c43-d73a23c6de8c
# ╟─efc6b5dc-f489-48dc-bdbb-a79cdd744e86
# ╟─8713cdf3-d880-48f6-8575-a39bb16acb01
# ╟─fd21bfc3-e448-4b2c-a723-257a8d1aa6ca
# ╟─06031333-5d1f-4cef-b1c4-c7f8f4093e66
# ╟─c586474d-9b33-4f1d-96fb-7f948724f4fe
# ╟─0aa3b638-ab24-4fbb-8dfc-f5bd5e931f0c
# ╟─749a5ee5-b336-43e7-9cea-fbe323eb478b
# ╠═4dc6ee64-1d24-11f1-b5d2-1b41a361acaa
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
