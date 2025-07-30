# NeuroSmart-EEG

Plataforma para la **monitorización cerebral basada en EEG** orientada al análisis de **enfermedades neurodegenerativas** y aplicaciones clínicas avanzadas. Este proyecto sigue una estructura modular inspirada en buenas prácticas y en la literatura (*EEG Signal Processing and Machine Learning*).

---

## 🎯 Objetivo
Desarrollar un framework reproducible para:
- Procesar y analizar señales EEG en contextos clínicos.
- Implementar métodos clásicos y avanzados (FFT, Wavelets, ICA, EMD, CEEMD).
- Incorporar métricas de **conectividad funcional** (Coherencia, PLV, ImCoh).
- Explorar biomarcadores para patologías como **Esclerosis Múltiple**, Alzheimer y Parkinson.
- Facilitar la transferencia metodológica desde **Structural Health Monitoring (SHM)**.

---

## 📂 Estructura del proyecto

NeuroSmart-EEG/
│
├── data/           # Datos EEG (NO se suben, confidencialidad)
├── notebooks/      # Análisis exploratorio (Jupyter/Pluto)
├── src/            # Código modular en Julia
├── scripts/        # Pipelines ejecutables
├── config/         # Parámetros y configuración
├── docs/           # Documentación científica y técnica
└── tests/          # Pruebas unitarias

### 🔒 **Carpeta data/**

**IMPORTANTE:** Esta carpeta **NO debe contener datos sensibles** en el repositorio público.

#### Estructura esperada:

- `raw/` → Datos EEG originales (.eeg, .vhdr, .vmrk) [EXCLUIDOS del repo]
- `preprocessed/` → Señales tras filtrado y segmentación
- `connectivity/` → Resultados de métricas funcionales
- `features/` → Características espectrales, EMD, etc.

#### Nota sobre confidencialidad:

Los datos clínicos asociados a este proyecto son **confidenciales** y **no deben ser subidos a GitHub**.  
Para pruebas, se pueden usar datasets públicos (ej. *PhysioNet EEG*, *BCI Competition*).

---

## 🔄 Flujo de trabajo del análisis EEG
Raw EEG (.eeg, .vhdr, .vmrk)
↓
Preprocesamiento (Filtros, ICA)
↓
Segmentación en epochs
↓
Análisis espectral (FFT, Wavelets)
↓
Conectividad funcional (Coherencia, PLV, ImCoh)
↓
Extracción de features
↓
Modelos estadísticos / ML

---

## ⚡ Instalación rápida

Clonar el repositorio y activar el entorno:
```bash
git clone git@github.com:me1catrr/NeuroSmart-EEG.git
cd NeuroSmart-EEG
julia --project=.
]
instantiate
```

## 🚀 Pipeline actual
	•	Preprocesamiento: conversión a µV, re-referenciación, submuestreo, filtrado IIR.
	•	ICA con revisión semiautomática para artefactos.
	•	Segmentación en epochs de 1s.
	•	Análisis espectral (FFT, Wavelets Morlet).
	•	Cálculo de potencia absoluta/relativa por bandas clásicas.

⸻

## 🔬 Análisis avanzados
	•	Conectividad funcional: Coherencia, PLV, ImCoh.
	•	Descomposición adaptativa: EMD, CEEMD.
	•	Comparación grupos: EM-FR vs Control.

⸻

## ⚠️ Datos y confidencialidad

Los datos clínicos NO se suben a este repositorio.
Solo se incluyen scripts y notebooks para reproducibilidad.
Para pruebas, se pueden usar datasets públicos (BCI Competition, PhysioNet EEG).

⸻

## 🛠 Tecnologías
	•	Lenguaje: Julia 1.x
	•	Paquetes: DSP.jl, Wavelets.jl, Flux.jl, Plots.jl
	•	Control de versiones: Git + GitHub (SSH)

⸻

## 👥 Equipo
	•	Investigador: Rafael Castro-Triguero
	•	Colaboradores: Grupo de investigación neurociencia computacional

⸻

## 📜 Licencia

Este proyecto se distribuye bajo licencia MIT modificada (ver LICENSE).