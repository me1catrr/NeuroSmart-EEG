# NeuroSmart-EEG

Plataforma para la **monitorización cerebral basada en EEG** orientada al análisis de **enfermedades neurodegenerativas** y aplicaciones clínicas avanzadas. Este proyecto sigue una estructura modular inspirada en repositorios similares y en alguna literatura existente al respecto (*EEG Signal Processing and Machine Learning* de *Saeid Sanei*).

---

## 🎯 Objetivo
Desarrollar una estructura de rutinas reproducible para:
- Procesar y analizar señales EEG en contextos clínicos.
- Implementar métodos clásicos y avanzados (FFT, Wavelets, ICA, EMD, CEEMD, etc…).
- Incorporar métricas de **conectividad funcional** (Coherencia, PLV, ImCoh, etc…).
- Explorar biomarcadores para patologías como **Esclerosis Múltiple**, **Alzheimer** y **Parkinson**.
- Facilitar la transferencia metodológica desde disciplinas de la ingeniería como el **Structural Health Monitoring (SHM)**.

---

## 📂 Estructura del proyecto

```
NeuroSmart-EEG/
│
├── data/        # Datos EEG (NO se suben por confidencialidad)
│   ├── raw/             # Señales EEG crudas (.eeg, .vhdr, .vmrk)
│   ├── preprocessed/    # Señales tras filtrado y segmentación
│   ├── features/        # Características extraídas (band power, conectividad)
│   └── connectivity/    # Métricas funcionales
│
├── notebooks/   # Análisis exploratorio (Jupyter/Pluto)
├── src/         # Código modular en Julia
│   ├── preprocessing/   # Preprocesamiento, segmentación, artefactos
│   ├── spectral/        # Análisis espectral (FFT, Wavelets)
│   ├── connectivity/    # Coherencia, PLV, ImCoh
│   ├── features/        # Band power, métricas resumen
│   ├── decomposition/   # ICA, EMD
│   ├── stats/           # Estadística descriptiva, inferencial
│   └── utils/           # Utilidades generales
│
├── scripts/     # Pipelines ejecutables (ej: run_preprocessing.jl)
├── config/      # Parámetros de análisis
├── docs/        # Documentación técnica y científica
└── tests/       # Pruebas unitarias (Julia Pkg.test)
```

### 🔒 **Carpeta data/**
**IMPORTANTE:** Este proyecto **NO contiene datos sensibles** en el repositorio público. La carpeta data/ no está disponible por confidencialidad.

#### Nota sobre confidencialidad:
Los datos clínicos asociados a este proyecto son **confidenciales** y **no deben ser subidos a GitHub**.
Para pruebas, se pueden usar datasets públicos (ej. *PhysioNet EEG*, *BCI Competition*).

---

## 🔄 Flujo de trabajo del análisis EEG
1. **Raw EEG** (.eeg, .vhdr, .vmrk)
2. **Preprocesamiento** (Filtros notch, band-pass, ICA, re-referenciación)
3. **Segmentación** en epochs (1 s)
4. **Análisis espectral** (FFT, Wavelets)
5. **Conectividad funcional** (Coherencia, PLV, ImCoh)
6. **Extracción de features**
7. **Modelos estadísticos / Machine Learning**

---

## ⚡ Instalación rápida

Clonar el repositorio y activar el entorno Julia:
```bash
git clone git@github.com:me1catrr/NeuroSmart-EEG.git
cd NeuroSmart-EEG
julia --project=.
]
instantiate
```

---

## 🚀 Pipeline actual
- ✅ Preprocesamiento: conversión a µV, filtrado IIR, notch 50 Hz, submuestreo (512 Hz).
- ✅ ICA: corrección semiautomática para artefactos.
- ✅ Segmentación en epochs de 1s.
- ✅ Análisis espectral: FFT y bandas clásicas (delta, theta, alfa, beta, gamma).
- ✅ Potencia absoluta y relativa por bandas.

Próximos pasos:
- Análisis avanzado: Wavelets, EMD, CEEMD.
- Conectividad: Coherencia, PLV.
- Modelos predictivos: SVM, Random Forest.

---

## 📊 Justificación técnica
- **Filtros:** Eliminan artefactos de baja y alta frecuencia y la interferencia eléctrica (50 Hz).
- **ICA:** Identifica y elimina fuentes independientes (parpadeos, artefactos musculares).
- **Segmentación:** Necesaria para análisis en ventanas temporales estandarizadas.
- **FFT:** Base del análisis espectral en neurociencia (potencia por banda).
- **Conectividad:** Relaciona regiones cerebrales funcionalmente.
- **Band power:** Indicador de estados cognitivos y patológicos.

---

## 🛠 Tecnologías
- Lenguaje: **Julia 1.x**
- Paquetes: DSP.jl, Wavelets.jl, Flux.jl, Plots.jl
- CI/CD: GitHub Actions
- Control de versiones: Git + GitHub (SSH)

---

## 👥 Equipo
- Investigadores: **Alejandro Galvao**, **Rafael Castro-Triguero**
- Grupo de investigación neurociencia computacional

---

## 📜 Licencia
Este proyecto se distribuye bajo licencia **MIT modificada** (ver LICENSE).
