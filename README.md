# NeuroSmart-EEG

Repositorio en **Julia** para organizar y analizar señales EEG en el contexto de enfermedades neurodegenerativas (como **Esclerosis Múltiple**) siguiendo el estándar **BIDS (Brain Imaging Data Structure)**.  

👉 **Importante:** este repositorio **no contiene datos EEG**. Solo incluye **código, configuraciones y plantillas**.

## 🔬 Información Técnica de Adquisición

### Sistema de Grabación EEG
- **Amplificador**: actiCHamp Base Unit (5001) + módulo 32 CH
- **Software**: BrainVision Recorder Professional v. 1.21.0303
- **Formato**: Brain Vision Data Exchange Header File v1.0

### Parámetros de Grabación
- **Canales**: 31 electrodos (sistema 10-20 internacional)
- **Frecuencia de muestreo**: 500 Hz
- **Resolución**: 0.0488281 µV por unidad digital
- **Duración**: 3-5 minutos por condición (ojos abiertos/cerrados)
- **Filtros hardware**: DC-140 Hz
- **Filtros software**: 0.63-70 Hz + notch 50 Hz

### Dataset
- **Sujetos**: 77 participantes
- **Sesiones**: 106 totales (44 completas, 62 parciales)
- **Condiciones**: ojos abiertos, ojos cerrados, tarea desconocida

📋 Para información técnica detallada, ver `bids/dataset_description.json`  

---

## 📂 Estructura principal

```
NeuroSmart-EEG/
├── src/          # código modular en Julia (funciones y módulos)
├── scripts/      # scripts ejecutables (pipelines, conversión a BIDS, análisis)
├── config/       # archivos de configuración y metadatos (sin datos)
├── sourcedata/   # datos crudos locales (no se versionan)
├── bids/         # dataset en formato BIDS (generado localmente, no se versiona)
└── derivatives/  # resultados de análisis (locales, no se versionan)
```

- **src/** → lógica principal, pensada como librería o módulo.  
- **scripts/** → programas listos para ejecutar, que llaman a funciones de `src/`.  
- **config/** → plantillas y metadatos necesarios (ej. `participants.tsv`, `demographics.csv`).  
- **sourcedata/** → almacenamiento local de los archivos originales (`.eeg`, `.vhdr`, `.vmrk`).  
- **bids/** → salida con la organización estándar BIDS.  
- **derivatives/** → resultados de preprocesado, análisis espectral, conectividad, etc.  

---

## 🚀 Uso básico

1. Clonar el repositorio y activar el entorno de Julia:

```bash
git clone git@github.com:me1catrr/NeuroSmart-EEG.git
cd NeuroSmart-EEG
julia --project=.
] instantiate
```

2. Colocar tus datos EEG en `sourcedata/` (no se suben a GitHub).  
3. Usar los scripts de `scripts/` para convertir y analizar en formato **BIDS**.  

---

## 🛠️ Scripts Disponibles

### Conversión y Organización
- `scripts/build_participants.jl` - Genera `participants.tsv` y `participants.json`
- `scripts/build_eeg_bids.jl` - Convierte archivos BrainVision a formato BIDS
- `scripts/validate_bids.jl` - Valida la estructura BIDS del dataset

### Análisis y Visualización
- `scripts/plot_raw_traces.jl` - Genera visualizaciones de trazas EEG sin procesar
  - Crea plots de control de calidad
  - Genera layout de electrodos
  - Incluye información demográfica en títulos

## 📊 Control de Calidad

El script `plot_raw_traces.jl` genera automáticamente:
- **Plots de trazas**: 10 segundos iniciales de cada registro
- **Layout de electrodos**: Visualización del montaje 10-20
- **Reporte de calidad**: Resumen del procesamiento
- **Metadatos**: Información demográfica (sexo, edad, grupo)

## 🔧 Dependencias

El proyecto utiliza las siguientes librerías de Julia:
- `GLMakie` - Visualización
- `DataFrames` - Manipulación de datos
- `CSV` - Lectura de archivos tabulares
- `JSON3` - Procesamiento JSON
- `IniFile` - Lectura de archivos .vhdr
- `EDF` - Soporte para archivos EDF (opcional)

## 📜 Licencia

Este proyecto se distribuye bajo licencia **MIT** (ver archivo `LICENSE`).  