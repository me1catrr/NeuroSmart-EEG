# NeuroSmart-EEG

Repositorio en **Julia** para organizar y analizar señales EEG en el contexto de enfermedades neurodegenerativas (como **Esclerosis Múltiple**) siguiendo el estándar **BIDS (Brain Imaging Data Structure)**.  

👉 **Importante:** este repositorio **no contiene datos EEG**. Solo incluye **código, configuraciones y plantillas**.  

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

## 📜 Licencia

Este proyecto se distribuye bajo licencia **MIT** (ver archivo `LICENSE`).  