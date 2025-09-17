# Carpeta `data/` — ignorada por Git (no versionada)

Esta carpeta está **excluida en `.gitignore`**, por lo que Git **no registra** su contenido ni lo sube a GitHub.

**Implicaciones**
- Los archivos que pongas aquí **no se comparten** con otras personas mediante el repositorio.
- **No quedan en el historial** de Git (no hay *diffs* ni recuperación desde Git si los borras).
- Haz copias/backup por otros medios (NAS, disco externo, DVC, etc.).

**Excepciones permitidas**
- Solo se incluyen en el repositorio estos archivos “marcadores” para documentar la estructura:
  - `data/**/README.md`
  - `data/**/.gitkeep`
  > Estos archivos **no contienen datos**.

**Comprobación rápida**
```bash
git status --ignored -s   # debería mostrar '!! data/' indicando que está ignorada

## Estructura esperada

data/
├── raw/                                    
    ├── all_data/                           # datos originales BrainVision (patients y controls) .eeg .vhdr .vmrk .ehst2 .hfinf2 
│   ├── patients/                           # Pacientes (.eeg .vhdr .vmrk .ehst2 .hfinf2)
│   ├── controls/                           # Controles (mismo formato que patients)
│   └── info/                               # Metadatos/diccionarios demográficos, etc.
│       ├── organize/                       # (opcional) utilidades locales para ordenar datos
│       │   ├── organize_from_csv.jl
│       │   └── organize_log.txt
│       ├── demograficos_t1_v2.csv
│       └── demograficos_controles_t1.csv
├── preprocessed/                           # Señal limpia, segmentada, etc. (.bson .txt)
│   ├── patients/
│   └── controls/
├── features/                               # Características ya calculadas (bandpower, etc.)
└── connectivity/                           # Métricas de conectividad

> Si tienes scripts genéricos para organizar datos, **mejor muévalos a** `scripts/` o `src/` para que queden versionados (esta carpeta está ignorada).

## Privacidad y cumplimiento

- Mantén esta carpeta **fuera del control de versiones**.
- Si compartes datasets internamente, usa almacenamiento seguro (NAS/Drive/SSH).  
  Considera DVC si quieres versionar *metadatos* de datasets sin subir archivos.

---

**Este `README.md` existe solo para documentar la estructura; no contiene datos.**