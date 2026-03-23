# README de GitHub Pages (NeuroSmart-EEG)

Guia principal para publicar notebooks de Pluto.jl como paginas HTML en GitHub Pages.

## 1) Principio de separacion

- **Fuente cientifica**: notebooks en `Pluto/...` (`.jl`).
- **Publicacion web**: HTML en `docs/Pluto/...`.
- La logica de publicacion vive en tooling (`tools/pages/`), no en `src/`.

## 2) Como publica GitHub Pages aqui

- GitHub Pages sirve contenido desde `docs/`.
- Portada principal: `docs/index.html`.
- Modulos publicables:
  - `docs/Pluto/BIDS/index.html`
  - `docs/Pluto/Connectivity/index.html`
  - `docs/Pluto/Preprocessing/index.html`
  - `docs/Pluto/Processing/index.html`
- Inventario de publicacion: `docs/pages_manifest.json`.

## 3) Flujo minimo recomendado

1. Editar notebook fuente en `Pluto/<Modulo>/<Modulo>.jl`.
2. Exportar notebook a HTML con tu flujo habitual de Pluto.
3. Ejecutar:
   - `julia tools/pages/build.jl <ruta_opcional_exportaciones>`
4. Confirmar que el modulo apunta a `docs/Pluto/<Modulo>/index.html`.
5. Commit + push.

> El script no inventa exportacion de Pluto: solo copia HTML existente o genera placeholders.

## 4) Automatizacion unificada

Script unico:

- `tools/pages/build.jl`

Responsabilidades:

- asegura carpetas `docs/Pluto/<Modulo>/`,
- copia HTML exportado si existe en una carpeta de entrada opcional,
- genera placeholder HTML si no hay exportado,
- regenera `docs/pages_manifest.json`,
- imprime resumen final (publicados y placeholders creados en la ejecucion).

## 5) Estructura recomendada

```text
EEG_JULIA/
├── Pluto/
│   ├── BIDS/BIDS.jl
│   ├── Connectivity/Connectivity.jl
│   ├── Preprocessing/Preprocessing.jl
│   └── Processing/Processing.jl
├── docs/
│   ├── .nojekyll
│   ├── index.html
│   ├── README_pages.md
│   ├── pages_manifest.json
│   └── Pluto/
│       ├── BIDS/index.html
│       ├── Connectivity/index.html
│       ├── Preprocessing/index.html
│       └── Processing/index.html
└── tools/
    └── pages/
        └── build.jl
```

## 6) Como anadir un modulo nuevo

1. Crear notebook fuente en `Pluto/NewModule/NewModule.jl`.
2. Exportar HTML y preparar destino en `docs/Pluto/NewModule/index.html`.
3. Actualizar la lista `MODULES` en `tools/pages/build.jl`.
4. Agregar tarjeta/enlace en `docs/index.html`.
5. Ejecutar `julia tools/pages/build.jl` para regenerar `pages_manifest.json`.
