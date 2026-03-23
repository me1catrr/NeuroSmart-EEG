# Publicacion GitHub Pages (`tools/pages`)

Guia unica del flujo de publicacion de Pages para este repositorio.

## Convencion del repositorio

- Rama principal unica: `main`
- Origen de publicacion en GitHub Pages: `main:/docs`
- Tooling de publicacion: `tools/pages/`
- Sin logica de publicacion en `src/`

## Scripts

- `tools/pages/build.jl`
  - Construye/actualiza `docs/` para Pages.
  - Copia HTML exportado si existe.
  - Genera placeholders donde no haya exportacion.
  - Regenera `docs/pages_manifest.json`.

- `tools/pages/publish.jl`
  - Ejecuta `build.jl`.
  - Muestra `git status --short`.
  - Si no hay cambios: termina correctamente.
  - Si hay cambios: `git add .` -> `git commit` -> `git push`.
  - Si `git push` falla por falta de upstream: reintenta con `git push --set-upstream origin <rama_actual>`.
  - Informa la rama actual y advierte (sin bloquear) si no es `main`.

## Comandos de uso

```bash
julia tools/pages/build.jl
```

```bash
julia tools/pages/publish.jl
```

```bash
julia tools/pages/publish.jl "mensaje personalizado"
```

## Flujo recomendado

1. Trabajar notebooks fuente en `Pluto/...`.
2. Exportar HTML de Pluto con tu flujo habitual.
3. Ejecutar `julia tools/pages/build.jl` para actualizar `docs/`.
4. Revisar localmente `docs/index.html` y `docs/Pluto/*/index.html`.
5. Ejecutar `julia tools/pages/publish.jl` para commit + push.

## Verificacion de actualizacion en GitHub Pages

Despues del push:

1. Esperar 1-3 minutos.
2. Revisar en GitHub el ultimo commit en `main`.
3. Abrir el sitio de Pages y confirmar que carga la version nueva de `docs/index.html`.
4. Verificar al menos un modulo en `docs/Pluto/.../index.html`.
