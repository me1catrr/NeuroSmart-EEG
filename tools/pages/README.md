# Publicacion GitHub Pages (`tools/pages`)

Guia principal y unica del flujo de publicacion web para este repositorio Julia + Pluto.

## Que se publica

- GitHub Pages publica desde `main:/docs`.
- Portada: `docs/index.html`.
- Modulos HTML: `docs/Pluto/<Modulo>/index.html`.

## Que NO se publica

- Carpetas de resultados internos o no publicables.
- `Javier_results` no forma parte de la web publicada.
- Notebooks fuente `.jl` no se publican directamente como pagina web.

## Ubicacion de cada tipo de archivo

- Fuente cientifica: `Pluto/<Modulo>/<Modulo>.jl`
- Staging de HTML exportado manualmente: `exports/Pluto/<Modulo>/index.html`
- Destino publicado en Pages: `docs/Pluto/<Modulo>/index.html`

## Scripts

- `tools/pages/build.jl`
  - Busca exportaciones reales en `exports/Pluto/...` (o en una ruta custom opcional).
  - Si existe HTML real del modulo, lo copia a `docs/Pluto/...`.
  - Si no existe, genera un placeholder claro con enlace de vuelta a portada.
  - Regenera `docs/pages_manifest.json`.

- `tools/pages/publish.jl`
  - Ejecuta `build.jl`.
  - Muestra `git status --short`.
  - Si no hay cambios: termina correctamente.
  - Si hay cambios: `git add .` -> `git commit` -> `git push`.
  - Fallback unico para upstream ausente: `git push --set-upstream origin <rama_actual>`.

## Comandos

```bash
julia tools/pages/build.jl
```

```bash
julia tools/pages/build.jl "/ruta/custom/exports/Pluto"
```

```bash
julia tools/pages/publish.jl
```

```bash
julia tools/pages/publish.jl "mensaje"
```

## Flujo recomendado

1. Editar notebook fuente en `Pluto/<Modulo>/<Modulo>.jl`.
2. Exportar HTML desde Pluto.jl manualmente.
3. Guardar el HTML exportado en `exports/Pluto/<Modulo>/index.html`.
4. Ejecutar `julia tools/pages/build.jl`.
5. Verificar `docs/index.html` y `docs/Pluto/<Modulo>/index.html`.
6. Ejecutar `julia tools/pages/publish.jl "mensaje"`.

## Nota sobre exportacion de Pluto

Este repositorio no fuerza una exportacion automatica de Pluto desde `build.jl`.
La exportacion HTML se considera un paso manual y verificable; `build.jl` solo sincroniza a `docs/` o genera placeholders cuando falta exportacion real.
