# Tools Pages

Tooling de publicacion para GitHub Pages en este repositorio.

## Flujo completo de publicacion

1. Construir/actualizar contenido de Pages (`docs/`) con placeholders o exportaciones reales.
2. Revisar cambios de git.
3. Hacer commit.
4. Hacer push al remoto.
5. Verificar que GitHub Pages se actualiza correctamente.

Este flujo se automatiza con `tools/pages/publish.jl`.

## Comandos de uso

Ejecucion con mensaje automatico:

```bash
julia tools/pages/publish.jl
```

Ejecucion con mensaje personalizado:

```bash
julia tools/pages/publish.jl "Update pages after Pluto export"
```

## Que hace internamente `publish.jl`

El script ejecuta estos pasos:

1. Lanza `tools/pages/build.jl`.
2. Consulta `git status --short` y muestra archivos modificados.
3. Ejecuta `git add .`.
4. Hace commit:
   - mensaje por defecto: `Update GitHub Pages (auto build)`,
   - o mensaje personalizado si se pasa por argumento.
5. Ejecuta `git push`.

## Comportamiento ante casos comunes

- Si **no hay cambios**, informa y termina sin error.
- Si falla un comando de git o build, muestra un error claro y detiene el proceso.
- Muestra en consola una traza amigable de cada paso para facilitar debug.

## Verificacion de GitHub Pages

Tras `git push`, espera 1-3 minutos y valida:

- que `docs/index.html` refleja la version esperada;
- que los enlaces de modulos `docs/Pluto/.../index.html` cargan correctamente;
- que `docs/pages_manifest.json` se actualizo con el estado correcto.
