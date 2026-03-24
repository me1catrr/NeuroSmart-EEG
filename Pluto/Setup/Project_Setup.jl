### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 808ac11e-1109-4938-9297-b6e828890816
begin
using Markdown
using InteractiveUtils
using PlutoUI
end

# ╔═╡ 8c34c76f-5d77-44d5-ab73-c0de35540fec
PlutoUI.TableOfContents(title = "Contenido")

# ╔═╡ c586474d-9b33-4f1d-96fb-7f948724f4fe
md"
# Configuración del proyecto

Este notebook describe:

- la configuración del entorno
- la organización del proyecto
- las herramientas utilizadas 

para desarrollar el pipeline de análisis **EEG** en el lenguaje de programación **Julia**. En particular, está orientado a la configuración del proyecto, el entorno de trabajo y la organización del repositorio para mantener un flujo reproducible y claro en la plataforma de notebooks **Pluto.jl** (entorno **Julia**).

---

## Entorno reproducible

El proyecto se gestiona con un entorno de **Julia** definido por dos archivos:

- `Project.toml`: declara las dependencias del proyecto.
- `Manifest.toml`: fija las versiones exactas del entorno.

Mantener ambos archivos versionados garantiza ejecución consistente y resultados reproducibles entre sistemas.

---

# Julia

**Julia** se eligió como lenguaje de programación para este pipeline de análisis **EEG** por su equilibrio entre legibilidad y rendimiento en cómputo numérico intensivo. En este proyecto, ese equilibrio es especialmente útil para operaciones frecuentes como filtrado, transformadas de Fourier, álgebra lineal y métricas de conectividad.

### Rendimiento y eficiencia

El modelo **JIT** (Just In Time) basado en **LLVM** permite ejecutar código científico con rendimiento cercano a lenguajes compilados, reduciendo la necesidad de depender de extensiones externas para partes críticas.

Esto resulta particularmente útil en tareas comunes en análisis **EEG** como:

- transformadas de Fourier  
- filtrado digital  
- operaciones matriciales de gran tamaño  
- cálculo de métricas de conectividad entre múltiples canales  
- simulaciones o generación de datos surrogate

Además, su gestor de paquetes y el uso de `Project.toml` + `Manifest.toml` facilitan reproducibilidad estricta del entorno en investigación.
"

# ╔═╡ c34fe44f-7595-46bf-b888-0624aae68bab
md"
---

# Pluto.jl

**Pluto.jl** se utiliza como interfaz principal para documentar y explorar el pipeline de análisis **EEG**. Su modelo reactivo mantiene consistencia entre celdas y evita estados desalineados durante el trabajo interactivo.

---

### ¿Qué es Pluto.jl?

**Pluto.jl** permite integrar en un mismo documento:

- código ejecutable  
- texto explicativo (Markdown)  
- ecuaciones matemáticas  
- visualizaciones y resultados  

en un notebook dinámico orientado a trazabilidad y comunicación técnica. **Pluto** detecta dependencias entre celdas y recalcula automáticamente el grafo afectado, manteniendo un estado coherente.

---

### Ejecución reactiva y coherencia del estado

Frente a flujos de ejecución manual por celdas, **Pluto** mantiene un modelo reactivo basado en dependencias:

- cada variable se define una sola vez  
- las dependencias entre celdas se detectan automáticamente  
- los cambios en una celda actualizan todas las celdas dependientes  

Esto mejora reproducibilidad, reduce errores por orden de ejecución y facilita la revisión del pipeline de análisis **EEG**.

---

## Cómo iniciar Pluto.jl

Para trabajar con **Pluto** es necesario tener **Julia** instalado y el paquete `Pluto` añadido al entorno.  Si es la primera vez que se utiliza, puede instalarse desde el gestor de paquetes de **Julia**:

```julia
using Pkg
Pkg.add(\"Pluto\")
```

Una vez instalado, **Pluto** se inicia desde el intérprete de **Julia**:

```julia
using Pluto
Pluto.run()
```

Al ejecutar este comando, **Pluto** abre un servidor local y la interfaz en el navegador.

Desde esta interfaz se pueden:

- crear notebooks nuevos
- abrir notebooks existentes
- gestionar el historial de notebooks recientes

---

## Abrir un notebook existente

Para abrir un notebook existente, puede seleccionarse el archivo `.jl` desde la interfaz de **Pluto**.

También es posible abrirlo directamente desde Julia:

```julia
using Pluto
Pluto.run(notebook=\"ruta/al/notebook.jl\")
```

**Pluto** cargará el notebook y ejecutará las celdas necesarias respetando sus dependencias.

---

## Cerrar Pluto.jl

Para cerrar Pluto:

- **Desde la interfaz web:** utilizando el botón *Shutdown* del servidor.
- **Desde el terminal o REPL de Julia:** presionando `Ctrl + C`.

Al detener el servidor, las sesiones activas se cierran y el navegador se desconecta del entorno de ejecución.

---

## Uso de Pluto en este proyecto

En este proyecto, **Pluto** se utiliza como **documentación interactiva y entorno de exploración** del pipeline de análisis **EEG**. Esto permite:

- describir cada fase del procesamiento **EEG** junto con el código correspondiente
- visualizar resultados intermedios y finales
- mantener una estructura reproducible del análisis

Además, su integración con **Julia** permite aprovechar de forma directa:

- bibliotecas de computación científica  
- procesamiento de señales  
- álgebra lineal de alto rendimiento  
- visualización avanzada

Esto resulta útil en el pipeline de análisis **EEG**, que incluye:

- procesamiento de grandes matrices de datos  
- transformadas de Fourier  
- filtrado de señales  
- cálculo de métricas de conectividad entre múltiples canales

---

### Beneficios para investigación reproducible

El uso de **Pluto.jl** aporta varias ventajas importantes para proyectos de investigación computacional:

- **reproducibilidad garantizada** del entorno y del flujo de ejecución  
- **transparencia del análisis**, al integrar código y explicación en el mismo documento  
- **facilidad para compartir notebooks completos** con otros investigadores  
- **exploración interactiva de datos y resultados**

En conjunto, **Pluto** permite que el pipeline de análisis EEG desarrollado en este proyecto sea **fácil de entender, reproducir y extender por otros investigadores**.

---

Una vez definido el entorno de trabajo, es necesario organizar el código y los datos de forma estructurada para mantener un pipeline claro y reproducible.
"

# ╔═╡ 3373ef23-bd16-46a7-8434-3c2caa851939
md"
# Estructura del proyecto

El código del proyecto está organizado como un único entorno de **Julia** denominado `EEG_JULIA`. Este entorno contiene todas las dependencias necesarias para ejecutar el pipeline completo de análisis.

La estructura general separa claramente los componentes del flujo de trabajo:

```text
EEG_JULIA/
├── config/        # configuración del pipeline
├── data/          # datos de entrada e intermedios
├── results/       # resultados generados
├── src/           # lógica principal del pipeline
├── script/        # ejecución automatizada
├── Pluto/         # notebooks interactivos
├── Project.toml
└── Manifest.toml
```

- **configuración**  
- **datos de entrada**  
- **resultados generados**  
- **módulos de código fuente**

Roles principales:

- `src/` concentra la lógica reutilizable del pipeline de análisis **EEG**.
- `script/` define puntos de entrada y ejecución automatizada.
- `Pluto/` se utiliza para documentación, exploración y visualización.

Esta separación de responsabilidades facilita mantenimiento y trazabilidad: los notebooks se enfocan en explicación y análisis interactivo, mientras `src/` y `script/` mantienen la ejecución reproducible del pipeline de análisis EEG.

Esta organización facilita la **reproducibilidad del análisis**, permitiendo ejecutar el pipeline completo a partir del mismo entorno de **Julia** y garantizando que todas las dependencias se encuentren en versiones controladas.

| Directorio | Contenido y propósito |
|---|---|
| `config/` | Archivos de configuración como `default_config.jl`: rutas, listas de sujetos y parámetros que controlan el pipeline sin necesidad de modificar el código. |
| `data/` | Datos de entrada e intermedios para cada etapa del procesamiento (por ejemplo `raw/`, `filtering/`, `ICA/`, `segmentation/`, `FFT/`, `CSD/`, `wPLI/`). Mantiene la organización del disco alineada con las etapas del pipeline de análisis EEG. |
| `results/` | Resultados generados agrupados por tipo: `figures/` (CSD, FFT, ICA_cleaning, wPLI), `logs/` y `tables/`. Incluye todo aquello que puede regenerarse a partir del código y los datos. |
| `src/` | Módulos principales de Julia, generalmente un archivo `.jl` por cada etapa del pipeline (IO, filtrado, ICA, segmentación, FFT, baseline, rechazo de artefactos y módulos de conectividad como `CSD.jl` o `wPLI.jl`). |
| `script/` | Scripts ejecutables (por ejemplo drivers del pipeline o scripts de análisis) que integran los módulos definidos en `src/`. |

---

## Generación de figuras (Makie)

Las figuras más complejas de resultados (matrices de conectividad, espectros, comparaciones entre grupos) se generan directamente en **Julia** utilizando (entre otras) la librería **Makie** y se exportan como archivos **PDF vectoriales** para su inclusión directa en este informe.

**Makie** fue elegida además de librerías como **Plots.jl** por varias razones técnicas.

**Plots.jl** proporciona una API unificada sobre múltiples backends gráficos (GR, PyPlot, entre otros) y es muy conveniente para visualización rápida o exploratoria. Sin embargo, cuando se trabaja con:

- matrices grandes (por ejemplo matrices canal × canal de conectividad),
- series temporales largas,
- o figuras complejas con múltiples paneles,

su rendimiento puede ser limitado y el comportamiento puede variar entre backends.

**Makie**, en cambio, es una **arquitectura gráfica moderna basada en un único stack**, lo que garantiza mayor consistencia y rendimiento. Sus principales ventajas son:

- **renderizado acelerado por GPU**
- **manejo eficiente de grandes matrices**
- **salida vectorial de alta calidad (PDF o SVG)**
- **control fino del diseño de figuras multi-panel**

Además, **Makie** utiliza un sistema de **scene graph y layouts jerárquicos**, que facilita la construcción de figuras complejas, como por ejemplo:

- una matriz de conectividad por banda de frecuencia  
- comparaciones entre grupos experimentales  
- paneles múltiples con escalas compartidas  

Para este proyecto se utiliza el backend **CairoMakie**, optimizado para la generación de figuras estáticas en formato vectorial.

---

## Ejemplo mínimo de generación de figura

El siguiente ejemplo muestra el flujo básico para crear una figura y exportarla a PDF.

```julia
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel = \"Channel\",
    ylabel = \"Channel\"
)

heatmap!(ax, rand(8,8))

save(\"connectivity-schematic.pdf\", fig)
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
|  | `git commit -m \"mensaje\"` | Crea un snapshot de los cambios preparados con un mensaje descriptivo. |
|  | `git push -u origin main` | Sube la rama `main` a GitHub y la establece como rama remota por defecto. |
| Clonar un repositorio existente de GitHub | `git clone <url>` | Descarga el repositorio remoto en una nueva carpeta local. |
|  | `cd <folder>` | Cambia al directorio del proyecto clonado. |
|  | `julia -e \"using Pkg; Pkg.instantiate()\"` | Instala el entorno de Julia definido por `Project.toml` y `Manifest.toml`. |
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

```text
data/
results/
```

se excluyen del repositorio mediante el archivo:

```text
.gitignore
```

Este archivo indica a Git qué archivos o carpetas deben ignorarse durante el seguimiento de cambios.

| Qué incluir en el repositorio | Por qué |
|---|---|
| `src/` | Código fuente principal (módulos y funciones del pipeline). |
| `script/` | Scripts ejecutables utilizados para lanzar o coordinar el pipeline de análisis EEG. |
| `config/` | Archivos de configuración que definen parámetros y rutas de forma portable. |
| `Project.toml` | Declara las dependencias del proyecto para que otros puedan reproducir el entorno. |
| `Manifest.toml` | Fija las versiones exactas de las dependencias y garantiza reproducibilidad estricta entre máquinas. |

### `.gitignore` recomendado

Un archivo `.gitignore` mínimo para un proyecto en **Julia** debería excluir artefactos generados automáticamente, configuraciones locales del editor y salidas de datos grandes.

| Patrón / ruta | Por qué ignorarlo |
|---|---|
| `*.jl.cov` | Archivos de cobertura generados durante las pruebas. |
| `data/` | Datos brutos o intermedios, que suelen ser grandes y no adecuados para control de versiones. |
| `results/` | Resultados generados automáticamente (figuras, exportaciones o cálculos cacheados). |
| `.DS_Store` | Archivos de metadatos generados por macOS. |
| `.vscode/`, `.cursor/` | Configuración local de editores o IDE (ignorar salvo que quieras compartir configuración de equipo). |

Ignorar estos archivos ayuda a mantener el repositorio pequeño y evita commits accidentales de archivos grandes o temporales.

En este proyecto se recomienda versionar siempre `Project.toml` y `Manifest.toml` para asegurar reproducibilidad exacta, especialmente en flujos científicos.

---

## Manejo de archivos grandes

En los casos en que ciertos archivos de datos deban mantenerse versionados (por ejemplo, modelos preentrenados o datasets pequeños pero importantes), puede utilizarse **Git Large File Storage (Git LFS)**.

Git LFS permite gestionar archivos grandes almacenando en el repositorio únicamente referencias ligeras, mientras que el contenido real se mantiene en un almacenamiento externo optimizado para este tipo de datos.

---

## Beneficios para la reproducibilidad científica

El uso de **Git** y **GitHub** contribuye de forma significativa a la **reproducibilidad del proyecto**:

- permite reconstruir exactamente qué versión del código produjo determinados resultados  
- facilita la revisión y auditoría del pipeline de análisis  
- simplifica la distribución del proyecto a otros investigadores  

En conjunto, el control de versiones garantiza que el desarrollo del pipeline de análisis EEG sea **transparente, trazable y reproducible**.
"

# ╔═╡ 749a5ee5-b336-43e7-9cea-fbe323eb478b
md"""
# Herramientas de apoyo al desarrollo

Algunas herramientas de **IA** se utilizaron como apoyo técnico para acelerar tareas de desarrollo y documentación, siempre con revisión manual del contenido generado.

---

## Cursor

El editor **Cursor** (<https://cursor.com>) se utilizó como entorno de desarrollo para el código en **Julia**. Principales usos:

- escritura inicial de fragmentos de código  
- refactorización y mejora de funciones existentes  
- ayuda en la depuración (*debugging*)  
- navegación y comprensión de la estructura del proyecto  
- generación y ajuste de fragmentos de documentación  

---

## ChatGPT

**ChatGPT (OpenAI)** se utilizó para apoyo conceptual y redacción técnica en etapas no críticas de implementación:

- clarificación de conceptos metodológicos (por ejemplo **BIDS**, **wPLI**, métodos **surrogate**, o aspectos de procesamiento de **EEG**)  
- generación de ideas preliminares para estructurar secciones del informe  
- redacción inicial de resúmenes o explicaciones técnicas  
- sugerencias de organización para la documentación del pipeline  

Las respuestas se utilizaron como borrador de trabajo y se validaron manualmente antes de su incorporación.

---

## Consideraciones de reproducibilidad

Para mantener trazabilidad y reproducibilidad del proyecto:

- las herramientas de **IA** se usaron como asistencia técnica  
- el diseño metodológico y la implementación final se verificaron manualmente  
- la documentación final se alinea con el código ejecutable del pipeline  

El objetivo fue mejorar eficiencia y claridad sin introducir cambios automáticos no revisados.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.80"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "497f6aa4b56a2c25732a286a0a968050c5860834"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "fbc875044d82c113a9dee6fc14e16cf01fd48872"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.80"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"
"""

# ╔═╡ Cell order:
# ╠═808ac11e-1109-4938-9297-b6e828890816
# ╠═8c34c76f-5d77-44d5-ab73-c0de35540fec
# ╠═c586474d-9b33-4f1d-96fb-7f948724f4fe
# ╠═c34fe44f-7595-46bf-b888-0624aae68bab
# ╟─3373ef23-bd16-46a7-8434-3c2caa851939
# ╟─0aa3b638-ab24-4fbb-8dfc-f5bd5e931f0c
# ╠═749a5ee5-b336-43e7-9cea-fbe323eb478b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
