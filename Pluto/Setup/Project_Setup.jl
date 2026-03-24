### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ c586474d-9b33-4f1d-96fb-7f948724f4fe
md"
# Implementación en Julia

## ¿Por qué Julia?

**Julia** es un lenguaje de programación de alto nivel diseñado específicamente para **computación científica y numérica de alto rendimiento**. Fue desarrollado con el objetivo de combinar la facilidad de uso de lenguajes dinámicos con la velocidad de los lenguajes compilados tradicionales.

En proyectos de análisis de datos neurofisiológicos —como el procesamiento de EEG y el análisis de conectividad— es habitual enfrentarse a grandes volúmenes de datos, transformaciones numéricas intensivas y pipelines complejos. En este contexto, Julia ofrece una combinación especialmente adecuada de **legibilidad, rendimiento y capacidad de desarrollo científico**.

### Rendimiento y eficiencia

Una de las características más destacadas de Julia es que su código se compila mediante **Just-In-Time (JIT compilation)** usando LLVM. Esto permite que muchas operaciones numéricas se ejecuten a velocidades comparables a lenguajes compilados como **C o Fortran**, sin sacrificar la simplicidad de un lenguaje de alto nivel.

En lenguajes como Python, alcanzar este nivel de rendimiento suele requerir escribir partes críticas del código en C, Cython o utilizar librerías externas altamente optimizadas. En Julia, por el contrario, es posible escribir directamente código científico de alto rendimiento en el propio lenguaje.

Esto resulta particularmente útil en tareas comunes en análisis EEG como:

- transformadas de Fourier  
- filtrado digital  
- operaciones matriciales de gran tamaño  
- cálculo de métricas de conectividad entre múltiples canales  
- simulaciones o generación de datos surrogate

### Comparación con Python

Python es actualmente uno de los lenguajes más utilizados en ciencia de datos y neurociencia computacional gracias a su gran ecosistema de librerías. Sin embargo, Julia ofrece varias ventajas relevantes:

- **Mayor rendimiento en código numérico puro**, evitando la necesidad de extensiones en C.
- **Tipado múltiple y despacho múltiple**, lo que permite escribir código genérico y altamente optimizado.
- **Sintaxis matemática clara**, especialmente adecuada para expresar algoritmos científicos.

En términos de experiencia de usuario, Julia mantiene muchas de las ventajas de Python:

- sintaxis legible  
- uso interactivo en notebooks  
- ecosistema creciente de paquetes científicos  

pero con un rendimiento significativamente mayor en muchos casos.

### Comparación con MATLAB

MATLAB ha sido tradicionalmente una herramienta muy utilizada en procesamiento de señales y neurociencia. Julia ofrece un entorno muy similar en cuanto a estilo de programación:

- operaciones vectorizadas y matriciales nativas  
- álgebra lineal integrada en el lenguaje  
- sintaxis matemática compacta  

Sin embargo, Julia presenta varias ventajas importantes:

- **software completamente libre y open source**
- **rendimiento comparable a lenguajes compilados**
- **gestión moderna de paquetes y entornos reproducibles**
- **mejor integración con herramientas modernas de computación científica**

En resumen, Julia combina características que tradicionalmente estaban separadas entre distintos lenguajes:  
es **tan fácil de usar como Python**, **tan natural para cálculo matricial como MATLAB**, y **tan rápido como C**.

---

## Pluto.jl como entorno 

El desarrollo, documentación y exploración interactiva del pipeline de análisis EEG se realizaron utilizando **Pluto.jl**, un entorno de notebooks reactivo diseñado específicamente para el lenguaje **Julia**. Pluto combina características de programación interactiva, documentación científica y reproducibilidad computacional en una única interfaz ligera.

A diferencia de los notebooks tradicionales, Pluto está diseñado desde el principio para garantizar **consistencia, trazabilidad y ejecución reproducible**, lo que lo convierte en una herramienta especialmente adecuada para proyectos de investigación científica.

---

### ¿Qué es Pluto.jl?

**Pluto.jl** es un sistema de notebooks interactivos basado en Julia que permite combinar:

- código ejecutable  
- texto explicativo (Markdown)  
- ecuaciones matemáticas  
- visualizaciones y resultados  

en un único documento dinámico.

Cada celda del notebook contiene una única expresión de Julia, y Pluto gestiona automáticamente las dependencias entre celdas. Cuando una variable cambia, **todas las celdas que dependen de ella se recalculan automáticamente**, manteniendo el notebook siempre en un estado consistente.

Este modelo reactivo elimina problemas comunes en otros entornos de notebooks, como el uso de variables obsoletas o ejecuciones fuera de orden.

---

### Ventajas frente a Jupyter Notebooks

Aunque **Jupyter Notebook** es una herramienta ampliamente utilizada en ciencia de datos, Pluto introduce varias mejoras importantes para el trabajo científico reproducible.

#### Ejecución reactiva y coherencia del estado

En Jupyter, las celdas pueden ejecutarse en cualquier orden, lo que puede generar inconsistencias entre el estado real del kernel y el código visible en el notebook. Esto puede dificultar la reproducción exacta de resultados.

Pluto evita este problema mediante un **modelo reactivo basado en dependencias**:

- cada variable solo puede definirse una vez  
- las dependencias entre celdas se detectan automáticamente  
- los cambios en una celda actualizan todas las celdas dependientes  

Esto garantiza que el notebook siempre represente un **estado válido y reproducible del análisis**.

---

## Cómo iniciar Pluto.jl

Para trabajar con Pluto es necesario tener **Julia instalado** y el paquete `Pluto` añadido al entorno.  
Si es la primera vez que se utiliza, puede instalarse desde el gestor de paquetes de Julia:

```
using Pkg
Pkg.add('Pluto')
```

Una vez instalado, Pluto se inicia desde el intérprete de Julia ejecutando:

```
using Pluto
Pluto.run()
```

Al ejecutar este comando, Pluto lanza automáticamente un **servidor local** y abre una pestaña en el navegador web con la interfaz de notebooks.

Desde esta interfaz se pueden:

- crear notebooks nuevos
- abrir notebooks existentes
- gestionar el historial de notebooks recientes

---

## Abrir un notebook existente

Para abrir un notebook previamente creado (por ejemplo uno del pipeline de análisis EEG), basta con seleccionar el archivo `.jl` correspondiente desde la interfaz de Pluto.

También es posible abrirlo directamente desde Julia:

```
using Pluto
Pluto.run(notebook='ruta/al/notebook.jl')
```

Pluto cargará el notebook y ejecutará automáticamente las celdas necesarias respetando las dependencias entre ellas.

---

## Cerrar Pluto.jl

Para cerrar Pluto existen varias opciones:

- **Desde la interfaz web:** utilizando el botón *Shutdown* del servidor.
- **Desde el terminal o REPL de Julia:** presionando `Ctrl + C`.

Al detener el servidor, todas las sesiones activas de Pluto se cierran y el navegador deja de estar conectado al entorno de ejecución.

---

## Ventajas para el desarrollo del pipeline

En este proyecto, Pluto se utiliza no solo como herramienta de ejecución sino también como **documentación interactiva del pipeline**. Esto permite:

- describir cada fase del procesamiento EEG junto con el código correspondiente
- visualizar resultados intermedios y finales
- mantener una estructura reproducible del análisis

De este modo, el notebook actúa simultáneamente como **documentación, entorno de experimentación y registro del análisis computacional**.

---

#### Reproducibilidad automática del entorno

Pluto integra directamente el sistema de gestión de paquetes de Julia. Cada notebook puede almacenar internamente la información sobre las dependencias utilizadas en el análisis.

Esto permite que otro usuario pueda abrir el notebook y reproducir automáticamente el entorno de trabajo sin necesidad de instalar manualmente librerías específicas.

Este enfoque reduce significativamente los problemas habituales de reproducibilidad asociados con diferencias de versiones de paquetes o configuraciones del sistema.

---

#### Integración natural con Julia

Pluto está diseñado específicamente para el ecosistema Julia, lo que permite aprovechar de forma directa:

- bibliotecas de computación científica  
- procesamiento de señales  
- álgebra lineal de alto rendimiento  
- visualización avanzada  

Esto resulta especialmente útil en pipelines de análisis de EEG que incluyen:

- procesamiento de grandes matrices de datos  
- transformadas de Fourier  
- filtrado de señales  
- cálculo de métricas de conectividad entre múltiples canales  

El alto rendimiento de Julia permite realizar estas operaciones directamente en el notebook sin necesidad de optimizaciones externas.

---

#### Documentación científica integrada

Pluto facilita la creación de documentos científicos interactivos que combinan:

- explicación teórica  
- código reproducible  
- resultados visuales  

Esto permite utilizar el notebook no solo como herramienta de análisis, sino también como **documentación ejecutable del pipeline completo**.

En el contexto de este proyecto, Pluto se utiliza para:

- documentar cada etapa del pipeline EEG  
- visualizar resultados intermedios y finales  
- explicar los fundamentos metodológicos del análisis  
- generar figuras utilizadas en el informe

---

### Beneficios para la investigación reproducible

El uso de Pluto.jl aporta varias ventajas importantes para proyectos de investigación computacional:

- **reproducibilidad garantizada** del entorno y del flujo de ejecución  
- **transparencia del análisis**, al integrar código y explicación en el mismo documento  
- **facilidad para compartir notebooks completos** con otros investigadores  
- **exploración interactiva de datos y resultados**

En conjunto, Pluto permite que el pipeline de análisis EEG desarrollado en este proyecto sea **fácil de entender, reproducir y extender por otros investigadores**.

---

## Estructura y organización del proyecto

El código del proyecto está organizado como un único entorno de Julia denominado `EEG_JULIA`. Este entorno contiene todas las dependencias necesarias para ejecutar el pipeline completo de análisis.

La estructura general del proyecto separa claramente los distintos componentes del flujo de trabajo:

- **configuración**  
- **datos de entrada**  
- **resultados generados**  
- **módulos de código fuente**

Esta organización facilita la **reproducibilidad del análisis**, permitiendo ejecutar el pipeline completo a partir del mismo entorno de Julia y garantizando que todas las dependencias se encuentren en versiones controladas.

| Directorio | Contenido y propósito |
|---|---|
| `config/` | Archivos de configuración como `default_config.jl`: rutas, listas de sujetos y parámetros que controlan el pipeline sin necesidad de modificar el código. |
| `data/` | Datos de entrada e intermedios para cada etapa del procesamiento (por ejemplo `raw/`, `filtering/`, `ICA/`, `segmentation/`, `FFT/`, `CSD/`, `wPLI/`). Mantiene la organización del disco alineada con las etapas del pipeline. |
| `results/` | Resultados generados agrupados por tipo: `figures/` (CSD, FFT, ICA_cleaning, wPLI), `logs/` y `tables/`. Incluye todo aquello que puede regenerarse a partir del código y los datos. |
| `src/` | Módulos principales de Julia, generalmente un archivo `.jl` por cada etapa del pipeline (IO, filtrado, ICA, segmentación, FFT, baseline, rechazo de artefactos y módulos de conectividad como `CSD.jl` o `wPLI.jl`). |
| `script/` | Scripts ejecutables (por ejemplo drivers del pipeline o scripts de análisis) que integran los módulos definidos en `src/`. |

---

## Generación de figuras (Makie)

Todas las figuras de resultados (matrices de conectividad, espectros, comparaciones entre grupos) se generan directamente en Julia utilizando la librería **Makie** y se exportan como archivos **PDF vectoriales** para su inclusión directa en este informe.

Makie fue elegida en lugar de librerías como **Plots.jl** por varias razones técnicas.

Plots.jl proporciona una API unificada sobre múltiples backends gráficos (GR, PyPlot, entre otros) y es muy conveniente para visualización rápida o exploratoria. Sin embargo, cuando se trabaja con:

- matrices grandes (por ejemplo matrices canal × canal de conectividad),
- series temporales largas,
- o figuras complejas con múltiples paneles,

su rendimiento puede ser limitado y el comportamiento puede variar entre backends.

Makie, en cambio, es una **arquitectura gráfica moderna basada en un único stack**, lo que garantiza mayor consistencia y rendimiento. Sus principales ventajas son:

- **renderizado acelerado por GPU**
- **manejo eficiente de grandes matrices**
- **salida vectorial de alta calidad (PDF o SVG)**
- **control fino del diseño de figuras multi-panel**

Además, Makie utiliza un sistema de **scene graph y layouts jerárquicos**, que facilita la construcción de figuras complejas, como por ejemplo:

- una matriz de conectividad por banda de frecuencia  
- comparaciones entre grupos experimentales  
- paneles múltiples con escalas compartidas  

Para este proyecto se utiliza el backend **CairoMakie**, optimizado para la generación de figuras estáticas en formato vectorial.

---

## Ejemplo mínimo de generación de figura

El siguiente ejemplo muestra el flujo básico para crear una figura y exportarla a PDF.

```
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1,1],
    xlabel = 'Channel',
    ylabel = 'Channel'
)

heatmap!(ax, rand(8,8))

save('connectivity-schematic.pdf', fig)
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
|  | `git commit -m 'mensaje'` | Crea un snapshot de los cambios preparados con un mensaje descriptivo. |
|  | `git push -u origin main` | Sube la rama `main` a GitHub y la establece como rama remota por defecto. |
| Clonar un repositorio existente de GitHub | `git clone <url>` | Descarga el repositorio remoto en una nueva carpeta local. |
|  | `cd <folder>` | Cambia al directorio del proyecto clonado. |
|  | `julia -e 'using Pkg; Pkg.instantiate()'` | Instala el entorno de Julia definido por `Project.toml` y `Manifest.toml`. |
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

```
data/
results/
```

se excluyen del repositorio mediante el archivo:

```
.gitignore
```

Este archivo indica a Git qué archivos o carpetas deben ignorarse durante el seguimiento de cambios.

| Qué incluir en el repositorio | Por qué |
|---|---|
| `src/` | Código fuente principal (módulos y funciones del pipeline). |
| `script/` | Scripts ejecutables utilizados para lanzar o coordinar el pipeline. |
| `config/` | Archivos de configuración que definen parámetros y rutas de forma portable. |
| `Project.toml` | Declara las dependencias del proyecto para que otros puedan reproducir el entorno. |
| `Manifest.toml` *(opcional)* | Fija las versiones exactas de las dependencias para garantizar reproducibilidad completa entre máquinas. |

### `.gitignore` recomendado

Un archivo `.gitignore` mínimo para un proyecto en **Julia** debería excluir artefactos generados automáticamente, configuraciones locales del editor y salidas de datos grandes.

| Patrón / ruta | Por qué ignorarlo |
|---|---|
| `*.jl.cov` | Archivos de cobertura generados durante las pruebas. |
| `Manifest.toml` | Opcional: se puede ignorar si se prefiere un entorno flexible sin versiones bloqueadas (inclúyelo si quieres reproducibilidad completa). |
| `data/` | Datos brutos o intermedios, que suelen ser grandes y no adecuados para control de versiones. |
| `results/` | Resultados generados automáticamente (figuras, exportaciones o cálculos cacheados). |
| `.DS_Store` | Archivos de metadatos generados por macOS. |
| `.vscode/`, `.cursor/` | Configuración local de editores o IDE (ignorar salvo que quieras compartir configuración de equipo). |

Ignorar estos archivos ayuda a mantener el repositorio pequeño y evita commits accidentales de archivos grandes o temporales.

---

## Manejo de archivos grandes

En los casos en que ciertos archivos de datos deban mantenerse versionados (por ejemplo, modelos preentrenados o datasets pequeños pero importantes), puede utilizarse **Git Large File Storage (Git LFS)**.

Git LFS permite gestionar archivos grandes almacenando en el repositorio únicamente referencias ligeras, mientras que el contenido real se mantiene en un almacenamiento externo optimizado para este tipo de datos.

---

## Beneficios para la reproducibilidad científica

El uso de Git y GitHub contribuye de forma significativa a la **reproducibilidad del proyecto**:

- permite reconstruir exactamente qué versión del código produjo determinados resultados  
- facilita la revisión y auditoría del pipeline de análisis  
- simplifica la distribución del proyecto a otros investigadores  

En conjunto, el control de versiones garantiza que el desarrollo del pipeline de análisis EEG sea **transparente, trazable y reproducible**.
"

# ╔═╡ 749a5ee5-b336-43e7-9cea-fbe323eb478b
md"""
# Uso responsable de herramientas IA

Durante la elaboración de este informe y el desarrollo del código en **Julia**, se utilizaron herramientas de **inteligencia artificial (IA)** como apoyo al proceso de trabajo. Estas herramientas se emplearon de forma **transparente, responsable y supervisada**, con el objetivo de mejorar la productividad, facilitar la exploración de ideas y apoyar la redacción técnica.

Es importante destacar que las herramientas de IA se utilizaron **como asistencia**, no como sustituto del razonamiento científico ni del trabajo personal. Todas las decisiones metodológicas, la implementación del pipeline y la interpretación de los resultados corresponden al autor del proyecto.

---

## Cursor

El editor **Cursor** (<https://cursor.com>) fue utilizado como entorno de desarrollo para el código en Julia. Cursor integra un asistente basado en inteligencia artificial que permite interactuar con el código del proyecto de forma contextual.

Este entorno facilitó varias tareas durante el desarrollo del pipeline:

- escritura inicial de fragmentos de código  
- refactorización y mejora de funciones existentes  
- ayuda en la depuración (*debugging*)  
- navegación y comprensión de la estructura del proyecto  
- generación o adaptación de pequeños fragmentos de documentación  

El asistente integrado permitió acelerar tareas repetitivas o mecánicas, manteniendo siempre la revisión manual del código generado.

---

## ChatGPT

**ChatGPT (OpenAI)** se utilizó como herramienta de apoyo conceptual y de redacción durante el desarrollo del proyecto.

Entre los usos principales se incluyen:

- clarificación de conceptos metodológicos (por ejemplo **BIDS**, **wPLI**, métodos **surrogate**, o aspectos de procesamiento de EEG)  
- generación de ideas preliminares para estructurar secciones del informe  
- redacción inicial de resúmenes o explicaciones técnicas  
- sugerencias de organización para la documentación del pipeline  

Las respuestas generadas por ChatGPT se utilizaron únicamente como **punto de partida para el trabajo personal**. Todo el contenido final fue revisado, editado y adaptado manualmente para garantizar su exactitud científica y coherencia con el proyecto.

---

## Consideraciones éticas y de reproducibilidad

El uso de herramientas de inteligencia artificial se realizó siguiendo principios de **transparencia y responsabilidad académica**. En particular:

- las herramientas de IA se emplearon como apoyo, no como sustituto del análisis científico  
- el diseño metodológico y la implementación final del código fueron realizados y verificados manualmente  
- las descripciones del pipeline se basan en el código y en la comprensión directa del proceso de análisis  

De este modo, el uso de IA contribuyó principalmente a **mejorar la claridad y eficiencia del proceso de desarrollo**, sin comprometer la integridad científica del trabajo.
"""

# ╔═╡ Cell order:
# ╟─c586474d-9b33-4f1d-96fb-7f948724f4fe
# ╟─0aa3b638-ab24-4fbb-8dfc-f5bd5e931f0c
# ╟─749a5ee5-b336-43e7-9cea-fbe323eb478b
