# Metodología EEG en NeuroSmart-EEG

## 1. Introducción
El análisis de EEG permite estudiar la dinámica eléctrica cerebral con alta resolución temporal. Este proyecto implementa un pipeline reproducible siguiendo estándares clínicos y recomendaciones de la literatura científica.

---

## 2. Flujo de Procesamiento
1. **Datos crudos (.vhdr, .eeg, .vmrk)**
2. **Preprocesamiento**
   - Filtrado Bandpass (0.5–45 Hz)
   - Filtro Notch (50 Hz)
   - Resampleo a 512 Hz
3. **Segmentación**
   - Epochs de 1 s
   - Corrección de línea base (0–100 ms)
4. **Rechazo de artefactos**
   - Umbral ±70 µV
   - Eliminación de epochs contaminados
5. **Análisis espectral**
   - FFT
   - Potencia absoluta y relativa en bandas:
     - Delta (0.5–4 Hz)
     - Theta (4–8 Hz)
     - Alpha (8–13 Hz)
     - Beta (13–30 Hz)
     - Gamma (>30 Hz)
6. **Conectividad funcional**
   - Coherencia
   - PLV
   - ImCoh
7. **Extracción de características**
   - Índices espectrales
   - Métricas de sincronización

---

## 3. Justificación de Parámetros
- **Filtro Bandpass (0.5–45 Hz):** Excluye artefactos lentos y ruido muscular.
- **Notch (50 Hz):** Elimina interferencia eléctrica.
- **Resampleo a 512 Hz:** Estándar en estudios clínicos y optimización FFT.
- **Baseline correction:** Reduce sesgo por variaciones lentas.
- **Artifact rejection ±70 µV:** Basado en rangos fisiológicos.
- **FFT:** Análisis en frecuencia robusto para biomarcadores.

---

## 4. Referencias
- Nunez PL, Srinivasan R (2006) *Electric Fields of the Brain*.
- Luck SJ (2014) *An Introduction to the ERP Technique*.
- Babiloni et al. (2016) *Brain connectivity in neurodegenerative diseases*.