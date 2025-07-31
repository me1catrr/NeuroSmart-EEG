#!/usr/bin/env julia
"""
Script: run_preprocessing.jl
Descripción: Ejecuta el pipeline de preprocesamiento completo para EEG.
"""

include("../src/preprocessing/preprocess_eeg.jl")

println("🚀 Ejecutando preprocesamiento completo...")
# Parámetros, rutas y logs están definidos en preprocess_eeg.jl