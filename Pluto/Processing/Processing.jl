### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 2ff7ee54-7bec-4e21-b804-e9fcc2f86a11
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido - Processing")
end

# ╔═╡ 49ca7fca-11dc-4f81-b74f-b54e8b252e10
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	using .PlutoTemplateBase
	using Serialization
end

# ╔═╡ c8aaf8df-4f1b-4c07-84fd-8a5eb18ab0a7
notebook_intro("PROCESSING")

# ╔═╡ 01b14950-89b8-4df1-981f-44a13f6026cb
md"""
Notebook por fase extraido de `Pluto/Notebook.jl`, centrado en:
- ICA / ICA_cleaning,
- segmentation / baseline / artifact_rejection / baseline_2st,
- FFT.
"""

# ╔═╡ 4b9838dc-bb6b-4bcc-8f51-e30f6f89be29
begin
	# Rutas por etapa sin hardcodeo
	dir_filtering = stage_dir(:filtering)
	dir_ica = stage_dir(:ICA)
	dir_ica_cleaning = stage_dir(:ICA_cleaning)
	dir_segmentation = stage_dir(:segmentation)
	dir_baseline = stage_dir(:baseline)
	dir_artifact = stage_dir(:artifact_rejection)
	dir_baseline_2st = stage_dir(:baseline_2st)
	dir_fft = stage_dir(:FFT)

	(
		dir_filtering = dir_filtering,
		dir_ica = dir_ica,
		dir_ica_cleaning = dir_ica_cleaning,
		dir_segmentation = dir_segmentation,
		dir_baseline = dir_baseline,
		dir_artifact = dir_artifact,
		dir_baseline_2st = dir_baseline_2st,
		dir_fft = dir_fft,
	)
end

# ╔═╡ 4d8d2716-9c55-4f39-98f0-a08bbfb661eb
begin
	# Equivalente al bloque de carga usado en la seccion ICA del notebook global
	path_dict_lowpass = joinpath(stage_dir(:filtering), "dict_EEG_Lowpass.bin")
	dict_EEG_filtered = if isfile(path_dict_lowpass)
		Serialization.deserialize(path_dict_lowpass)
	else
		@warn "No se encontro dict_EEG_Lowpass.bin; usando demo." path_dict_lowpass
		Dict("Cz" => randn(2000), "Pz" => randn(2000), "Fz" => randn(2000), "Oz" => randn(2000))
	end
	dict_EEG_filtered
end

# ╔═╡ 4b069ffe-1b16-492d-acf3-83018ef16a43
begin
	channels_ICA = sort(collect(keys(dict_EEG_filtered)))
	n_channels = length(channels_ICA)
	n_samples = n_channels > 0 ? length(dict_EEG_filtered[first(channels_ICA)]) : 0
	(
		n_channels = n_channels,
		n_samples = n_samples,
		channels_ICA = channels_ICA,
	)
end

# ╔═╡ 5353f914-a0eb-4e3b-bf5f-cbe90c63df88
md"""
## Salidas esperadas por fase (Processing)

- `ICA`: `$(stage_dir(:ICA))`
- `ICA_cleaning`: `$(stage_dir(:ICA_cleaning))`
- `segmentation`: `$(stage_dir(:segmentation))`
- `baseline`: `$(stage_dir(:baseline))`
- `artifact_rejection`: `$(stage_dir(:artifact_rejection))`
- `baseline_2st`: `$(stage_dir(:baseline_2st))`
- `FFT`: `$(stage_dir(:FFT))`
"""

# ╔═╡ Cell order:
# ╟─2ff7ee54-7bec-4e21-b804-e9fcc2f86a11
# ╠═49ca7fca-11dc-4f81-b74f-b54e8b252e10
# ╟─c8aaf8df-4f1b-4c07-84fd-8a5eb18ab0a7
# ╟─01b14950-89b8-4df1-981f-44a13f6026cb
# ╠═4b9838dc-bb6b-4bcc-8f51-e30f6f89be29
# ╠═4d8d2716-9c55-4f39-98f0-a08bbfb661eb
# ╠═4b069ffe-1b16-492d-acf3-83018ef16a43
# ╟─5353f914-a0eb-4e3b-bf5f-cbe90c63df88
