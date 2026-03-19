### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ bbd7a65f-132a-4092-bf17-5638bd4eb8ac
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido - Connectivity")
end

# ╔═╡ 99d7f8d4-05e5-47fd-a8ff-37b42d10b5af
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	using .PlutoTemplateBase
	using Serialization
end

# ╔═╡ c970cf66-ae42-403f-aec5-82f6dbb63069
notebook_intro("CONNECTIVITY")

# ╔═╡ c970cf66-ae42-403f-aec5-82f6dbb63068
md"""
Notebook por fase extraido de `Pluto/Notebook.jl`, centrado en:
- CSD,
- wPLI,
- y rutas de figuras/tablas para resultados finales.
"""

# ╔═╡ c2f60335-16cf-4902-ad20-f673d13f9fee
begin
	dir_csd = stage_dir(:CSD)
	dir_wpli = stage_dir(:wPLI)
	dir_fft = stage_dir(:FFT)

	(
		dir_csd = dir_csd,
		dir_wpli = dir_wpli,
		dir_fft = dir_fft,
		fig_csd = stage_dir(:CSD; kind = :figures),
		fig_wpli = stage_dir(:wPLI; kind = :figures),
		tab_wpli = stage_dir(:wPLI; kind = :tables),
	)
end

# ╔═╡ b0413220-3929-4266-b9d3-540ec57f191b
begin
	# Entrada tipica para conectividad: salida de FFT
	path_fft = joinpath(stage_dir(:FFT), "dict_EEG_FFT.bin")
	fft_data = if isfile(path_fft)
		Serialization.deserialize(path_fft)
	else
		@warn "No se encontro dict_EEG_FFT.bin; usando placeholder." path_fft
		Dict{String, Any}()
	end
	fft_data
end

# ╔═╡ f64af7c4-6a3c-4efc-a30a-6cbd9120be43
md"""
## Salidas esperadas (Connectivity)

- Datos CSD: `$(stage_dir(:CSD))`
- Datos wPLI: `$(stage_dir(:wPLI))`
- Figuras CSD: `$(stage_dir(:CSD; kind = :figures))`
- Figuras wPLI: `$(stage_dir(:wPLI; kind = :figures))`
- Tablas wPLI: `$(stage_dir(:wPLI; kind = :tables))`
"""

# ╔═╡ Cell order:
# ╟─bbd7a65f-132a-4092-bf17-5638bd4eb8ac
# ╠═99d7f8d4-05e5-47fd-a8ff-37b42d10b5af
# ╟─c970cf66-ae42-403f-aec5-82f6dbb63069
# ╟─c970cf66-ae42-403f-aec5-82f6dbb63068
# ╠═c2f60335-16cf-4902-ad20-f673d13f9fee
# ╠═b0413220-3929-4266-b9d3-540ec57f191b
# ╟─f64af7c4-6a3c-4efc-a30a-6cbd9120be43
