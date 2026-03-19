### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 4d8ca8b0-149f-4f6f-93f4-4e2d7f1a0b01
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido - Preprocessing")
end

# ╔═╡ d249d3f0-cf95-44e1-aec6-8ca4ad0fe102
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	using .PlutoTemplateBase
	using DataFrames, CSV, Serialization
end

# ╔═╡ e857112f-df27-4f7f-bcc7-34421e0c3103
notebook_intro("PREPROCESSING")

# ╔═╡ 204da7a3-9f60-4918-97f4-06a2d621dc15
md"""
Notebook por fase extraido de `Pluto/Notebook.jl`, centrado en:
- carga de datos (`raw` y `electrodes`),
- construccion de `dict_EEG.bin` (IO),
- y rutas de salida para `filtering`.
"""

# ╔═╡ b1a100f5-2f4f-4d30-b4c4-1dca20c50e04
begin
	# Rutas sin hardcodeo (segun src/modules/paths.jl)
	dir_raw = raw_dir()
	dir_electrodes = electrodes_dir()
	dir_io = stage_dir(:IO)
	dir_filtering = stage_dir(:filtering)

	(
		dir_raw = dir_raw,
		dir_electrodes = dir_electrodes,
		dir_io = dir_io,
		dir_filtering = dir_filtering,
	)
end

# ╔═╡ c5a95b03-b0c6-4bf8-a5ab-c5a31af34205
begin
	raw_path = joinpath(raw_dir(), "sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv")
	data_raw = if isfile(raw_path)
		CSV.read(raw_path, DataFrame; delim = '\t')
	else
		@warn "No se encontro raw TSV; usando DataFrame vacio (modo publico/export)." raw_path
		DataFrame(Channel = String[])
	end
	data_raw
end

# ╔═╡ 11b0a622-c1e1-4e05-98c3-cf3a2f72a306
begin
	electrodes_path = joinpath(electrodes_dir(), "sub-M05_ses-T2_electrodes.tsv")
	electrodes = if isfile(electrodes_path)
		CSV.read(electrodes_path, DataFrame; delim = '\t')
	else
		@warn "No se encontro electrodes.tsv; usando tabla vacia (modo publico/export)." electrodes_path
		DataFrame(name = String[], x = Float64[], y = Float64[], z = Float64[], type = String[])
	end
	electrodes
end

# ╔═╡ 3249aa95-e924-4e96-bf72-30c9a02c8607
begin
	dict_EEG = if nrow(data_raw) > 0 && (:Channel in names(data_raw)) && (size(data_raw, 2) > 1)
		mat = Matrix(data_raw[:, Not(:Channel)])
		Dict(data_raw.Channel[i] => vec(mat[i, :]) for i in 1:size(data_raw, 1))
	else
		Dict{String, Vector{Float64}}()
	end
	dict_EEG
end

# ╔═╡ d359220a-c31c-4ce4-91e6-b58a7a9e7d08
begin
	isdir(dir_io) || mkpath(dir_io)
	path_dict = joinpath(dir_io, "dict_EEG.bin")
	Serialization.serialize(path_dict, dict_EEG)
	"dict_EEG guardado en: $(abspath(path_dict))"
end

# ╔═╡ 8de2f84f-2f33-44b7-b9cb-4a88fd6f8f09
md"""
## Siguiente etapa

Los resultados de filtrado deben guardarse en:
`$(stage_dir(:filtering))`
"""

# ╔═╡ Cell order:
# ╟─4d8ca8b0-149f-4f6f-93f4-4e2d7f1a0b01
# ╠═d249d3f0-cf95-44e1-aec6-8ca4ad0fe102
# ╟─e857112f-df27-4f7f-bcc7-34421e0c3103
# ╟─204da7a3-9f60-4918-97f4-06a2d621dc15
# ╠═b1a100f5-2f4f-4d30-b4c4-1dca20c50e04
# ╠═c5a95b03-b0c6-4bf8-a5ab-c5a31af34205
# ╠═11b0a622-c1e1-4e05-98c3-cf3a2f72a306
# ╠═3249aa95-e924-4e96-bf72-30c9a02c8607
# ╠═d359220a-c31c-4ce4-91e6-b58a7a9e7d08
# ╟─8de2f84f-2f33-44b7-b9cb-4a88fd6f8f09
