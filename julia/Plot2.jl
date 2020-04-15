
__precompile__()

module Plot2

export figure, RGB
export line_plot, scatter_plot, bar_plot, violin_plot, histogram_plot
export area_plot, matrix_plot, histogram_plot, genome_scatter_plot
export beeswarm_plot, survival_arrow_plot, stacked_bar_plot, box_plot
export rectangle, circle, errorbar
export xlim, ylim, xticks, yticks, subplot, ylabel, xlabel, title

using Helpers, Printf, Survival, DelimitedFiles
import PyCall

global PyPlot = nothing;
global drawing_figure = false;

function figure(block::Function, path::String; size=(5, 3))
	global drawing_figure, PyPlot;
	if drawing_figure == true
		block()
	else
		drawing_figure = true
		path = expanduser(path)
		if !isdir(dirname(abspath(path)))
			error("Cannot create figure $(path). Directory does not exist.")
		end

		# Import matplotlib.pyplot if we have not done so yet
		if PyPlot == nothing
			#matplotlib = PyCall.pyimport("matplotlib")
			#matplotlib.use("PDF")
			PyPlot = PyCall.pyimport("matplotlib.pyplot")
		end

		pyplot_fig = PyPlot.figure(figsize=size)
		try
			block()
			PyPlot.tight_layout()
			PyPlot.savefig(path)
			println(stderr, "Figure saved at $(path)")
		finally
			PyPlot.close()
			drawing_figure = false
		end
	end
	return nothing
end

figure(block::Function; kwargs...) = figure(block, expanduser("~/plot.pdf"); kwargs...)

xlabel(label::String) = PyPlot.xlabel(label)
ylabel(label::String) = PyPlot.ylabel(label)
title(text::String) = PyPlot.title(text)

function xticks(labels::AbstractVector{String}; rotation=90)
	PyPlot.xticks(1:length(labels), labels, rotation=rotation)
	PyPlot.tick_params(axis="x", length=0)
end

function yticks(labels::AbstractVector{String})
	PyPlot.yticks(1:length(labels), labels)
	PyPlot.tick_params(axis="y", length=0)
end

xticks(positions::AbstractVector) = PyPlot.xticks(positions)
yticks(positions::AbstractVector) = PyPlot.yticks(positions)

function xlim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.xlim(left=low); end
	if isfinite(high); PyPlot.xlim(right=high); end
	if log > 1; PyPlot.semilogx(basex=log); end
end

function ylim(low::Real, high::Real; log=NaN)
	if isfinite(low); PyPlot.ylim(bottom=low); end
	if isfinite(high); PyPlot.ylim(top=high); end
	if log > 1; PyPlot.semilogy(basey=log); end
end

function subplot(panel::Integer, rows::Integer)
	PyPlot.subplot(rows, 1, panel)
end

# function figure(path; size=(5, 3), panels=1, bottom_margin=NaN)
# 	pyplot_fig = PyPlot.figure(figsize=size)
# 	if 0 <= bottom_margin <= 1
# 		pyplot_fig[:subplots_adjust](bottom=bottom_margin)
# 	end
# 	global current_figure = Figure(path, size, 0, panels, NaN, NaN, pyplot_fig);
# end
# figure(; kwargs...) = figure("~/plot.pdf"; kwargs...)


struct RGB; r::UInt8; g::UInt8; b::UInt8; end
RGB(gray::Integer) = RGB(gray, gray, gray)
Base.broadcastable(r::RGB) = Ref(r)
hex(c::RGB) = @sprintf("#%02X%02X%02X", c.r, c.g, c.b)
hex(colors::Vector{RGB}) = hex.(colors)

function line_plot(x::AbstractVector, y::AbstractVector; color=RGB(0,0,0), line_width=1)
	figure() do
		PyPlot.plot(x, y, color=hex(color), linewidth=line_width)
		PyPlot.grid(true, axis="both", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function errorbar(args...; kwargs...)
	figure() do
		PyPlot.errorbar(args...; kwargs...)
	end
end

function genome_scatter_plot(chromosome::AbstractVector{String},
	position::AbstractVector, value::AbstractVector;
	chr_sizes="~/homo_sapiens/hg38.chrom.sizes", color=RGB(0))

	d = readdlm(expanduser(chr_sizes))
	chr_names = d[:, 1]
	chr_lens = Int32.(d[:, 2])
	chr_start = map(c -> sum(chr_lens[1:c-1]), 1:length(chr_names))
	gpos = chr_start[indexin(chromosome, chr_names)] .+ position

	chr_labels = [r"(19|21|M)$" in chr ? "" : replace(chr, "chr", "")
		for chr in chr_names]

	figure() do
		scatter_plot(gpos, value, color=color)
		PyPlot.grid(true, axis="both", which="major", linestyle="dashed")
		PyPlot.xlim(0, sum(chr_lens))
		PyPlot.tick_params(axis="x", which="both", length=0)
		PyPlot.gca().set_xticks(chr_start[2:end])
		PyPlot.gca().set_xticklabels([])
		PyPlot.gca().set_xticks(chr_start .+ chr_lens ./ 2, minor=true)
		PyPlot.gca().set_xticklabels(chr_labels, minor=true)
	end
end

function scatter_plot(x::AbstractVector, y::AbstractVector;
	size=2, color=RGB(0,0,0), border_color=nothing, marker="o", xerror=nothing, yerror=nothing, error_line_width=1, error_color=RGB(0,0,0), zorder=1)

	@assert(length(x) == length(y))

	# This is a workaround for the fact that Matplotlib misindexes the
	# "error_color" vector if some error ranges are NaN.
	if xerror != nothing
		@assert(xerror isa AbstractMatrix && Base.size(xerror, 1) == length(x))
		for k in 1:Base.size(xerror, 1)
			if any(isnan, xerror[k, :]); xerror[k, :] .= 0; end
		end
	end
	if yerror != nothing
		@assert(yerror isa AbstractMatrix && Base.size(yerror, 1) == length(x))
		for k in 1:Base.size(yerror, 1)
			if any(isnan, yerror[k, :]); yerror[k, :] .= 0; end
		end
	end

	figure() do
		if xerror isa AbstractMatrix
			PyPlot.hlines(y, xmin=xerror[:, 1], xmax=xerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end
		if yerror isa AbstractMatrix
			PyPlot.vlines(x, ymin=yerror[:, 1], ymax=yerror[:, 2], colors=hex(error_color), linewidth=error_line_width, zorder=0)
		end

		PyPlot.grid(true, axis="both", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		PyPlot.scatter(x, y, s=size, marker=marker, c=hex(color),
			edgecolors=(isa(border_color, RGB) ? hex(border_color) : "none"), zorder=zorder)
	end
end

function bar_plot(values::AbstractVector; box_width=0.8, color=RGB(0,0,0))
	figure() do
		PyPlot.bar(1:length(values), values, color=hex(color),
			width=box_width, linewidth=0,
			bottom=0.0000000000000001, zorder=10)
		PyPlot.xlim(0, length(values) + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
	end
end

function stacked_bar_plot(values::AbstractMatrix; box_width=0.8, colors=[])
	B = size(values, 1)    # Number of bars
	L = size(values, 2)    # Number of levels

	if isempty(colors)
		colors = [RGB(round(Int, (l - 1) / L)) for l in 1:L]
	end

	@assert(colors isa AbstractVector{RGB})
	if length(colors) != L
		error("Bar plot has $L levels but $(length(colors)) colors were provided in 'colors' keyword argument.")
	end

	figure() do
		bottoms = zeros(B)
		PyPlot.xlim(0, B + 1)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for level in 1:L
			PyPlot.bar(1:B, values[:, level], bottom=bottoms, color=hex(colors[level]), zorder=10)
			bottoms .+= values[:, level]
		end
	end
end

function box_plot(groups...; show_means=false)
	figure() do
		PyPlot.boxplot(groups, showmeans=show_means)
	end
end

function histogram_plot(bins, counts::Array; yrange=[],
	log_y=0, xlabel="", ylabel="", title="")
	# g = start_plot()
	# boxw = bins[2] - bins[1]
	# x = [bins[:] - boxw/2 bins[:] + boxw/2]'[:]
	# y = [counts[:] counts[:]]'[:]
	# data = write_data(g, [x[1] 0; x y; x[end] 0])
	# write(g, "set style fill solid 1 noborder\n")
	# write(g, "set xrange [$(bins[1]-boxw/2):$(bins[end]+boxw/2)]\n")
	# if log_y > 1
	# 	write(g, "set logscale y $(log_y)\n")
	# end
	# if !isempty(yrange)
	# 	write(g, "set yrange [$(yrange[1]):$(yrange[end])]\n")
	# elseif all(counts .>= 0)
	# 	write(g, "set yrange [0:*]\n")
	# end
	# write(g, "set xtics offset 0,0.5\n")
	# write(g, "set ytics offset 0.5,0\n")
	# if !isempty(title); write(g, "set title '$(title)' offset 0,-1.5\n"); end
	# if !isempty(xlabel); write(g, "set xlabel '$(xlabel)' offset 0,1\n"); end
	# if !isempty(ylabel); write(g, "set ylabel '$(ylabel)' offset 2,0\n"); end
	# write(g, "plot $data using 1:2 with filledcurves\n")
	# end_plot()
end

function area_plot(x::AbstractVector, y::AbstractVector; xbottom=nothing, ybottom=nothing, color=RGB(0))

	figure() do
		if xbottom == nothing && ybottom == nothing
			PyPlot.fill_between(x, y, facecolor=hex(color))
		end
		if xbottom isa Number || xbottom isa AbstractVector
			PyPlot.fill_betweenx(y, xbottom, x, facecolor=hex(color))
		end
		if ybottom isa Number || ybottom isa AbstractVector
			PyPlot.fill_between(x, ybottom, y, facecolor=hex(color))
		end
	end
end

#function histogram_plot(values::Array; path="~/plot.pdf")
	#centers =
	#_, h = hist(values, edges); centers = edges[1:end-1] + step(edges) / 2
	#histogram_plot(values, path=path, box_width=1.0)
#end

function violin_plot(bins::AbstractVector, density::AbstractArray; colors=[], halved=false)

	V = size(density, 2);   # How many violins?

	bins = Float64.(bins); density = Float64.(density)
	for v in 1:V
		density[:, v] /= maximum(density[:, v]) * 2.05
	end

	if isempty(colors)
		colors = map(v -> RGB(0, 0, 0), 1:V)
	end

	figure() do
		PyPlot.xlim(0.25, V + 0.75)
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for v in 1:V
			if halved
				x = vcat(v .- density[:, v], v, v, v .- density[1, v])
				y = vcat(bins, bins[end], bins[1], bins[1])
			else
				x = vcat(v .- density[:, v], v .+ density[end:-1:1, v])
				y = vcat(bins, bins[end:-1:1])
			end
			PyPlot.fill(x, y, hex(colors[v]), zorder=10)
		end
	end
end

function rectangle(x1::Real, x2::Real, y1::Real, y2::Real; edge_width=0, color=RGB(0, 0, 0), opacity=1.0)
	width = abs(x2 - x1); height = abs(y2 - y1);
	x = min(x1, x2); y = min(y1, y2);
	rect = PyPlot.matplotlib.patches.Rectangle((x, y), width, height, linewidth=edge_width, facecolor=hex(color), alpha=opacity)
	PyPlot.gca().add_patch(rect)
end

function ellipse(x::Real, y::Real, rx::Real, ry::Real; color=RGB(0, 0, 0))
	patch = PyPlot.matplotlib.patches.Ellipse(x, y, 2*rx, 2*ry,
		facecolor=hex(color))
	PyPlot.gca().add_patch(patch)
end

function beeswarm_plot(groups...; point_size=2, color=RGB(0,0,0))
	x = zeros(0); y = zeros(0);
	for (group, values) in enumerate(groups)
		x = vcat(x, group .+ clamp.(randn(length(values)) / 10, -0.4, 0.4))
		y = vcat(y, values)
	end
	figure() do
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)
		PyPlot.scatter(x, y, s=point_size, c=hex(color), marker=".", zorder=10)
	end
end

function survival_arrow_plot(survival::Vector{SurvivalTime}; color=RGB(0,0,0), start_time=[], head_length=NaN)

	S = length(survival)
	max_time = maximum(abs.(survival))
	if isnan(head_length); head_length = max_time / 20; end

	if isa(color, RGB); color = fill(color, S); end

	if isempty(start_time)
		start_time = zeros(S)
	end

	figure() do
		PyPlot.xlim(0.5, S + 0.5)
		PyPlot.ylim(0, max_time * 1.1)
		PyPlot.tick_params(axis="x", length=0)
		PyPlot.grid(true, axis="y", linestyle="dashed")
		PyPlot.gca().spines["top"].set_visible(false)
		PyPlot.gca().spines["right"].set_visible(false)

		for s in 1:S
			if survival[s] == event(0); continue; end
			PyPlot.arrow(s, start_time[s], 0, abs(survival[s]),
				color=hex(color[s]), width=0.3,
				head_width=(is_censored(survival[s]) ? 0.8 : 0),
				head_length=(is_censored(survival[s]) ? head_length : 0))
		end
	end
end


export MatrixLayer, matrix_plot
export GLYPH_NONE, GLYPH_TILE, GLYPH_SQUARE, GLYPH_SQUARE_BELOW, GLYPH_SQUARE_ABOVE, GLYPH_HORIZONTAL_BOX

const GLYPH_NONE = 0
const GLYPH_TILE = 1
const GLYPH_SQUARE = 2
const GLYPH_SQUARE_BELOW = 3
const GLYPH_SQUARE_ABOVE = 4
const GLYPH_HORIZONTAL_BOX = 5

struct MatrixLayer
	glyph::Array{Int8}
	color::Array{RGB}
end
MatrixLayer(rows::Integer, cols::Integer; glyph=GLYPH_NONE) =
	MatrixLayer(fill(glyph, rows, cols), fill(RGB(255, 255, 255), rows, cols))
Base.getindex(layer::MatrixLayer, rows, cols) = MatrixLayer(layer.glyph[rows, cols], layer.color[rows, cols])

function render_glyph(svg::IO, glyph::Int8, color::RGB,
	x::Int64, y::Int64, cell_width::Int64, cell_height::Int64)

	rgb = "rgb($(color.r), $(color.g), $(color.b))"
	if glyph == GLYPH_TILE
		write(svg, """<rect x="$(x)" y="$(y)" width="$(cell_width)" height="$(cell_height)" style="stroke: rgb(255, 255, 255); fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE
		write(svg, """<rect x="$(x + cell_width/4)" y="$(y + cell_height/4)" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE_BELOW
		write(svg, """<rect x="$(x + cell_width/6)" y="$(y + cell_height/6)" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_SQUARE_ABOVE
		write(svg, """<rect x="$(x + cell_width*2/6)" y="$(y + cell_height*2/6))" width="$(cell_width/2)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	elseif glyph == GLYPH_HORIZONTAL_BOX
		write(svg, """<rect x="$(x)" y="$(y + cell_height/4)" width="$(cell_width)" height="$(cell_height/2)" style="stroke: none; fill: $rgb" />""")
	end
end

function matrix_plot(layers::Array{MatrixLayer}; path="~/plot.svg", cell_width=50, cell_height=50)
	R, C = size(layers[1].glyph)
	canvas_w = C * cell_width
	canvas_h = R * cell_height
	svg = open(expanduser(path), "w")
	write(svg, """<svg width="$(canvas_w)" height="$(canvas_h)">""")
	for layer in layers
		@assert(size(layer.glyph) == (R, C))
		for r in 1:R
			for c in 1:C
				x = (c - 1) * cell_width
				y = (r - 1) * cell_height
				render_glyph(svg, layer.glyph[r, c], layer.color[r, c], x, y, cell_width, cell_height)
			end
		end
	end
	write(svg, "</svg>")
	close(svg)
end

matrix_plot(layer::MatrixLayer; kwargs...) = matrix_plot([layer]; kwargs...)

end
