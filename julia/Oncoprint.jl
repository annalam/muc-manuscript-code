
__precompile__()

module Oncoprint

using Plot2

export GeneStatus, oncoprint

mutable struct GeneStatus
	missing::Bool
	cn::Int8      # Between -2 and +2
	truncating::Int8
	missense::Int8
	rearrangements::Int8
	promoter_mutations::Int8
end
GeneStatus() = GeneStatus(false, 0, 0, 0, 0, 0)

cna_styles = [0 80 255; 150 200 255; 230 230 230; 255 150 150; 255 30 0]

function oncoprint(cells::Matrix{GeneStatus}; kwargs...)
	R, C = size(cells)
	tiles = MatrixLayer(R, C, glyph=GLYPH_TILE)
	layers = (MatrixLayer(R, C), MatrixLayer(R, C))
	for r in 1:R
		for c in 1:C
			gs = deepcopy(cells[r, c])   # Copy because we mutate it
			if gs.missing
				tiles.color[r, c] = RGB(255, 255, 255)
				continue
			end
			
			@assert(-2 <= gs.cn <= 2)
			tiles.color[r, c] = RGB(cna_styles[gs.cn + 3, :]...)

			hits = gs.truncating + gs.missense + gs.rearrangements + gs.promoter_mutations
			if hits == 0; continue; end

			for level in 1:min(hits, 2)
				glyph = hits == 1 ? GLYPH_SQUARE : (level == 1 ? GLYPH_SQUARE_BELOW : GLYPH_SQUARE_ABOVE)
				if gs.rearrangements > 0
					layers[level].glyph[r, c] = glyph
					layers[level].color[r, c] = RGB(0, 0, 0)
					gs.rearrangements -= 1
				elseif gs.truncating > 0
					layers[level].glyph[r, c] = glyph
					layers[level].color[r, c] = RGB(255, 200, 0)
					gs.truncating -= 1
				elseif gs.promoter_mutations > 0
					layers[level].glyph[r, c] = glyph
					layers[level].color[r, c] = RGB(162, 0, 172)
					gs.promoter_mutations -= 1
				elseif gs.missense > 0
					layers[level].glyph[r, c] = glyph
					layers[level].color[r, c] = RGB(120, 180, 60)
					gs.missense -= 1
				end
			end
		end
	end
	matrix_plot([tiles, layers...]; kwargs...)
end

end