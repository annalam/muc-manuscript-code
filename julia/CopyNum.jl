
__precompile__()

module CopyNum

using Printf, Statistics, DelimitedFiles
export CopyNumber, writetsv, median_decimate

mutable struct CopyNumber
	logratio::Array{Float32}
	median_af::Array{Float32}
	call::Array{Int16}
	hetz_snps::Array{Int32}    # Number of hetz SNPs used for median AF calc
	num_probes::Vector{Int32}   # Targeted seq amplicons or microarray probes
	sample::Vector{String}
	sample_noise::Vector{Float32}
	gene::Vector{String}
	chromosome::Vector{String}
	position::Array{Int32, 2}
end

Base.size(cn::CopyNumber, args...) = Base.size(cn.logratio, args...)
Base.getindex(cn::CopyNumber, rows, cols) = CopyNumber(
	cn.logratio[rows, cols], cn.median_af[rows, cols], cn.call[rows, cols],
	cn.hetz_snps[rows, cols], cn.num_probes[rows],
	cn.sample[cols], cn.sample_noise[cols], 
	cn.gene[rows], cn.chromosome[rows], cn.position[rows, :])

#function Base.getindex(cn::CopyNumber, gene_name::AbstractString, cols)
#	row = get(cn.gene_name_to_row, gene_name, 0)
#	if row == 0; error("Gene '$(gene_name)' not found."); end
#	return cn[row, cols]
#end

function CopyNumber(genes::AbstractVector, samples::AbstractVector)
	S = length(samples); G = length(genes)
	#gene_name_to_row = Dict(g => row for (row, g) in enumerate(genes))
	return CopyNumber(zeros(Float32, G, S), fill(Float32(NaN), G, S),
		zeros(Int16, G, S), zeros(Int32, G, S), zeros(Int32, G),
		samples, fill(Float32(NaN), S), genes, fill("", G), zeros(Int32, G, 2))
end

function CopyNumber(tsv_path::AbstractString)
	d = readdlm(tsv_path, '\t'); headers = d[1, :][:];
	genes = d[2:end, 1][:]; samples = headers[5:end];
	cn = CopyNumber(genes, samples)
	for g in 1:length(genes)
		cn.chromosome[g] = d[1+g, 2]
		cn.position[g, 1] = d[1+g, 3]
		cn.position[g, 2] = d[1+g, 4]
		for s in 1:length(samples)
			fields = split(d[1+g, 4+s], ':')
			cn.call[g, s] = parse(Int, fields[1])
			cn.logratio[g, s] = parse(Float32, fields[2])
			cn.median_af[g, s] = parse(Float32, fields[3])
			cn.hetz_snps[g, s] = parse(Int, fields[4])
		end
	end
	return cn
end

function writetsv(out::IO, cn::CopyNumber)
	write(out, "GENE\tCHROMOSOME\tSTART\tEND")
	for s in cn.sample; write(out, "\t$(s)"); end; write(out, '\n')

	for g in 1:length(cn.gene)
		@printf(out, "%s\t%s\t%d\t%d", cn.gene[g], cn.chromosome[g],
			cn.position[g, 1], cn.position[g, 2])
		for s in 1:length(cn.sample)
			@printf(out, "\t%d:%.2f:%.3f:%d", cn.call[g, s],
				cn.logratio[g, s], cn.median_af[g, s], cn.hetz_snps[g, s])
		end
		write(out, '\n')
	end
end
writetsv(out_path::AbstractString, cn::CopyNumber; kwargs...) =
	open(fd -> writetsv(fd, cn; kwargs...), expanduser(out_path), "w")

function median_decimate(values::AbstractVector, fold::Integer)
	if fold == 1; return copy(values); end
	if isempty(values); return zeros(eltype(values), 0); end
	starts = 1:fold:length(values)
	decimated = zeros(eltype(values), length(starts))
	for k in 1:length(starts)-1     # Last window handled as special case
		decimated[k] = median(values[(0:fold-1) .+ starts[k]])
	end
	decimated[length(starts)] = median(values[starts[end]:end])
	return decimated
end

end