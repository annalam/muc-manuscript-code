
"Robust estimator for function ax + b."
function theil_sen(x::Array, y::Array)
	assert(length(x) == length(y))
	pairs = rand(1:length(x), 1000, 2)
	a = fill(NaN, size(pairs, 1)); b = fill(NaN, size(pairs, 1));
	for p in 1:size(pairs, 1)
		dx = x[pairs[p, 2]] - x[pairs[p, 1]]
		if dx == 0; continue; end
		dy = y[pairs[p, 2]] - y[pairs[p, 1]]
		a[p] = dy / dx
		b[p] = y[pairs[p, 1]] - x[pairs[p, 1]] * a[p]
	end
	return (median(remove_nans(a)), median(remove_nans(b)))
end

function true_runs(values::AbstractArray{Bool})
	runs = Array{UnitRange{Int64}, 1}()
	start = 0
	for k in 1:length(values)
		if start == 0 && values[k] == true; start = k; end
		if start > 0 && values[k] == false
			push!(runs, start:(k-1)); start = 0
		end
	end
	if values[end] == true; push!(runs, start:length(values)); end
	return runs
end

function partition(keys::Array{ASCIIString})
	parts = Dict{ASCIIString, Array{Int}}()
	for k in 1:length(keys)
		push!(get!(parts, keys[k], []), k)
	end
	return parts
end

function runs(values::Array)
	runs = Array{UnitRange{Int64}, 1}()
	start = 1
	for k in 2:length(values)
		if values[k] != values[k-1]
			push!(runs, start:(k-1)); start = k
		end
	end
	push!(runs, start:length(values))
	return runs
end

function binned_mode(data, edges)
	counts = hist(logratio[chr_ranges[1, c]:chr_ranges[2, c]], edges)
	mode_idx = indmax(counts)
	return (edges[mode_idx] + edges[mode_idx+1]) / 2
end

#function zopen(path::AbstractString, mode)
#	path = expanduser(path)
#	if path == "-"; return STDIN; end
#	if ismatch(r"\.gz$", path)
#		return GZip.open(path, mode)
#	else
#		return open(path, mode)
#	end
#end
#zopen(path::AbstractString) = zopen(path, "r")

function complement(base::Char)
	if base == 'A'; return 'T'
	elseif base == 'C'; return 'G'
	elseif base == 'G'; return 'C'
	elseif base == 'T'; return 'A'
	else; return '*'; end
end
complement(dna::AbstractString) = map(c -> complement(c), dna)
revcomplement(dna::AbstractString) = reverse(complement(dna))

const human_chr = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
	"11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
	"21", "22", "X", "Y", "M"]
const human_chr2num = Dict(chr => c for (c, chr) in enumerate(human_chr))

chrom2num(chr::AbstractString) =
	get(human_chr2num, replace(chr, "chr", ""), 0)
chrom2num(chr::Integer) = chr
chrom2num(chrs::Array) = convert(Array{Int8}, map(chrom2num, chrs))

