
__precompile__()

module TSV

using DelimitedFiles, Dates

export read_tsv

# Helpers for reading tab-delimited files
mutable struct DelimitedFile
	headers::Vector{String}
	header_index::Dict{String, Int}
	data::Matrix{Any}
end
Base.size(f::DelimitedFile, args...) = Base.size(f.data, args...)
Base.lastindex(f::DelimitedFile, dim::Integer) = Base.size(f.data, dim)
Base.getindex(f::DelimitedFile, rows, col::AbstractString) = 
	Base.getindex(f.data, rows, get(f.header_index, col, 0))
Base.getindex(f::DelimitedFile, rows, col_regex::Regex) =
	Base.getindex(f.data, rows, col_regex .∈ f.headers)
Base.getindex(f::DelimitedFile, rows, cols) = Base.getindex(f.data, rows, cols)
Base.show(io::IO, f::DelimitedFile) = Base.show(io, f.data)  # FIXME

function read_tsv(tsv_file::IO; header=true, text=false)
	d = readdlm(tsv_file, '\t', text ? String : Any)
	headers = header ? d[1, :][:] : []
	header_index = Dict(zip(headers, 1:length(headers)))
	return DelimitedFile(headers, header_index, header ? d[2:end, :] : d)
end
read_tsv(cmd::Base.AbstractCmd; kwargs...) =
	open(f -> read_tsv(f; kwargs...), cmd)
read_tsv(tsv_path::AbstractString; kwargs...) =
	tsv_path == "-" ? read_tsv(STDIN; kwargs...) : open(f -> read_tsv(f; kwargs...), expanduser(tsv_path))

# TODO: Maybe put this here?
#function parse_date(x::AbstractString, century::Int)
#	if x == "" || x == "."; return nothing; end
#	if r"[0-9]+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-[0-9]+" in x
#		d = try Date(x, "d-u-y") catch; error("Invalid date: $(x)") end
#	elseif r"\d+/\d+/\d+" in x
#		d = try Date(x, "m/d/y") catch; error("Invalid date: $(x)") end
#	else
#		error("Invalid date: $(x)")
#	end
#	return d + Dates.Year(century)
#end

end