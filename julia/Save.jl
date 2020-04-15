
using Serialization

function save(filename::String, data)
	open(filename, "w") do f
		serialize(f, data)
	end
end

# Code is based on:
# https://github.com/JuliaIO/JLD2.jl/blob/master/src/loadsave.jl
macro save(filename, vars...)
	dict_push_exprs = Vector{Expr}(undef, length(vars))
	for i in 1:length(vars)
		dict_push_exprs[i] = :(data[$(string(vars[i]))] = $(esc(vars[i])))
	end

	quote
		data = Dict{String, Any}()
		$(Expr(:block, dict_push_exprs...))
		open($(esc(filename)), "w") do f
			serialize(f, data)
		end
	end
end

function load(filename::String)
	data = open(deserialize, filename, "r")
	for (key, value) in data
		eval(:($(Symbol(key)) = $(value)))
	end
end

macro load(filename)
	quote
		load($(esc(filename)))
	end
end

