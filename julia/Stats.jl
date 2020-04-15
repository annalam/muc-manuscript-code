
__precompile__()

module Stats

using Helpers, Statistics, HypothesisTests, StatsBase

export fisher_exact, chisq_test, ranksum_test, proportion_confidence
export spearman, median, mean, quantile

function fisher_exact(a::AbstractVector{Bool}, b::AbstractVector{Bool})
	@assert(length(a) == length(b))
	println(stderr, [sum(a .& b) sum(a .& .!b); sum(.!a .& b) sum(.!a .& .!b)])
	# Here method=:minlike ensures that results are compatible with R
	return pvalue(FisherExactTest(sum(a .& b), sum(a .& .!b), sum(.!a .& b), sum(.!a .& .!b)), method=:minlike)
end

fisher_exact(a::Integer, b::Integer, c::Integer, d::Integer) =
	pvalue(FisherExactTest(Int64(a), Int64(b), Int64(c), Int64(d)))

function chisq_test(a::AbstractVector{Bool}, b::AbstractVector{Bool})
	@assert(length(a) == length(b))
	return pvalue(ChisqTest([sum(a .& b) sum(a .& .!b); sum(.!a .& b) sum(.!a .& .!b)]))
end

function proportion_confidence(x::Integer, n::Integer; alpha=0.05,
	method=:wilson)
	return confint(BinomialTest(x, n), alpha, method=method)
end

function ranksum_test(a::AbstractVector, b::AbstractVector; verbose=false)
	a_no_nans = filter(isfinite, a)
	b_no_nans = filter(isfinite, b)
	if verbose
		info("1st group: n = $(length(a_no_nans)), median = $(median(a_no_nans))")
		info("2nd group: n = $(length(b_no_nans)), median = $(median(b_no_nans))")
	end
	p = pvalue(MannWhitneyUTest(a_no_nans, b_no_nans))
	return p
end

spearman(a::AbstractVector, b::AbstractVector) = corspearman(a, b)

end