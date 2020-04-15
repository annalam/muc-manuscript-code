
__precompile__()

module QualityMetrics

export plot_quality_metrics

using Plot2, Printf, DelimitedFiles, Statistics, Helpers

function plot_quality_metrics(samples; labels=[], output="", pdf_filepath="~/plot.pdf")
	S = length(samples)
	if isempty(labels)
		labels = copy(samples)
	else
		if !isa(labels, AbstractVector{String})
			error("Labels must be given as a vector of strings.")
		end
		if length(labels) != length(samples)
			error("Number of labels and files must match.")
		end
	end

	total = zeros(S); aligned = zeros(S); duplicate = zeros(S)
	on_target = zeros(S)
	coverage = zeros(countlines(open("$(samples[1]).coverage")), S)
	coverage_bins = readdlm("$(samples[1]).coverage", '\t')[:, 1]

	for (s, sample) in enumerate(samples)
		d = Dict{AbstractString, Int64}()
		for line in eachline(open("$(sample).general"))
			if ": " in line == false; continue; end
			key, value = split(line, ": ")
			d[key] = parse(Int, split(value, ' ')[1])
		end
		total[s] = d["Total reads"]
		duplicate[s] = d["Duplicate reads"]
		aligned[s] = d["Aligned reads"]
		on_target[s] = get(d, "On-target reads", 0)
		coverage[:, s] = readdlm("$(sample).coverage", '\t')[:, 2]
	end

	# Calculate median coverage of target regions for each sample
	sample_median_coverage = zeros(S)
	sample_10_percentile_coverage = zeros(S)
	for s in 1:S
		if sum(coverage[:, s]) == 0
			warn("Sample $(samples[s]) has no coverage.")
			continue
		end
		cum_fraction = cumsum(coverage[:, s]) ./ sum(coverage[:, s])
		sample_median_coverage[s] = findfirst(cum_fraction .>= 0.5)
		sample_10_percentile_coverage[s] = findfirst(cum_fraction .>= 0.1)
	end
	quartiles = quantile(sample_median_coverage, [0.25, 0.5, 0.75])
	println(stderr, "Median sample coverage: $(quartiles[2]) (IQR $(quartiles[1]) - $(quartiles[3]))")
	
	# Calculate y-axis range for coverage violin plot
	ymax = [x for x in vcat(100:100:1000, 1500:500:5000, 6000:1000:20000)
		if all(sample_median_coverage * 1.5 .< x)][1]

	max_label_chars = maximum(length.(labels))

	figure(pdf_filepath, size=(1 + 0.15 * S, 6 + 0.25 * max_label_chars)) do
		subplot(1, 3)
		stacked_bar_plot(
			hcat(aligned - duplicate, duplicate, total - aligned) ./ 1_000_000,
			colors=[RGB(143, 178, 88), RGB(254, 216, 51), RGB(246, 126, 32)])
		xticks(labels); ylabel("Million reads")

		subplot(2, 3)
		bar_plot(on_target ./ aligned * 100)
		ylim(0, 100); xticks(labels); ylabel("On-target (%)")

		subplot(3, 3)
		violin_plot(coverage_bins[10:10:end], coverage[10:10:end, :])
		ylim(0, ymax); xticks(labels); ylabel("Target region coverage")

		#violin_plot(1:10:5000, fragment_lens[1:10:end, :], labels=samples, ylabel="Fragment length", yrange=0:300)
	end

	# Export quality metrics as a table
	out = output == "" ? stdin : open(expanduser(output), "w")
	println(out, "SAMPLE\tTOTAL READS (MILLION)\tCOVERAGE (MEDIAN)\tCOVERAGE (10%)\tON-TARGET %\tDUPLICATE%")
	for s in 1:S
		@printf(out, "%s\t%.1f\t%d\t%d\t%.1f\t%.1f\n",
			labels[s], total[s] / 1_000_000,
			sample_median_coverage[s], sample_10_percentile_coverage[s],
			on_target[s] ./ aligned[s] * 100, duplicate[s] ./ aligned[s] * 100)
	end
end

end
