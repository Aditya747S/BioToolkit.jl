module BioToolkitPlotsExt

using BioToolkit
using Plots

import Plots: plot

import BioToolkit: GenomeBrowser, GeneRenderPlan, CoverageRenderPlan, AlignmentRenderPlan, GeneTrack, CoverageTrack, AlignmentTrack, PhyloTree
import BioToolkit: render_browser, coordinates, get_terminals
import BioToolkit.BioPlotting: VolcanoPlotResult, MAPlotResult, ClusteredHeatmapResult, ManhattanPlotResult, QQPlotResult, ForestPlotResult
import BioToolkit.BioPlotting: volcano_data, ma_data, clustered_heatmap, manhattan_data, qq_data
import BioToolkit.Clinical: KaplanMeierResult, CoxResult, ROCResult, CIFResult, OncoprintResult, forest_plot, oncoprint_plot
import BioToolkit.GWAS: GWASResult, MetaAnalysisResult

const BioPlotting = BioToolkit.BioPlotting

function _save_plot(path::Union{Nothing,AbstractString}, figure)
    path === nothing && return nothing
    savefig(figure, path)
    return path
end

function _km_step_data(result::KaplanMeierResult)
    x = Float64[0.0]
    y = Float64[1.0]
    current = 1.0
    for (time, survival) in zip(result.time, result.survival)
        push!(x, time)
        push!(y, current)
        push!(x, time)
        push!(y, survival)
        current = survival
    end
    return x, y
end

function _browser_track_plot(plan::GeneRenderPlan, track; title::AbstractString="")
    figure = Plots.plot(; title=title, xlabel="genomic position", ylabel="gene row", legend=false, background_color=:white)
    if plan.lod == :chromosome
        x = collect(plan.bin_edges)
        isempty(x) || plot!(figure, x, plan.density; color=track.color, linewidth=2, label=false)
        isempty(x) || scatter!(figure, x, plan.density; color=track.color, markersize=3, label=false)
    else
        for placement in plan.placements
            for segment in placement.segments
                plot!(figure, [segment.left, segment.right], fill(placement.row, 2); color=track.color, linewidth=2, label=false)
            end
            placement.label === nothing || annotate!(figure, placement.interval.left, placement.row + 0.2, text(placement.label, 8, :black, :left))
        end
        yticks!(figure, (1:max(plan.row_count, 1), string.(1:max(plan.row_count, 1))))
    end
    return figure
end

function _browser_track_plot(plan::CoverageRenderPlan, track; title::AbstractString="")
    figure = Plots.plot(; title=title, xlabel="genomic position", ylabel="coverage", legend=false, background_color=:white)
    if !isempty(plan.bins)
        x = Float64[]
        y = Float64[]
        for bin in plan.bins
            push!(x, bin.left)
            push!(x, bin.right)
            push!(y, bin.mean_value)
            push!(y, bin.mean_value)
        end
        plot!(figure, x, y; color=track.color, linewidth=2, label=false)
        scatter!(figure, [bin.left for bin in plan.bins], [bin.mean_value for bin in plan.bins]; color=track.color, markersize=3, label=false)
    end
    return figure
end

function _browser_track_plot(plan::AlignmentRenderPlan, track; title::AbstractString="")
    figure = Plots.plot(; title=title, xlabel="genomic position", ylabel="read row", legend=false, background_color=:white)
    for read in plan.reads
        for (left, right) in read.blocks
            plot!(figure, [left, right], fill(read.row, 2); color=track.color, linewidth=3, label=false)
        end
        for mismatch in read.mismatches
            scatter!(figure, [mismatch.position], [read.row]; color=:red, markersize=4, label=false)
        end
        annotate!(figure, read.left, read.row + 0.2, text(read.qname, 8, :black, :left))
    end
    yticks!(figure, (1:max(length(plan.reads), 1), string.(1:max(length(plan.reads), 1))))
    return figure
end

function volcano_plot(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, label_top::Integer=10, kwargs...)
    points, summary = volcano_data(results; lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff, label_top=label_top)
    figure = Plots.plot(; title="Volcano plot", xlabel="log2 fold change", ylabel="-log10 p-value", size=(980, 700), legend=:outerright, background_color=:white, left_margin=15Plots.mm, bottom_margin=15Plots.mm, top_margin=10Plots.mm, right_margin=28Plots.mm)
    palette = Dict(:up => :firebrick, :down => :royalblue, :neutral => :gray)
    labels = Dict(:up => "Upregulated", :down => "Downregulated", :neutral => "Neutral")
    for category in (:neutral, :up, :down)
        selected = [point for point in points if point.category == category]
        isempty(selected) && continue
        scatter!(figure, [point.log2_fold_change for point in selected], [point.negative_log10_pvalue for point in selected]; color=palette[category], markersize=7, markerstrokewidth=0, label=labels[category])
    end
    vline!(figure, [-Float64(lfc_cutoff), Float64(lfc_cutoff)]; color=:gray, linestyle=:dash, label=false)
    hline!(figure, [-log10(clamp(Float64(fdr_cutoff), eps(Float64), 1.0))]; color=:gray, linestyle=:dash, label=false)
    labeled = sort([point for point in points if point.labeled]; by = point -> point.negative_log10_pvalue, rev=true)
    for (index, point) in enumerate(Iterators.take(labeled, 8))
        offset = 0.25 + 0.15 * ((index - 1) % 3)
        annotate!(figure, point.log2_fold_change, point.negative_log10_pvalue + offset, text(point.gene_id, 7, :black, :center, :bottom; rotation=90))
    end
    return BioPlotting.VolcanoPlotResult(points, summary, figure)
end

function ma_plot(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, kwargs...)
    points, summary = ma_data(results; lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff)
    figure = Plots.plot(; title="MA plot", xlabel="log2 abundance", ylabel="log2 fold change", size=(900, 660), legend=:topright, background_color=:white, left_margin=15Plots.mm, bottom_margin=15Plots.mm, top_margin=10Plots.mm)
    palette = Dict(:significant => :darkmagenta, :background => :slategray)
    labels = Dict(:significant => "Significant", :background => "Background")
    for category in (:background, :significant)
        selected = [point for point in points if point.category == category]
        isempty(selected) && continue
        scatter!(figure, [point.abundance for point in selected], [point.log2_fold_change for point in selected]; color=palette[category], markersize=7, markerstrokewidth=0, label=labels[category])
    end
    hline!(figure, [0.0]; color=:gray, linestyle=:dash, label=false)
    return BioPlotting.MAPlotResult(points, summary, figure)
end

function manhattan_plot(result::GWASResult; pvalue_threshold::Real=5e-8, label_top::Integer=10, kwargs...)
    points, summary = manhattan_data(result; pvalue_threshold=pvalue_threshold, label_top=label_top)
    order = unique(result.chromosomes)
    sort!(order; by = chromosome -> begin
        stripped = replace(lowercase(String(chromosome)), "chr" => "")
        parsed = tryparse(Int, stripped)
        parsed === nothing ? typemax(Int) : parsed
    end)
    _, _, centers, total_span = BioPlotting._chromosome_layout(result)
    figure = Plots.plot(; size=(1200, 620), background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=6Plots.mm, right_margin=4Plots.mm, top_margin=5Plots.mm, bottom_margin=6Plots.mm)
    xs = [point.cumulative_position for point in points]
    ys = [point.negative_log10_pvalue for point in points]
    alternating = Dict(order[index] => (isodd(index) ? "#4f6d7a" : "#9fb3c8") for index in eachindex(order))
    colors = [point.category == :significant ? "#d1495b" : alternating[point.chromosome] for point in points]
    scatter!(figure, xs, ys; color=colors, markersize=4.8, markerstrokewidth=0, xlabel="genomic position", ylabel="-log10(p)", title="Manhattan plot")
    xticks!(figure, ([get(centers, chromosome, 0.0) for chromosome in order], order))
    xlims!(figure, (0.0, max(total_span, 1.0)))
    ypeak = isempty(ys) ? 0.0 : maximum(ys)
    ylims!(figure, (0.0, max(ypeak * 1.08, -log10(clamp(Float64(pvalue_threshold), eps(Float64), 1.0)) * 1.15)))
    if any(point.labeled for point in points)
        labeled = [point for point in points if point.labeled]
        annotate!(figure, [(point.cumulative_position, point.negative_log10_pvalue, text(point.snp_id, 7, :black, :left)) for point in labeled])
    end
    hline!(figure, [-log10(clamp(Float64(pvalue_threshold), eps(Float64), 1.0))]; linestyle=:dash, color="#7d1f2f", linewidth=1.5)
    return BioPlotting.ManhattanPlotResult(points, summary, figure)
end

function qq_plot(result::GWASResult; kwargs...)
    points, summary = qq_data(result)
    figure = Plots.plot(; size=(840, 760), background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=6Plots.mm, right_margin=4Plots.mm, top_margin=5Plots.mm, bottom_margin=6Plots.mm, aspect_ratio=:equal)
    xs = [point.expected for point in points]
    ys = [point.observed for point in points]
    if !isempty(points)
        n = length(points)
        confidence_low = [max(-log10(quantile(Beta(index, n - index + 1), 0.025)), 0.0) for index in 1:n]
        confidence_high = [max(-log10(quantile(Beta(index, n - index + 1), 0.975)), 0.0) for index in 1:n]
        plot!(figure, xs, confidence_high; fillrange=confidence_low, fillalpha=0.18, linecolor=:transparent, fillcolor="#9db4c0")
    end
    scatter!(figure, xs, ys; color="#173f5f", markersize=4.2, markerstrokewidth=0, xlabel="expected -log10(p)", ylabel="observed -log10(p)", title="QQ plot")
    max_value = isempty(points) ? 1.0 : max(maximum(xs), maximum(ys))
    plot!(figure, [0.0, max_value], [0.0, max_value]; color="#7a7a7a", linestyle=:dash, linewidth=1.4)
    return BioPlotting.QQPlotResult(points, summary, figure)
end

function gwas_forest_plot(result::MetaAnalysisResult; kwargs...)
    points = BioPlotting.ForestPoint[]
    for index in eachindex(result.snp_ids)
        beta = result.beta[index]
        stderr = result.standard_error[index]
        lower = beta - 1.96 * stderr
        upper = beta + 1.96 * stderr
        push!(points, BioPlotting.ForestPoint(result.snp_ids[index], beta, stderr, lower, upper, result.pvalue[index], result.qvalue[index]))
    end
    order = sortperm([point.qvalue for point in points])
    selected = points[order]
    figure = Plots.plot(; size=(920, max(320, 42 * max(length(selected), 1))), background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=7Plots.mm, right_margin=5Plots.mm, top_margin=6Plots.mm, bottom_margin=6Plots.mm)
    if !isempty(selected)
        y_positions = collect(length(selected):-1:1)
        scatter!(figure, [point.beta for point in selected], y_positions; xerror=[1.96 * point.standard_error for point in selected], color="#173f5f", markersize=5.5, markerstrokewidth=0)
        for (index, point) in enumerate(selected)
            plot!(figure, [point.lower, point.upper], [y_positions[index], y_positions[index]]; color="#173f5f", linewidth=2.2)
        end
        yticks!(figure, (y_positions, [point.label for point in selected]))
        vline!(figure, [0.0]; color="#7d1f2f", linestyle=:dash, linewidth=1.3)
        xlims!(figure, (minimum([point.lower for point in selected]) - 0.1, maximum([point.upper for point in selected]) + 0.1))
        ylims!(figure, (0.5, length(selected) + 0.5))
    end
    summary = (
        total = length(selected),
        study_count = result.study_count,
        random_effects = result.method == "random_effects",
    )
    return BioPlotting.ForestPlotResult(selected, summary, figure)
end

function _roc_figure(result::ROCResult; title::AbstractString="Survival ROC", size=(820, 660), kwargs...)
    order = sortperm(result.fpr)
    fpr = vcat(0.0, result.fpr[order], 1.0)
    tpr = vcat(0.0, result.tpr[order], 1.0)
    figure = Plots.plot(fpr, tpr; seriestype=:steppost, linewidth=2.8, color=:navy, label="ROC curve", xlabel="False positive rate", ylabel="True positive rate", title="$(title) (AUC=$(round(result.auc; digits=3)))", legend=:bottomright, size=size, left_margin=15Plots.mm, bottom_margin=15Plots.mm, top_margin=10Plots.mm)
    plot!(figure, [0.0, 1.0], [0.0, 1.0]; linestyle=:dash, color=:gray, linewidth=1.4, label="No-skill line")
    return figure
end

function _cif_figure(result::CIFResult; title::AbstractString="Competing risks cumulative incidence", size=(820, 620), kwargs...)
    figure = Plots.plot(; xlabel="Time", ylabel="Cumulative incidence", title=title, legend=:bottomright, size=size)
    for (cause, values) in sort(collect(result.cumulative_incidence); by = first)
        plot!(figure, result.time, values; linewidth=2.4, label="Cause $(cause)")
    end
    if !isempty(result.censoring)
        plot!(figure, result.time, result.censoring; linewidth=1.8, linestyle=:dash, color=:black, label="Censoring")
    end
    return figure
end

function Plots.plot(browser::BioToolkit.GenomeBrowser; kwargs...)
    plans = render_browser(browser)
    figures = [_browser_track_plot(item.plan, item.track; title=index == 1 ? browser.title : "") for (index, item) in enumerate(plans)]
    width = browser.viewport.pixel_width
    track_heights = sum(hasproperty(item.track, :height) ? getproperty(item.track, :height) : 100 for item in plans)
    height = max(300, track_heights + browser.gap * max(length(plans) - 1, 0))
    return Plots.plot(figures...; layout=(length(figures), 1), size=(width, height))
end

function Plots.plot(plan::BioToolkit.GeneRenderPlan; track::Union{Nothing,BioToolkit.GeneTrack}=nothing, kwargs...)
    plot_track = track === nothing ? BioToolkit.GeneTrack(BioToolkit.GenomicInterval[]; height=120) : track
    return _browser_track_plot(plan, plot_track)
end

function Plots.plot(plan::BioToolkit.CoverageRenderPlan; track::Union{Nothing,BioToolkit.CoverageTrack}=nothing, kwargs...)
    plot_track = track === nothing ? BioToolkit.CoverageTrack(nothing; height=150) : track
    return _browser_track_plot(plan, plot_track)
end

function Plots.plot(plan::BioToolkit.AlignmentRenderPlan; track::Union{Nothing,BioToolkit.AlignmentTrack}=nothing, kwargs...)
    plot_track = track === nothing ? BioToolkit.AlignmentTrack(nothing; height=300) : track
    return _browser_track_plot(plan, plot_track)
end

function Plots.plot(tree::BioToolkit.PhyloTree; title::AbstractString="Phylogenetic tree", kwargs...)
    positions = coordinates(tree)
    leaves = [node.name for node in get_terminals(tree)]
    figure = Plots.plot(; xlabel="branch length", ylabel="leaf order", title=title, legend=false, size=(900, 600), yflip=true, background_color=:white)

    function _draw(node::PhyloTree)
        x1, y1 = positions[node]
        for child in node.children
            x2, y2 = positions[child]
            plot!(figure, [x1, x2], [y1, y1]; color=:black, linewidth=1.6, label=false)
            plot!(figure, [x2, x2], [y1, y2]; color=:black, linewidth=1.6, label=false)
            _draw(child)
        end
        isempty(node.children) || return
        annotate!(figure, x1 + 0.15, y1, text(isempty(node.name) ? "leaf" : node.name, 8, :black, :left))
    end

    _draw(tree)
    yticks!(figure, (1:length(leaves), leaves))
    max_x = maximum(position[1] for position in values(positions))
    xlims!(figure, (0.0, max_x + 1.5))
    return figure
end

function pair_plot(data::AbstractMatrix{<:Real}; labels::Vector{String}, groups::Vector{String}, kwargs...)
    nvars = size(data, 2)
    unique_groups = unique(groups)
    palette = [:navy, :firebrick, :seagreen, :darkorange, :purple]
    colors = Dict(group => palette[mod1(index, length(palette))] for (index, group) in enumerate(unique_groups))
    figure = Plots.plot(layout=(nvars, nvars), size=(950, 950), legend=false, background_color=:white)
    for row in 1:nvars, col in 1:nvars
        subplot = (row - 1) * nvars + col
        if row == col
            for group in unique_groups
                selected = groups .== group
                histogram!(figure[subplot], data[selected, col]; bins=18, color=colors[group], alpha=0.55, xlabel=labels[col], ylabel="count", label=false)
            end
        else
            for group in unique_groups
                selected = groups .== group
                scatter!(figure[subplot], data[selected, col], data[selected, row]; color=colors[group], markersize=2.5, markerstrokewidth=0, xlabel=labels[col], ylabel=labels[row], label=false)
            end
        end
    end
    return figure
end

function violin_plot(values::AbstractVector{<:Real}, groups::AbstractVector{<:AbstractString}; kwargs...)
    figure = Plots.plot(; title="Violin plot", xlabel="Group", ylabel="Value", legend=false, size=(860, 620), background_color=:white)
    unique_groups = unique(groups)
    for (index, group) in enumerate(unique_groups)
        selected = values[groups .== group]
        isempty(selected) && continue
        edges = range(minimum(selected), maximum(selected), length=28)
        counts = zeros(Float64, length(edges) - 1)
        for value in selected
            bin = clamp(searchsortedlast(edges, value), 1, length(counts))
            counts[bin] += 1
        end
        scaled = maximum(counts) > 0 ? counts ./ maximum(counts) .* 0.35 : counts
        centers = [0.5 * (edges[i] + edges[i + 1]) for i in 1:length(counts)]
        left = vcat(fill(index, length(centers)) .- scaled, reverse(fill(index, length(centers)) .+ scaled))
        right = vcat(centers, reverse(centers))
        plot!(figure, Shape(left, right); color=:cornflowerblue, alpha=0.65, linecolor=:cornflowerblue, label=false)
        scatter!(figure, fill(index, length(selected)), selected; color=:black, markersize=2.5, markerstrokewidth=0, alpha=0.35, label=false)
    end
    xticks!(figure, (1:length(unique_groups), unique_groups))
    return figure
end

function mosaic_plot(table::AbstractMatrix{<:Real}; row_labels::Vector{String}, col_labels::Vector{String}, kwargs...)
    figure = Plots.plot(; title="Mosaic plot", xaxis=false, yaxis=false, aspect_ratio=:equal, legend=false, size=(900, 600), background_color=:white)
    row_totals = sum(table, dims=2)[:, 1]
    total = sum(row_totals)
    x0 = 0.0
    colors = [:steelblue, :darkorange, :seagreen, :firebrick, :orchid]
    for (row_index, row_label) in enumerate(row_labels)
        width = total == 0 ? 0.0 : row_totals[row_index] / total
        y0 = 0.0
        row_sum = sum(table[row_index, :])
        for (col_index, col_label) in enumerate(col_labels)
            height = row_sum == 0 ? 0.0 : table[row_index, col_index] / row_sum
            if height > 0
                xs = [x0, x0 + width, x0 + width, x0]
                ys = [y0, y0, y0 + height, y0 + height]
                plot!(figure, Shape(xs, ys); color=colors[mod1(col_index, length(colors))], linecolor=:white, alpha=0.8, label=false)
                annotate!(figure, x0 + width / 2, y0 + height / 2, text(string(row_label, " / ", col_label), 7, :white, :center))
            end
            y0 += height
        end
        annotate!(figure, x0 + width / 2, 1.03, text(row_label, 8, :black, :center))
        x0 += width
    end
    return figure
end

function circos_plot(links::Vector{Tuple{String,String,Float64}}; kwargs...)
    nodes = unique(vcat([link[1] for link in links], [link[2] for link in links]))
    n = max(length(nodes), 1)
    angles = Dict(node => 2π * (index - 1) / n for (index, node) in enumerate(nodes))
    xs = [cos(angles[node]) for node in nodes]
    ys = [sin(angles[node]) for node in nodes]
    figure = Plots.plot(; aspect_ratio=:equal, xlabel="", ylabel="", title="Circos-style chord plot", legend=:outerright, size=(900, 900), xaxis=false, yaxis=false, background_color=:white)
    plot!(figure, vcat(xs, xs[1]), vcat(ys, ys[1]); color=:gray, linewidth=2, label=false)
    for (index, node) in enumerate(nodes)
        annotate!(figure, 1.08 * xs[index], 1.08 * ys[index], text(node, 8, :black, :center))
    end
    for (left, right, weight) in links
        t = range(0, 1, length=60)
        x = (1 .- t).^2 .* cos(angles[left]) .+ t.^2 .* cos(angles[right])
        y = (1 .- t).^2 .* sin(angles[left]) .+ t.^2 .* sin(angles[right])
        plot!(figure, x, y; color=weight > 0.6 ? :firebrick : :steelblue, alpha=clamp(weight, 0.2, 0.9), linewidth=2, label=false)
    end
    return figure
end

Plots.plot(result::VolcanoPlotResult; kwargs...) = result.figure
Plots.plot(result::MAPlotResult; kwargs...) = result.figure
Plots.plot(result::ClusteredHeatmapResult; kwargs...) = result.figure
Plots.plot(result::ManhattanPlotResult; kwargs...) = result.figure
Plots.plot(result::QQPlotResult; kwargs...) = result.figure
Plots.plot(result::ForestPlotResult; kwargs...) = result.figure
Plots.plot(result::KaplanMeierResult; kwargs...) = BioToolkit.Clinical.kaplan_meier_plot(result; kwargs...)
Plots.plot(results::AbstractVector{<:KaplanMeierResult}; kwargs...) = _kaplan_meier_figure(results; kwargs...)
Plots.plot(result::ROCResult; kwargs...) = _roc_figure(result; kwargs...)
Plots.plot(result::CIFResult; kwargs...) = _cif_figure(result; kwargs...)
Plots.plot(result::CoxResult; kwargs...) = forest_plot(result)
Plots.plot(result::OncoprintResult; kwargs...) = oncoprint_plot(result; kwargs...)

function _kaplan_meier_figure(results::AbstractVector{<:KaplanMeierResult}; labels=nothing, title::AbstractString="Kaplan-Meier survival", size=(900, 640), show_censors::Bool=true, kwargs...)
    figure = Plots.plot(; xlabel="Time", ylabel="Survival probability", title=title, legend=:bottomleft, size=size, left_margin=18Plots.mm, bottom_margin=18Plots.mm, top_margin=10Plots.mm)
    palette = [:firebrick, :royalblue, :darkgreen, :darkorange, :purple]
    curve_labels = labels === nothing ? ["Group $(index)" for index in eachindex(results)] : String.(labels)
    for (index, result) in enumerate(results)
        x, y = _km_step_data(result)
        color = palette[mod1(index, length(palette))]
        plot!(figure, x, y; seriestype=:steppost, linewidth=2.4, color=color, label=curve_labels[index])
        if show_censors && !isempty(result.censor_times)
            censor_y = [begin
                survival = 1.0
                for (time, current_survival) in zip(result.time, result.survival)
                    time <= t || break
                    survival = current_survival
                end
                survival
            end for t in result.censor_times]
            scatter!(figure, result.censor_times, censor_y; markershape=:x, color=color, markersize=5, label=false)
        end
    end
    return figure
end

end
