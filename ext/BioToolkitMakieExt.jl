module BioToolkitMakieExt

using BioToolkit
using Makie

import BioToolkit: render!, _render_axis!, GenomeBrowser, GenomeViewport, GeneTrack, CoverageTrack, AlignmentTrack, SingleCellViewer, interactive_singlecell_viewer, cell_hover_text, recluster_singlecell_viewer!

function _browser_figure(browser::GenomeBrowser)
    rows = max(length(browser.tracks), 1)
    height = max(1, sum(track -> track.height, browser.tracks) + browser.gap * max(rows - 1, 0))
    figure = Makie.Figure(resolution=(browser.viewport.pixel_width, height), fontsize=14)
    for (index, track) in enumerate(browser.tracks)
        axis = Makie.Axis(figure[index, 1]; title=index == 1 ? browser.title : "", xlabel=index == rows ? "genomic position" : "", ylabel="track $(index)")
        render!(axis, track, browser.viewport)
    end
    return figure
end

function render!(axis, track::GeneTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(Makie, axis, plan, track)
end

function render!(axis, track::CoverageTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(Makie, axis, plan, track)
end

function render!(axis, track::AlignmentTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(Makie, axis, plan, track)
end

function Base.display(browser::GenomeBrowser)
    display(_browser_figure(browser))
end

function _viewer_figure(viewer::SingleCellViewer)
    coordinates = viewer.embedding[:, 1:2]
    label_state = Observable(copy(viewer.cluster_labels))
    figure = Makie.Figure(resolution=(900, 700), fontsize=14)
    axis = Makie.Axis(figure[1, 1]; title=viewer.title, xlabel="embedding 1", ylabel="embedding 2")
    scatter = Makie.scatter!(axis, coordinates[:, 1], coordinates[:, 2]; color=label_state, markersize=8, inspectable=true, inspector_label = (plot, index, position) -> cell_hover_text(viewer, index))
    slider_label = Makie.Label(figure[2, 1], "resolution k = $(viewer.k)")
    slider = Makie.Slider(figure[3, 1], range=5:1:max(5, min(50, size(viewer.embedding, 1))), startvalue=viewer.k)
    Makie.on(slider.value) do value
        k = Int(round(value))
        recluster_singlecell_viewer!(viewer; k=k)
        label_state[] = copy(viewer.cluster_labels)
        slider_label.text = "resolution k = $(k)"
    end
    return figure
end

function Base.display(viewer::SingleCellViewer)
    display(_viewer_figure(viewer))
end

end