module BioToolkitCairoMakieExt

using BioToolkit
using CairoMakie

import BioToolkit: export_figure, render!, _render_axis!, GenomeBrowser, GenomeViewport, GeneTrack, CoverageTrack, AlignmentTrack

function _browser_figure(browser::GenomeBrowser)
    rows = max(length(browser.tracks), 1)
    height = max(1, sum(track -> track.height, browser.tracks) + browser.gap * max(rows - 1, 0))
    figure = CairoMakie.Figure(resolution=(browser.viewport.pixel_width, height), fontsize=14)
    for (index, track) in enumerate(browser.tracks)
        axis = CairoMakie.Axis(figure[index, 1]; title=index == 1 ? browser.title : "", xlabel=index == rows ? "genomic position" : "", ylabel="track $(index)")
        render!(axis, track, browser.viewport)
    end
    return figure
end

function render!(axis, track::GeneTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(CairoMakie, axis, plan, track)
end

function render!(axis, track::CoverageTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(CairoMakie, axis, plan, track)
end

function render!(axis, track::AlignmentTrack, viewport::GenomeViewport)
    plan = BioToolkit.render_track(track, viewport)
    _render_axis!(CairoMakie, axis, plan, track)
end

function export_figure(browser::GenomeBrowser, path::AbstractString; dpi::Integer=300)
    figure = _browser_figure(browser)
    CairoMakie.save(path, figure)
    return path
end

end