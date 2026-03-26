# BioToolkitMakieExt.jl

## Purpose
This extension provides Makie-based rendering for the genome-browser layer. It lets the browser render plans be displayed through Makie without making Makie a hard dependency of the core package.

## Main responsibilities
- Build a Makie `Figure` for a `GenomeBrowser`.
- Dispatch `render!` for gene, coverage, and alignment tracks.
- Reuse the core `render_track` logic from BioToolkit while delegating drawing to Makie axes.
- Keep interactive rendering in a thin backend layer so the browser planner can remain backend-neutral.

## Public behavior added by the extension
- `render!(axis, track::GeneTrack, viewport::GenomeViewport)` renders gene-track plans.
- `render!(axis, track::CoverageTrack, viewport::GenomeViewport)` renders coverage plans.
- `render!(axis, track::AlignmentTrack, viewport::GenomeViewport)` renders alignment plans.
- `Base.display(browser::GenomeBrowser)` displays the browser as a Makie figure.
- The extension uses the same track-plan objects that the core browser module computes, so the visual output stays consistent with the non-Makie render path.

## Internal helper
- `_browser_figure(browser)` creates a figure sized to the browser viewport and track stack.

## How it is used
When Makie is installed, a browser object can be displayed directly and the extension will render each track into its own axis. The core browser planning code still decides what should be drawn; this extension only handles the Makie surface.

Typical usage is to construct a `GenomeBrowser` from one or more tracks, let the core module compute the render plan, and then simply display the browser object.

## Why it matters
This keeps the browser backend-agnostic. BioToolkit can prepare the same render plan once and then expose it through Makie or other supported backends.