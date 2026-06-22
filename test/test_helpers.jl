# Shared test helpers — included exactly once in runtests.jl before any other
# test file.  Centralising helpers here prevents the "Method definition
# overwritten" warning that occurs when multiple test files each define the
# same helper at top level in Main.

using Test
using BioToolkit

"""
    _node_by_operation(ctx, operation)

Return the first `ProvenanceNode` in `ctx` whose `operation` field equals
`operation`.  Asserts (via `@test`) that at least one such node exists.
"""
function _node_by_operation(ctx::ProvenanceContext, operation::AbstractString)
    matches = [node for node in values(ctx.nodes) if node.operation == operation]
    @test !isempty(matches)
    return first(matches)
end

"""
    _write_text_file(path, text)

Write `text` to `path` and return `path`.
"""
function _write_text_file(path::AbstractString, text::AbstractString)
    open(path, "w") do io
        write(io, text)
    end
    return path
end
