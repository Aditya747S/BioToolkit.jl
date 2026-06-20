using Test, DataFrames, BioToolkit

println("=== Thorough Type-Safety & Edge Case Analysis ===")

# 1. register_container_provenance! with Symbol-keyed metadata and PROVENANCE_HASH_KEY
println("\n1. register_container_provenance! Symbol dict + hash key")
try
    struct SymMetaContainer
        metadata::Dict{Symbol,Any}
    end
    c = SymMetaContainer(Dict{Symbol,Any}())
    ctx = BioToolkit.ProvenanceContext()
    id = BioToolkit.register_container_provenance!(ctx, c, "test_op";
        provenance_hash="hashval")
    @assert id !== nothing "Should return an id"
    # Check the hash was stored with the right key type
    @assert haskey(c.metadata, :provenance_hash) "Should have Symbol key"
    @assert c.metadata[:provenance_hash] == "hashval"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 2. metadata_provenance with Symbol vs String key
println("\n2. metadata_provenance key type handling")
try
    # Symbol key
    meta1 = Dict{Symbol,Any}(:provenance => BioToolkit.provenance_record("label1", "src1"))
    p1 = BioToolkit.metadata_provenance(meta1)
    @assert p1 !== nothing "Should find Symbol key"
    @assert p1.label == "label1"
    
    # String key
    meta2 = Dict{String,Any}("provenance" => BioToolkit.provenance_record("label2", "src2"))
    p2 = BioToolkit.metadata_provenance(meta2)
    # PROVENANCE_METADATA_KEY is "provenance" (a String)
    # The function checks both PROVENANCE_METADATA_KEY and Symbol(PROVENANCE_METADATA_KEY)
    # For Dict{String,Any}, haskey(meta2, "provenance") should work
    # But haskey(meta2, :provenance) should return false
    @assert p2 !== nothing "Should find String key"
    @assert p2.label == "label2"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 3. _metadata_key with various dict key types
println("\n3. _metadata_key correctness")
try
    sym_dict = Dict{Symbol,Any}()
    str_dict = Dict{String,Any}()
    any_dict = Dict{Any,Any}()
    
    # With Symbol input
    @assert BioToolkit._metadata_key(sym_dict, :test) === :test
    @assert BioToolkit._metadata_key(str_dict, :test) == "test"
    
    # With String input
    @assert BioToolkit._metadata_key(sym_dict, "test") === :test
    @assert BioToolkit._metadata_key(str_dict, "test") == "test"
    
    # With Any key type - should use Symbol since Any <: Symbol is false
    @assert BioToolkit._metadata_key(any_dict, :test) == "test"
    @assert BioToolkit._metadata_key(any_dict, "test") == "test"
    
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 4. container_provenance_id with both Symbol and String PROVENANCE_ID_KEY attempts
println("\n4. container_provenance_id with Dict{Any,Any}")
try
    struct AnyMetaContainer
        metadata::Dict{Any,Any}
    end
    # Test with Symbol key
    c1 = AnyMetaContainer(Dict{Any,Any}(:provenance_id => "id_sym"))
    id1 = BioToolkit.container_provenance_id(c1)
    println("  Symbol key result: $id1")
    
    # Test with String key
    c2 = AnyMetaContainer(Dict{Any,Any}("provenance_id" => "id_str"))
    id2 = BioToolkit.container_provenance_id(c2)
    println("  String key result: $id2")
    
    @assert id1 == "id_sym" || id2 == "id_str" "At least one should work"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 5. ensure_provenance_id! with Dict{String,Any} that has Symbol key (impossible but good to check)
println("\n5. ensure_provenance_id! with empty Dict generates new ID")
try
    meta = Dict{String,Any}()
    id = BioToolkit.ensure_provenance_id!(meta)
    @assert id !== nothing
    @assert length(id) == 64 "Should be SHA-256 hex"
    @assert haskey(meta, "provenance_id") "Key should be String type"
    println("  PASS: id=$(first(id, 16))...")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 6. provenance_record from Dict with non-standard status types
println("\n6. provenance_record from Dict with various value types")
try
    # Status as string instead of symbol
    payload = Dict{String,Any}(
        "label" => "test",
        "source" => "src",
        "status" => "warn",  # String, not Symbol
        "warnings" => "single warning string",  # String, not Vector
        "parameters" => Dict("k" => "v")
    )
    prov = BioToolkit.provenance_record(payload)
    @assert prov.label == "test"
    @assert prov.status === :warn "Should convert String to Symbol"
    @assert length(prov.warnings) == 1 "Should wrap single string in vector"
    @assert prov.warnings[1] == "single warning string"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 7. _provenance_parameters with vector of pairs
println("\n7. _provenance_parameters_dict with Vector{Pair}")
try
    pairs_vec = ["method" => "zscore", "n" => 5]
    d = BioToolkit._provenance_parameters_dict(pairs_vec)
    @assert d["method"] == "zscore"
    @assert d["n"] == 5
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 8. provenance_structural_hash with Nothing
println("\n8. provenance_structural_hash(nothing)")
try
    h = BioToolkit.provenance_structural_hash(nothing)
    @assert length(h) == 64
    # nothing is not a Number, not a String, not Array, not Dict, not DataFrame
    # Falls through to `string(typeof(data), ':', objectid(data))`
    println("  PASS: hash=$(first(h, 16))...")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 9. provenance_content_hash(nothing)
println("\n9. provenance_content_hash(nothing)")
try
    h = BioToolkit.provenance_content_hash(nothing)
    @assert length(h) == 64
    println("  PASS: hash=$(first(h, 16))...")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 10. Mermaid label with extremely long node ID (orphan parent)
println("\n10. provenance_to_mermaid with very long orphan parent ID")
try
    ctx = BioToolkit.ProvenanceContext()
    long_id = repeat("a", 100)
    n = BioToolkit.ProvenanceNode("child", "child_op", Dict{String,Any}(), [long_id], "t1")
    ctx.nodes["child"] = n
    mermaid = BioToolkit.provenance_to_mermaid(ctx)
    # Should truncate to first 16 chars + "…"
    @assert occursin("…", mermaid) "Should have truncation marker"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 11. with_provenance on a Vector of items with some having metadata
println("\n11. with_provenance on mixed Vector")
try
    BioToolkit.disable_provenance!()
    struct MetaItem
        metadata::Dict{Symbol,Any}
        value::Int
    end
    items = [MetaItem(Dict{Symbol,Any}(), 1), MetaItem(Dict{Symbol,Any}(), 2)]
    BioToolkit.with_provenance(items, "batch", "src"; parameters=(n=2,))
    # Each item should get its own provenance
    p1 = BioToolkit.metadata_provenance(items[1].metadata)
    p2 = BioToolkit.metadata_provenance(items[2].metadata)
    @assert p1 !== nothing "Item 1 should have provenance"
    @assert p2 !== nothing "Item 2 should have provenance"
    # IDs should be different (each gets its own)
    id1 = BioToolkit.container_provenance_id(items[1])
    id2 = BioToolkit.container_provenance_id(items[2])
    @assert id1 != id2 "Each item should have unique ID"
    println("  PASS: id1=$(first(id1, 8))..., id2=$(first(id2, 8))...")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 12. flatten_provenance! on DataFrame without provenance
println("\n12. flatten_provenance! on untracked DataFrame")
try
    df = DataFrame(x=[1,2])
    BioToolkit.flatten_provenance!(df)
    # Should not add columns since there's no provenance
    @assert !("provenance_id" in names(df)) "Should not add provenance_id"
    @assert !("provenance_hash" in names(df)) "Should not add provenance_hash"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 13. Multiple calls to stamp_provenance! should update, not duplicate
println("\n13. stamp_provenance! idempotency")
try
    meta = Dict{Symbol,Any}()
    BioToolkit.stamp_provenance!(meta; label="first", source="s1")
    id1 = BioToolkit.container_provenance_id(meta)
    
    BioToolkit.stamp_provenance!(meta; label="second", source="s2")
    id2 = BioToolkit.container_provenance_id(meta)
    
    prov = BioToolkit.metadata_provenance(meta)
    @assert prov.label == "second" "Should overwrite label"
    @assert prov.source == "s2" "Should overwrite source"
    println("  PASS: id1=$(first(id1, 8))..., id2=$(first(id2, 8))...")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 14. _scoped_provenance_stack returns same vector across calls in same task
println("\n14. _scoped_provenance_stack consistency")
try
    s1 = BioToolkit._scoped_provenance_stack()
    s2 = BioToolkit._scoped_provenance_stack()
    @assert s1 === s2 "Should be the same Vector"
    println("  PASS")
catch e
    println("  FAIL: ", sprint(showerror, e))
end

# 15. Copy semantics — ProvenanceNode is immutable struct, copying context should be safe
println("\n15. ProvenanceNode immutability verification")
try
    n = BioToolkit.ProvenanceNode("id", "op", Dict{String,Any}("k" => "v"), ["p1"], "t")
    # ProvenanceNode fields are immutable (except parent_ids which is a Vector)
    # But parent_ids can be mutated in-place (Vector is mutable)
    ctx = BioToolkit.ProvenanceContext()
    ctx.nodes["id"] = n
    copied = copy(ctx)
    
    # Mutating parent_ids in original should NOT affect copy... wait,
    # copy() uses OrderedDict copy which shallow-copies values
    # So the ProvenanceNode in copied.nodes["id"] is the SAME object
    push!(n.parent_ids, "p2")
    
    # Check if the copy was affected
    copied_node = copied.nodes["id"]
    if "p2" in copied_node.parent_ids
        println("  WARNING: copy() does shallow copy of ProvenanceNode - parent_ids mutation leaks!")
        println("  This is a BUG in copy(ProvenanceContext)")
    else
        println("  PASS: copy() properly isolates ProvenanceNodes")
    end
catch e
    println("  FAIL: ", sprint(showerror, e))
end

println("\n=== All thorough analysis tests completed ===")
