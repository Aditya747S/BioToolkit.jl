module LongReadBenchmarkFixtures

using BioToolkit

export longread_fixture, evaluate_sv_recovery

@inline _safe_div(num::Real, den::Real) = den == 0 ? 0.0 : Float64(num) / Float64(den)

function _f1(precision::Real, recall::Real)
    denom = precision + recall
    return denom == 0 ? 0.0 : 2.0 * precision * recall / denom
end

function _match_calls(predicted::Vector{StructuralVariant}, truth::AbstractVector{<:NamedTuple}; position_tolerance::Integer=60)
    truth_used = falses(length(truth))
    matched_pairs = Tuple{Int,Int}[]
    false_positive_indices = Int[]

    for (pred_idx, call) in enumerate(predicted)
        best_truth_index = 0
        best_delta = typemax(Int)

        for (truth_idx, target) in enumerate(truth)
            truth_used[truth_idx] && continue
            call.sv_type == target.sv_type || continue
            call.chrom1 == target.chrom1 || continue
            call.chrom2 == target.chrom2 || continue

            delta1 = abs(call.pos1 - target.pos1)
            delta2 = abs(call.pos2 - target.pos2)
            if delta1 <= position_tolerance && delta2 <= position_tolerance
                combined_delta = delta1 + delta2
                if combined_delta < best_delta
                    best_delta = combined_delta
                    best_truth_index = truth_idx
                end
            end
        end

        if best_truth_index > 0
            truth_used[best_truth_index] = true
            push!(matched_pairs, (pred_idx, best_truth_index))
        else
            push!(false_positive_indices, pred_idx)
        end
    end

    false_negative_indices = findall(!, truth_used)
    return matched_pairs, false_positive_indices, false_negative_indices
end

function _type_rows(predicted::Vector{StructuralVariant}, truth::AbstractVector{<:NamedTuple}, matched_pairs::Vector{Tuple{Int,Int}}, false_positive_indices::Vector{Int}, false_negative_indices::Vector{Int})
    types = sort!(collect(union(Symbol[call.sv_type for call in predicted], Symbol[target.sv_type for target in truth])))

    rows = NamedTuple[]
    for sv_type in types
        tp = count(pair -> predicted[pair[1]].sv_type == sv_type, matched_pairs)
        fp = count(idx -> predicted[idx].sv_type == sv_type, false_positive_indices)
        fn = count(idx -> truth[idx].sv_type == sv_type, false_negative_indices)

        precision = _safe_div(tp, tp + fp)
        recall = _safe_div(tp, tp + fn)

        push!(rows, (
            sv_type=String(sv_type),
            tp=tp,
            fp=fp,
            fn=fn,
            precision=precision,
            recall=recall,
            f1=_f1(precision, recall),
        ))
    end

    overall_tp = length(matched_pairs)
    overall_fp = length(false_positive_indices)
    overall_fn = length(false_negative_indices)
    overall_precision = _safe_div(overall_tp, overall_tp + overall_fp)
    overall_recall = _safe_div(overall_tp, overall_tp + overall_fn)

    pushfirst!(rows, (
        sv_type="ALL",
        tp=overall_tp,
        fp=overall_fp,
        fn=overall_fn,
        precision=overall_precision,
        recall=overall_recall,
        f1=_f1(overall_precision, overall_recall),
    ))

    return rows
end

function evaluate_sv_recovery(predicted::Vector{StructuralVariant}, truth::AbstractVector{<:NamedTuple}; position_tolerance::Integer=60)
    matched_pairs, false_positive_indices, false_negative_indices = _match_calls(predicted, truth; position_tolerance=position_tolerance)
    return _type_rows(predicted, truth, matched_pairs, false_positive_indices, false_negative_indices)
end

function longread_fixture()
    records = BamRecord[]

    # Deletion evidence cluster around chr1:~305-385.
    push!(records, BamRecord(
        "del_a",
        "chr1",
        250,
        [BamCigarOp(50, 'M'), BamCigarOp(80, 'D'), BamCigarOp(50, 'M')],
        repeat("A", 100);
        mapq=60,
        mate_refname="chr1",
        mate_pos=610,
        template_length=420,
    ))
    push!(records, BamRecord(
        "del_b",
        "chr1",
        257,
        [BamCigarOp(50, 'M'), BamCigarOp(80, 'D'), BamCigarOp(50, 'M')],
        repeat("A", 100);
        mapq=60,
        mate_refname="chr1",
        mate_pos=617,
        template_length=430,
    ))

    # Inversion evidence via same orientation mates on chr1.
    push!(records, BamRecord(
        "inv_a",
        "chr1",
        700,
        [BamCigarOp(100, 'M')],
        repeat("C", 100);
        flag=0x30,
        mapq=60,
        mate_refname="chr1",
        mate_pos=840,
        template_length=220,
    ))
    push!(records, BamRecord(
        "inv_b",
        "chr1",
        740,
        [BamCigarOp(100, 'M')],
        repeat("C", 100);
        flag=0x30,
        mapq=60,
        mate_refname="chr1",
        mate_pos=880,
        template_length=240,
    ))

    # Translocation evidence across chr1 -> chr2.
    push!(records, BamRecord(
        "tra_a",
        "chr1",
        1200,
        [BamCigarOp(100, 'M')],
        repeat("G", 100);
        mapq=55,
        mate_refname="chr2",
        mate_pos=400,
        template_length=0,
    ))
    push!(records, BamRecord(
        "tra_b",
        "chr1",
        1208,
        [BamCigarOp(100, 'M')],
        repeat("G", 100);
        mapq=55,
        mate_refname="chr2",
        mate_pos=408,
        template_length=0,
    ))

    truth = [
        (sv_type=:DEL, chrom1="chr1", pos1=305, chrom2="chr1", pos2=385),
        (sv_type=:INV, chrom1="chr1", pos1=721, chrom2="chr1", pos2=861),
        (sv_type=:TRA, chrom1="chr1", pos1=1205, chrom2="chr2", pos2=405),
    ]

    return (
        records=records,
        truth=truth,
        call_kwargs=(
            min_sv_size=60,
            min_softclip=30,
            insert_size_mean=500.0,
            insert_size_std=100.0,
            z_threshold=3.0,
            cluster_window=90,
            min_support=2,
        ),
    )
end

end
