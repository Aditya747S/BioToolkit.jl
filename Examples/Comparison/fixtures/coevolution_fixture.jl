module CoevolutionBenchmarkFixtures

using Random
using Statistics
using BioToolkit

export coevolution_fixture, evaluate_contact_recovery

function _canonical_pair(i::Integer, j::Integer)
    return i <= j ? (Int(i), Int(j)) : (Int(j), Int(i))
end

function _collect_scored_pairs(scores::AbstractMatrix{<:Real}; min_separation::Integer=5)
    pairs = Tuple{Int,Int,Float64}[]
    n_positions = size(scores, 1)

    for i in 1:n_positions-1
        for j in i+1:n_positions
            abs(i - j) < min_separation && continue
            score = Float64(scores[i, j])
            score <= 0 && continue
            push!(pairs, (i, j, score))
        end
    end

    sort!(pairs; by=pair -> pair[3], rev=true)
    return pairs
end

function _truth_rank_recovery(ranked_pairs::Vector{Tuple{Int,Int,Float64}}, truth_pairs::Vector{Tuple{Int,Int}})
    ranked_lookup = Dict{Tuple{Int,Int},Int}()
    for (rank, pair) in enumerate(ranked_pairs)
        ranked_lookup[_canonical_pair(pair[1], pair[2])] = rank
    end

    n_ranked = max(length(ranked_pairs), 1)
    rank_scores = Float64[]

    for pair in truth_pairs
        key = _canonical_pair(pair[1], pair[2])
        rank = get(ranked_lookup, key, n_ranked)
        percentile = 1.0 - (rank - 1) / n_ranked
        push!(rank_scores, percentile)
    end

    return isempty(rank_scores) ? 0.0 : mean(rank_scores)
end

function _mean_true_pair_distance(structure::Structure, truth_pairs::Vector{Tuple{Int,Int}})
    model = structure.models[1]
    chain = model.chains[1]

    distances = Float64[]
    for (left, right) in truth_pairs
        atom_left = chain.residues[left].atoms[1]
        atom_right = chain.residues[right].atoms[1]
        distance = sqrt((atom_left.x - atom_right.x)^2 + (atom_left.y - atom_right.y)^2 + (atom_left.z - atom_right.z)^2)
        push!(distances, distance)
    end

    return isempty(distances) ? Inf : mean(distances)
end

function evaluate_contact_recovery(contact_map::ContactMap, truth_pairs::Vector{Tuple{Int,Int}}; top_n::Integer=length(truth_pairs), min_separation::Integer=5)
    truth_set = Set(_canonical_pair(pair[1], pair[2]) for pair in truth_pairs)
    ranked_pairs = _collect_scored_pairs(contact_map.scores; min_separation=min_separation)

    top_n_effective = min(Int(top_n), length(ranked_pairs))
    predicted_top = ranked_pairs[1:top_n_effective]
    predicted_set = Set(_canonical_pair(pair[1], pair[2]) for pair in predicted_top)

    tp = length(intersect(truth_set, predicted_set))
    fp = max(length(predicted_set) - tp, 0)
    fn = max(length(truth_set) - tp, 0)

    precision = length(predicted_set) == 0 ? 0.0 : tp / length(predicted_set)
    recall = length(truth_set) == 0 ? 0.0 : tp / length(truth_set)
    f1 = (precision + recall) == 0 ? 0.0 : 2.0 * precision * recall / (precision + recall)

    true_scores = Float64[]
    background_scores = Float64[]
    for pair in ranked_pairs
        key = _canonical_pair(pair[1], pair[2])
        if key in truth_set
            push!(true_scores, pair[3])
        else
            push!(background_scores, pair[3])
        end
    end

    signal_gap = (isempty(true_scores) || isempty(background_scores)) ? 0.0 : mean(true_scores) - mean(background_scores)
    rank_recovery = _truth_rank_recovery(ranked_pairs, truth_pairs)

    structure = fold_from_contacts(contact_map; top_n=max(top_n_effective, length(truth_pairs) * 2), iterations=900, learning_rate=0.02, seed=13)
    mean_true_pair_distance = _mean_true_pair_distance(structure, truth_pairs)

    return (
        tp=tp,
        fp=fp,
        fn=fn,
        precision=precision,
        recall=recall,
        f1=f1,
        signal_gap=signal_gap,
        rank_recovery=rank_recovery,
        mean_true_pair_distance=mean_true_pair_distance,
    )
end

function coevolution_fixture(; seed::Integer=2026, n_sequences::Integer=240, sequence_length::Integer=36, noise_rate::Real=0.12)
    Random.seed!(seed)

    alphabet = collect("ACDEFGHIKLMNPQRSTVWY")
    truth_pairs = [(5, 24), (10, 28), (14, 33)]

    coupling_states = [
        [('A', 'V'), ('G', 'L'), ('S', 'T'), ('D', 'E')],
        [('F', 'Y'), ('I', 'M'), ('N', 'Q'), ('R', 'K')],
        [('C', 'W'), ('H', 'Y'), ('P', 'L'), ('T', 'S')],
    ]

    sequences = String[]
    for _ in 1:n_sequences
        chars = [rand(alphabet) for _ in 1:sequence_length]

        for pair_index in eachindex(truth_pairs)
            left, right = truth_pairs[pair_index]
            state_idx = rand(eachindex(coupling_states[pair_index]))
            expected_left, expected_right = coupling_states[pair_index][state_idx]

            if rand() < noise_rate
                chars[left] = rand(alphabet)
                chars[right] = rand(alphabet)
            else
                chars[left] = expected_left
                chars[right] = expected_right
            end
        end

        push!(sequences, String(chars))
    end

    msa = MultipleSequenceAlignment(sequences)
    return (
        alignment=msa,
        truth_pairs=truth_pairs,
    )
end

end
