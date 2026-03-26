const _CUDA_SEQUENCE_LOADED = Ref(false)
const _CUDA_SEQUENCE_HAMMING_IMPL = Ref{Any}(nothing)
const _CUDA_SEQUENCE_KMER_IMPL = Ref{Any}(nothing)
const _CUDA_SEQUENCE_DOTMATRIX_IMPL = Ref{Any}(nothing)

const _CUDA_QUERY_LOADED = Ref(false)
const _CUDA_QUERY_BIN_IMPL = Ref{Any}(nothing)
const _CUDA_QUERY_WINDOW_IMPL = Ref{Any}(nothing)

const _CUDA_PHYLO_LOADED = Ref(false)
const _CUDA_PHYLO_DISTANCE_MATRIX_IMPL = Ref{Any}(nothing)

_is_cuda_backed_array(x) = nameof(typeof(x)) == :CuArray && nameof(parentmodule(typeof(x))) == :CUDA

"""
    _ensure_cuda_sequence!()

Load and cache the CUDA implementations used by sequence helpers.
"""
function _ensure_cuda_sequence!()
    if !_CUDA_SEQUENCE_LOADED[]
        @eval function _cuda_sequence_hamming_impl(left::CUDA.CuArray{UInt8,1}, right::CUDA.CuArray{UInt8,1})
            length(left) == length(right) || throw(ArgumentError("sequences must have the same length"))
            return Int(sum(left .!= right))
        end

        @eval function _cuda_sequence_kmer_impl(bytes::CUDA.CuArray{UInt8,1}, k::Integer)
            k <= 0 && throw(ArgumentError("k must be positive"))
            length(bytes) < k && return Dict{String,Int}()
            k <= 10 || throw(ArgumentError("CUDA k-mer frequency supports k <= 10"))

            total_kmers = Int(5)^Int(k)
            counts = CUDA.zeros(Int32, total_kmers)
            threads = 256
            blocks = cld(length(bytes) - Int(k) + 1, threads)

            function _kmer_frequency_kernel!(counts, bytes, k::Int32)
                index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
                last_start = length(bytes) - Int(k) + 1
                if index <= last_start
                    current_key = 0
                    @inbounds for offset in 0:(Int(k) - 1)
                        code = _dna_code_byte(bytes[index + offset])
                        code == UInt8(255) && return
                        current_key = current_key * 5 + Int(code)
                    end
                    CUDA.@atomic counts[current_key + 1] += Int32(1)
                end
                return
            end

            CUDA.@sync @cuda threads=threads blocks=blocks _kmer_frequency_kernel!(counts, bytes, Int32(k))

            host_counts = Array(counts)
            result = Dict{String,Int}()
            @inbounds for (index, count) in pairs(host_counts)
                count == 0 && continue
                result[_decode_kmer_key(index, Int(k))] = Int(count)
            end

            return result
        end

        @eval function _cuda_sequence_dotmatrix_impl(seq1::CUDA.CuArray{UInt8,1}, seq2::CUDA.CuArray{UInt8,1})
            m = Int32(length(seq1))
            n = Int32(length(seq2))
            matrix = CUDA.zeros(Int8, m, n)
            total = m * n
            threads = 256
            blocks = cld(Int(total), threads)

            function _dotmatrix_kernel!(matrix, seq1, seq2, m::Int32, n::Int32)
                idx = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
                if idx <= m * n
                    j = div(idx - Int32(1), m) + Int32(1)
                    i = rem(idx - Int32(1), m) + Int32(1)
                    matrix[i, j] = Int8((seq1[i] | 0x20) == (seq2[j] | 0x20))
                end
                return
            end

            CUDA.@sync @cuda threads=threads blocks=blocks _dotmatrix_kernel!(matrix, seq1, seq2, m, n)
            return matrix
        end

        _CUDA_SEQUENCE_HAMMING_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_sequence_hamming_impl"))
        _CUDA_SEQUENCE_KMER_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_sequence_kmer_impl"))
        _CUDA_SEQUENCE_DOTMATRIX_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_sequence_dotmatrix_impl"))
        _CUDA_SEQUENCE_LOADED[] = true
    end
    return nothing
end

"""
    _ensure_cuda_query!()

Load and cache the CUDA implementations used by query helpers.
"""
function _ensure_cuda_query!()
    if !_CUDA_QUERY_LOADED[]
        @eval function _cuda_query_bin_impl(positions::CUDA.CuArray{<:Integer,1}, bin_size::Integer)
            isempty(positions) && return Dict{Int,Int}()

            bin_count = Int(fld(Int(maximum(positions)), Int(bin_size))) + 1
            histogram = CUDA.zeros(Int32, bin_count)
            threads = 256
            blocks = cld(length(positions), threads)

            function _bin_positions_kernel!(histogram, positions, bin_size::Int32)
                index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
                if index <= length(positions)
                    position = Int32(positions[index])
                    bin_index = Int32(fld(position, bin_size)) + Int32(1)
                    CUDA.@atomic histogram[bin_index] += Int32(1)
                end
                return
            end

            CUDA.@sync @cuda threads=threads blocks=blocks _bin_positions_kernel!(histogram, positions, Int32(bin_size))

            counts = Dict{Int,Int}()
            host_histogram = Array(histogram)
            @inbounds for (index, count) in pairs(host_histogram)
                count == 0 && continue
                counts[index - 1] = Int(count)
            end

            return counts
        end

        @eval function _cuda_query_window_impl(starts::CUDA.CuArray{<:Integer,1}, stops::CUDA.CuArray{<:Integer,1}, window_size::Integer)
            length(starts) == length(stops) || throw(ArgumentError("starts and stops must have the same length"))
            isempty(starts) && return Dict{Int,Int}()

            max_stop = Int(maximum(stops))
            max_window = fld(max_stop - 1, window_size)
            coverage = CUDA.zeros(Int32, max_window + 1)
            threads = 256
            blocks = cld(length(starts), threads)

            function _window_coverage_kernel!(coverage, starts, stops, window_size::Int32)
                index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
                if index <= length(starts)
                    start = Int(starts[index])
                    stop = Int(stops[index])
                    if stop > start
                        first_window = fld(start, Int(window_size))
                        last_window = fld(stop - 1, Int(window_size))

                        for window_index in first_window:last_window
                            window_start = window_index * Int(window_size)
                            window_stop = window_start + Int(window_size)
                            overlap_start = max(start, window_start)
                            overlap_stop = min(stop, window_stop)
                            overlap = overlap_stop - overlap_start
                            if overlap > 0
                                CUDA.@atomic coverage[window_index + 1] += Int32(overlap)
                            end
                        end
                    end
                end
                return
            end

            CUDA.@sync @cuda threads=threads blocks=blocks _window_coverage_kernel!(coverage, starts, stops, Int32(window_size))

            host_coverage = Array(coverage)
            result = Dict{Int,Int}()
            @inbounds for (index, count) in pairs(host_coverage)
                count == 0 && continue
                result[index - 1] = Int(count)
            end

            return result
        end

        _CUDA_QUERY_BIN_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_query_bin_impl"))
        _CUDA_QUERY_WINDOW_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_query_window_impl"))
        _CUDA_QUERY_LOADED[] = true
    end
    return nothing
end

"""
    _ensure_cuda_phylo!()

Load and cache the CUDA implementation used by phylogeny helpers.
"""
function _ensure_cuda_phylo!()
    if !_CUDA_PHYLO_LOADED[]
        @eval function _cuda_phylo_distance_matrix_impl(sequences::Vector{<:AbstractString}; method=:hamming)
            N = length(sequences)
            D = zeros(Float64, N, N)

            if (method == :hamming || method == :p_distance) && N > 0
                len = length(sequences[1])
                devices = collect(CUDA.devices())
                num_devices = length(devices)
                if num_devices > 0
                    Threads.@threads for dev_idx in 1:num_devices
                        CUDA.device!(dev_idx - 1)
                        chunk_size = cld(N, num_devices)
                        start_i = (dev_idx - 1) * chunk_size + 1
                        end_i = min(dev_idx * chunk_size, N)

                        for i in start_i:end_i
                            seq_i_cu = CUDA.CuArray{UInt8}(Vector{UInt8}(sequences[i]))
                            for j in (i+1):N
                                seq_j_cu = CUDA.CuArray{UInt8}(Vector{UInt8}(sequences[j]))
                                mismatches = CUDA.sum(seq_i_cu .!= seq_j_cu)
                                val = Float64(mismatches) / len
                                D[i, j] = val
                                D[j, i] = val
                            end
                        end
                    end
                    return D
                end
            end

            return D
        end

        _CUDA_PHYLO_DISTANCE_MATRIX_IMPL[] = Base.invokelatest(getfield, @__MODULE__, Symbol("_cuda_phylo_distance_matrix_impl"))
        _CUDA_PHYLO_LOADED[] = true
    end
    return nothing
end