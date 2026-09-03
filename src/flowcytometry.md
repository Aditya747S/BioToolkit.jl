# FlowCytometry Module

Comprehensive flow cytometry and CyTOF analysis toolkit with FCS I/O, gating, compensation, transformation, clustering, and statistical analysis.

## Overview

The `FlowCytometry` module provides a complete pipeline for flow cytometry data analysis:

- **FCS I/O**: Read FCS 2.0/3.0/3.1 files with full metadata support
- **Gating**: Rectangle, quadrant, polygon, and ellipsoid gates with vectorized operations
- **Compensation**: Spillover matrix compensation using numerically stable linear algebra
- **Transformations**: arcsinh, Logicle, inverse Logicle, and hlog transforms
- **Clustering**: FlowSOM (self-organizing maps) and XShift (density-based)
- **Statistics**: Summary statistics, correlation matrices, and provenance tracking

All operations integrate with BioToolkit's provenance system for reproducibility.

---

## Installation

```julia
using BioToolkit
using BioToolkit.FlowCytometry
```

---

## Data Structures

### FlowExperiment

```julia
struct FlowExperiment <: AbstractAnalysisResult
    events::Matrix{Float64}       # events × channels
    channels::Vector{String}      # channel names
    metadata::Dict{String,Any}    # FCS metadata + analysis info
    provenance::ResultProvenance  # full provenance chain
end
```

- `events[i, j]` = intensity of event `i` in channel `j`
- `metadata` includes original FCS keywords plus transformation history

### GateResult

```julia
struct GateResult <: AbstractAnalysisResult
    gate_type::Symbol             # :rectangle, :polygon, :ellipsoid
    indices::Vector{Int}          # event indices passing the gate
    provenance::ResultProvenance
end
```

### QuadrantGateResult

```julia
struct QuadrantGateResult <: AbstractAnalysisResult
    quadrants::NTuple{4, Vector{Int}}  # (Q1, Q2, Q3, Q4) indices
    indices::Vector{Int}               # all gated events
    provenance::ResultProvenance
end
```

Quadrant order:
- Q1 (LL): x < threshold1, y < threshold2
- Q2 (UL): x < threshold1, y ≥ threshold2
- Q3 (UR): x ≥ threshold1, y ≥ threshold2
- Q4 (LR): x ≥ threshold1, y < threshold2

### FlowSOMResult

```julia
struct FlowSOMResult <: AbstractAnalysisResult
    metaclusters::Vector{Int}     # cluster assignment per event
    codes::Matrix{Float64}        # SOM codebook vectors (n_codes × n_features)
    cluster_centers::Matrix{Float64}  # mean expression per metacluster
    provenance::ResultProvenance
end
```

---

## FCS I/O

### read_fcs

```julia
read_fcs(path::String; on_error::Symbol=:throw, mock::Bool=false) -> FlowExperiment
```

Read an FCS file (2.0, 3.0, or 3.1). Supports all datatypes: A (ASCII), F (float32), D (float64), I (integer).

```julia
# Read real file
fcs = read_fcs("sample.fcs")

# Or generate mock data for testing
fcs = read_fcs("nonexistent.fcs"; mock=true)
# Equivalent to:
fcs = mock_flow_experiment(n_events=1000, seed=42)
```

**Parameters:**
- `on_error`: `:throw` (default), `:mock` (return mock data on error)
- `mock`: if `true`, return mock data instead of reading file

### mock_flow_experiment

```julia
mock_flow_experiment(; n_events::Int=1000, 
                     channels::Vector{String}=["FSC-A", "SSC-A", "FITC-A", "PE-A", "APC-A"],
                     seed::Int=1) -> FlowExperiment
```

Generate synthetic flow data with realistic population structure.

---

## Compensation

### compensate_fcs

```julia
compensate_fcs(fcs::FlowExperiment, spillover_matrix::Matrix{Float64}) -> FlowExperiment
```

Apply spillover compensation using numerically stable linear solve (`A \ B` instead of `inv(A) * B`).

```julia
# 5-channel example with 5% spillover
spillover = Matrix(Diagonal(ones(5))) + 0.05 * ones(5, 5)
compensated = compensate_fcs(fcs, spillover)
```

The function validates that the spillover matrix dimensions match the number of channels.

---

## Gating

All gates return a `GateResult` or `QuadrantGateResult` with provenance tracking. Gates are fully vectorized for performance.

### rectangle_gate

```julia
rectangle_gate(fcs::FlowExperiment, 
               channel1::String, min1::Real, max1::Real,
               channel2::String, min2::Real, max2::Real) -> GateResult
```

Select events within rectangular bounds on two channels.

```julia
gate = rectangle_gate(fcs, "FSC-A", 8.0, 12.0, "SSC-A", 6.0, 10.0)
```

### quadrant_gate

```julia
quadrant_gate(fcs::FlowExperiment,
              channel1::String, threshold1::Real,
              channel2::String, threshold2::Real) -> QuadrantGateResult
```

Divide events into four quadrants based on two thresholds.

```julia
qgate = quadrant_gate(fcs, "FSC-A", 10.0, "SSC-A", 8.0)
# qgate.quadrants[1] = Q1 (LL) indices
# qgate.quadrants[2] = Q2 (UL) indices
# qgate.quadrants[3] = Q3 (UR) indices
# qgate.quadrants[4] = Q4 (LR) indices
# qgate.indices = all gated events
```

### polygon_gate

```julia
polygon_gate(fcs::FlowExperiment,
             channel1::String, channel2::String,
             vertices::Vector{Tuple{Float64, Float64}}) -> GateResult
```

Select events within an arbitrary polygon using ray-casting algorithm.

```julia
vertices = [(8.0, 6.0), (12.0, 6.0), (12.0, 10.0), (8.0, 10.0)]
gate = polygon_gate(fcs, "FSC-A", "SSC-A", vertices)
```

Vertices must be ordered (clockwise or counter-clockwise). Minimum 3 vertices required.

### ellipsoid_gate

```julia
ellipsoid_gate(fcs::FlowExperiment,
               channels::Vector{String},
               center::Vector{Float64},
               covariance::Matrix{Float64};
               confidence::Float64=0.95) -> GateResult
```

Select events within a multivariate ellipsoid using Mahalanobis distance. Uses the χ² quantile for the confidence threshold.

```julia
center = [10.0, 8.0, 5.0, 5.0, 5.0]
cov = Matrix(Diagonal([2.0, 2.0, 2.0, 2.0, 2.0]))
gate = ellipsoid_gate(fcs, fcs.channels, center, cov; confidence=0.95)
```

**Parameters:**
- `channels`: vector of channel names defining the ellipsoid dimensions
- `center`: mean vector (length = n_channels)
- `covariance`: covariance matrix (n_channels × n_channels)
- `confidence`: confidence level for χ² threshold (default 0.95)

### apply_gate

```julia
apply_gate(fcs::FlowExperiment, gate::Union{GateResult, QuadrantGateResult}) -> FlowExperiment
```

Apply a gate to create a new FlowExperiment with only the gated events.

```julia
gated_fcs = apply_gate(fcs, gate)
```

---

## Transformations

All transformations return a new `FlowExperiment` with transformed data and updated metadata.

### arcsinh_transform

```julia
arcsinh_transform(fcs::FlowExperiment; 
                  cofactor::Float64=5.0,
                  channels::Vector{String}=fcs.channels) -> FlowExperiment
```

Apply inverse hyperbolic sine (arcsinh) transformation: `asinh(x / cofactor)`.

Standard for CyTOF/mass cytometry data. The cofactor controls the linear region around zero.

```julia
t_fcs = arcsinh_transform(fcs; cofactor=5.0)
```

### logicle_transform

```julia
logicle_transform(fcs::FlowExperiment;
                  channels::Vector{String}=fcs.channels,
                  T::Float64=262144.0,  # top of scale
                  W::Float64=0.5,       # width of linear region
                  M::Float64=4.5,       # decades of log space
                  A::Float64=0.0)       # additional negative range
                  -> FlowExperiment
```

Apply Logicle transformation (Taylor et al., 2001). Standard for traditional flow cytometry.

```julia
l_fcs = logicle_transform(fcs; T=262144.0, W=0.5, M=4.5, A=0.0)
```

Parameters follow BD/FCS convention. The transformation handles negative values smoothly.

### inverse_logicle_transform

```julia
inverse_logicle_transform(fcs::FlowExperiment;
                          channels::Vector{String}=fcs.channels,
                          T::Float64=262144.0,
                          W::Float64=0.5,
                          M::Float64=4.5,
                          A::Float64=0.0) -> FlowExperiment
```

Apply inverse Logicle to recover original scale. Must use same parameters as forward transform.

```julia
original = inverse_logicle_transform(logicle_fcs; T=262144.0, W=0.5, M=4.5, A=0.0)
```

### hlog_transform

```julia
hlog_transform(fcs::FlowExperiment;
               channels::Vector{String}=fcs.channels,
               b::Float64=1.0) -> FlowExperiment
```

Apply hyperbolic log transformation: `asinh(x / (2b)) + log(2b)`.

Smoothly handles negative values and zero. Parameter `b` controls the transition point.

```julia
h_fcs = hlog_transform(fcs; b=1.0)
```

---

## Clustering

### flowsom_cluster

```julia
flowsom_cluster(fcs::FlowExperiment;
                n_metaclusters::Int=10,
                grid_size::Int=10,
                n_epochs::Int=10,
                learning_rate::Float64=0.5,
                initial_radius::Float64=5.0,
                seed::Int=42) -> FlowSOMResult
```

FlowSOM implementation: Self-Organizing Map followed by hierarchical metaclustering.

**Algorithm:**
1. Normalize events (z-score per channel)
2. Initialize grid_size × grid_size codebook vectors
3. Train SOM for n_epochs with decaying learning rate and neighborhood radius
4. Compute pairwise distances between codebook vectors
5. Hierarchically merge nearest clusters until n_metaclusters reached
6. Assign each event to its BMU's metacluster

```julia
som = flowsom_cluster(fcs; 
    n_metaclusters=10,    # final number of metaclusters
    grid_size=10,         # SOM grid dimension (100 codes)
    n_epochs=10,          # training epochs
    learning_rate=0.5,    # initial learning rate
    initial_radius=5.0,   # initial neighborhood radius
    seed=42)              # reproducibility
```

**Returns:** `FlowSOMResult` with:
- `metaclusters`: cluster assignment per event (Int vector)
- `codes`: final SOM codebook (n_codes × n_features)
- `cluster_centers`: mean expression per metacluster

### xshift_cluster

```julia
xshift_cluster(fcs::FlowExperiment;
               K::Int=20,
               approximate::Bool=true,
               seed::Int=42) -> Vector{Int}
```

XShift density-based clustering with approximate kNN for scalability.

**Algorithm:**
1. Build approximate kNN graph using random projection trees (or exact for n ≤ 5000)
2. Compute local density: inverse mean distance to k nearest neighbors
3. Build parent tree: each point links to highest-density neighbor (including self)
4. Find peaks (roots where parent == self)
5. Assign clusters by following parent links to peaks

```julia
# Approximate (recommended for n > 5000)
clusters = xshift_cluster(fcs; K=20, approximate=true, seed=42)

# Exact (only for small datasets)
clusters = xshift_cluster(fcs; K=20, approximate=false, seed=42)
```

**Parameters:**
- `K`: number of nearest neighbors for density estimation
- `approximate`: use random projection trees (O(n log n)) vs exact O(n²)
- `seed`: reproducibility

---

## Statistics

### fcs_stats

```julia
fcs_stats(fcs::FlowExperiment; channels::Vector{String}=fcs.channels) -> Dict{String, Dict}
```

Compute summary statistics per channel.

```julia
stats = fcs_stats(fcs)
# stats["FSC-A"] = Dict(
#   "n" => 1000,
#   "mean" => 14.96,
#   "std" => 2.01,
#   "median" => 14.85,
#   "min" => 8.2,
#   "max" => 22.1,
#   "q25" => 13.5,
#   "q75" => 16.2,
#   "cv" => 0.134
# )
```

### fcs_correlation

```julia
fcs_correlation(fcs::FlowExperiment; channels::Vector{String}=fcs.channels) -> Matrix{Float64}
```

Compute Pearson correlation matrix between channels.

```julia
cor_mat = fcs_correlation(fcs)
# 5×5 Matrix{Float64}
```

---

## Provenance Integration

All functions integrate with BioToolkit's provenance system:

```julia
using BioToolkit: active_provenance_context, get_provenance_context

# Enable automatic tracking
ctx = active_provenance_context() !== nothing ? nothing : BioToolkit.enable_provenance!()

# All operations automatically recorded
fcs = mock_flow_experiment()
gate = rectangle_gate(fcs, "FSC-A", 8.0, 12.0, "SSC-A", 6.0, 10.0)
gated = apply_gate(fcs, gate)
som = flowsom_cluster(fcs)

# Inspect provenance
ctx = get_provenance_context()
for node in values(ctx.nodes)
    println("$(node.operation): $(node.parameters)")
end
```

Each result carries a `ResultProvenance` with:
- Unique ID
- Operation label and source
- Input parameters
- Parent provenance IDs
- Timestamp

---

## Complete Example

```julia
using BioToolkit
using BioToolkit.FlowCytometry
using LinearAlgebra

# 1. Load data (or generate mock)
fcs = mock_flow_experiment(n_events=50000, seed=123)

# 2. Compensate
spillover = Matrix(Diagonal(ones(5))) + 0.05 * ones(5, 5)
fcs_comp = compensate_fcs(fcs, spillover)

# 3. Transform (arcsinh for CyTOF, Logicle for flow)
fcs_trans = arcsinh_transform(fcs_comp; cofactor=5.0)

# 4. Gate lymphocytes (FSC vs SSC)
lymph_gate = rectangle_gate(fcs_trans, "FSC-A", 5.0, 10.0, "SSC-A", 3.0, 8.0)
lymph = apply_gate(fcs_trans, lymph_gate)

# 5. Gate singlets (FSC-H vs FSC-A)
singlet_gate = rectangle_gate(lymph, "FSC-H", 8.0, 12.0, "FSC-A", 5.0, 10.0)
singlets = apply_gate(lymph, singlet_gate)

# 6. Cluster with FlowSOM
som = flowsom_cluster(singlets; 
    n_metaclusters=15, grid_size=10, n_epochs=10, seed=42)

# 7. Analyze clusters
for cl in 1:15
    cl_events = singlets.events[som.metaclusters .== cl, :]
    println("Cluster $cl: $(size(cl_events, 1)) events")
    println("  Medians: ", median(cl_events, dims=1))
end

# 8. Export results
using DataFrames
df = DataFrame(som.metaclusters, :cluster)
for (i, ch) in enumerate(singlets.channels)
    df[!, ch] = singlets.events[:, i]
end
CSV.write("clusters.csv", df)

# 9. Provenance report
ctx = BioToolkit.get_provenance_context()
println(BioToolkit.generate_methods_section(ctx))
```

---

## Performance Tips

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Rectangle gate | O(n) | Fully vectorized |
| Quadrant gate | O(n) | Fully vectorized |
| Polygon gate | O(n·v) | v = vertices |
| Ellipsoid gate | O(n·d²) | d = dimensions |
| Compensation | O(n·d²) | Uses `\` (LU factorization) |
| arcsinh/Logicle | O(n·d) | Vectorized |
| FlowSOM | O(epochs·n·codes) | Grid size affects codes |
| XShift (approx) | O(n log n) | Random projection trees |
| XShift (exact) | O(n²) | Only for n ≤ 5000 |

**Recommendations:**
- Use `approximate=true` for XShift on datasets > 5000 events
- For FlowSOM, `grid_size=10` (100 codes) works well for 5-20 features
- Chain gates: `apply_gate` creates new arrays - consider combining gates
- Transform before clustering (arcsinh/Logicle normalizes dynamic range)

---

## References

- **FCS Format**: ISAC FCS 3.1 Specification
- **Logicle**: Taylor et al., "A unified approach to fluorescence compensation and transformation" (2001)
- **FlowSOM**: Van Gassen et al., "FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data" (2015)
- **XShift**: X-Shift: A general framework for density-based clustering (2015)
- **arcsinh**: Bendall et al., "Single-cell mass cytometry of differential immune and drug responses" (2011)

---

## See Also

- `BioToolkit.read_fcs` - FCS file I/O
- `BioToolkit.compensate_fcs` - Spillover compensation
- `BioToolkit.FlowCytometry.arcsinh_transform` - CyTOF transformation
- `BioToolkit.FlowCytometry.logicle_transform` - Flow transformation
- `BioToolkit.FlowCytometry.flowsom_cluster` - FlowSOM clustering
- `BioToolkit.FlowCytometry.xshift_cluster` - Density-based clustering
- `BioToolkit.provenance` - Provenance tracking system