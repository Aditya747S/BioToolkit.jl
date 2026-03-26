# BioToolkitTuringExt.jl

## Purpose
This extension binds BioToolkit's source-tracking APIs to Turing.jl models. It supplies the probabilistic model definitions that the core microbiome and metabolomics modules call when Bayesian source tracking is requested.

## Main responsibilities
- Define the Turing model for metabolomics source tracking.
- Define the Turing model for microbiome source tracking.
- Populate the lazy model hooks in the core modules during `__init__`.
- Keep the probabilistic layer optional so the core package can load without Turing.

## Public behavior added by the extension
- `METABOLOMICS_SOURCE_TRACKING_MODEL_IMPL` stores the metabolomics model constructor.
- `SOURCE_TRACKING_MODEL_IMPL` stores the microbiome model constructor.
- `__init__()` registers both model constructors with `BioToolkit.Metabolomics` and `BioToolkit.Microbiome`.
- The extension does not run inference itself; it only provides the model constructors that the core APIs call.

## Model structure
Both models use the same basic pattern:
- infer source proportions with a Dirichlet prior,
- mix the source profiles into one expected composition,
- normalize the mixture,
- observe counts with a multinomial likelihood.

## How it is used
When the extension is loaded, calls to `source_tracking` or `metabolomics_source_tracking` in the core modules can invoke the cached Turing model implementations and sample posterior chains.

The core modules then turn the resulting chain into posterior summaries with their own result containers.

## Why it matters
The Bayesian source-tracking APIs stay optional and isolated from the core package load path. This extension keeps the probabilistic layer modular while still making it available to the modules that need it.