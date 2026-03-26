# clinical.jl

## Purpose
This file implements the clinical genomics layer of BioToolkit. It combines patient cohort management, survival analysis, Cox modeling, ROC and competing-risks summaries, MAF handling, and several clinically oriented plotting helpers.

## Main structs
- PatientCohort: links a clinical table with genomic data and stable patient IDs.
- KaplanMeierResult: a survival curve with event times, survival probabilities, risk counts, and censoring information.
- CoxTermResult: one term in a Cox proportional hazards model.
- CoxResult: the complete Cox model output with all terms and baseline hazard data.
- MAFRecord: one mutation annotation format record.
- MAFSummary: aggregate mutation counts by sample, gene, and variant class.
- ROCResult: survival ROC summary with thresholds, TPR, FPR, and AUC.
- CIFResult: cumulative incidence function output for competing risks.
- NeuralCoxResult: neural-network-style Cox model output.
- DoseResponseResult: dose-response fitting output with Emax, EC50, Hill coefficient, and fitted values.
- OncoprintResult: mutation matrix and labels for oncoprint-style display.

## Public functions
- read_maf and summarize_maf: parse MAF records and generate summaries.
- tcga_query, tcga_download_files, merge_tcga_count_files, tcga_ingest: data acquisition and ingestion helpers for TCGA-style workflows.
- kaplan_meier: build a KaplanMeierResult from time and status vectors.
- logrank_test: compare survival curves between groups.
- cox_ph: fit Cox proportional hazards models.
- forest_plot: produce a forest plot from model terms.
- survival_roc: compute a time-dependent survival ROC curve.
- cif_curve: compute cumulative incidence curves.
- neural_cox: fit a neural Cox model.
- dose_response_curve: fit dose-response data.
- oncoprint: build a mutation heatmap-style summary.

## What the module does
The module turns clinical and genomic data into interpretable statistical objects. It validates cohort structure, handles survival event coding, supports modeling across patient groups, and returns compact result types that can be plotted or further summarized.

## Typical usage
1. Create a PatientCohort from a clinical DataFrame and a matching genomic matrix.
2. Use kaplan_meier or cox_ph to build survival results.
3. Use logrank_test to compare groups and survival_roc or cif_curve to inspect predictive performance or competing risks.
4. Use read_maf and summarize_maf to inspect mutation burden and variant distribution.
5. Use dose_response_curve or oncoprint for response and alteration summaries.

## Plotting helpers
This file also contains the main clinical plotting API. Kaplan-Meier, forest, ROC, CIF, dose-response, and oncoprint figures are meant to be generated from the corresponding result types, and the plotting extension layer can reuse those results directly.

## Important implementation details
- _clamp_probability keeps p-values and probabilities numerically stable.
- _clinical_index and _cohort_view support patient-level lookup and subsetting.
- _sorted_event_data, _validated_group_levels, and _count_tie_events support survival calculations.
- _km_survival_at and _km_plot_data support Kaplan-Meier plotting and censor alignment.

## Why this file matters
This module is the clinical analysis hub for BioToolkit. It keeps patient-level data, survival models, mutation summaries, and clinical plots aligned so the package can move from cohort input to publication-ready clinical interpretation.
