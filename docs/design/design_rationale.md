# Design Rationale: Neuro-Preprocessing Pipeline

This document captures the key design decisions, their justifications, and important context for the neuro-preprocessing pipeline. It is intended to serve as a companion to the technical documentation and data dictionaries, providing the "why" behind the "what."

---

## Table of Contents

1. [Overall Architecture](#overall-architecture)
2. [Spike-Centric Design Philosophy](#spike-centric-design-philosophy)
3. [Manifest-Driven Workflow](#manifest-driven-workflow)
4. [Channel Ordering and Probe Validation](#channel-ordering-and-probe-validation)
5. [Trial Alignment: The Anchor-and-Step Algorithm](#trial-alignment-the-anchor-and-step-algorithm)
6. [Timestamp Correction for Legacy Tasks](#timestamp-correction-for-legacy-tasks)
7. [Dynamic Field Discovery from PLDAPS](#dynamic-field-discovery-from-pldaps)
8. [Waveform Extraction Strategy](#waveform-extraction-strategy)
9. [Error Handling Philosophy](#error-handling-philosophy)
10. [Known Issues and Future Considerations](#known-issues-and-future-considerations)

---

## Overall Architecture

### The Five-Stage Pipeline

The preprocessing pipeline transforms raw electrophysiology data into analysis-ready `session_data.mat` files through five sequential stages:

| Stage | Function | Input | Output |
|-------|----------|-------|--------|
| 1. Spike Data Prep | `prep.prepare_spikes_for_kilosort` | `.ns5` (raw broadband) | `.dat` (Kilosort-ready binary) |
| 2. Behavioral Data Prep | `prep.prepare_behavioral_data` | `.nev` + PLDAPS `.mat` | `_intermediate_data.mat` |
| 3. Spike Sorting | Kilosort + Phy (manual) | `.dat` | `spike_times.npy`, `spike_clusters.npy`, etc. |
| 4. Waveform Extraction | `consolidate.extract_waveforms` | `.dat` + Kilosort output | Preliminary `_session_data.mat` |
| 5. Final Consolidation | `consolidate.consolidate_data` | All intermediate files | Final `_session_data.mat` |

### Rationale

This staged approach provides several benefits:

- **Modularity**: Each stage can be re-run independently if issues are discovered.
- **Checkpointing**: The manifest tracks completion status, allowing the pipeline to resume after interruption.
- **Transparency**: Intermediate files allow inspection at each stage for debugging.
- **Accommodation of manual steps**: Stage 3 (spike sorting) requires human curation; the pipeline gracefully waits for this by checking for Kilosort output files.

---

## Spike-Centric Design Philosophy

### Decision

The pipeline processes only spiking activity. Local Field Potentials (LFPs) are not extracted, processed, or included in the output.

### Rationale

Spikes are the fundamental unit of neural communication—the "lingua franca" of the central nervous system. When well-isolated single-unit activity is available (as it is with modern high-density probes and careful spike sorting), it provides more direct and interpretable information about neural computation than LFPs.

LFPs represent aggregate activity and are most useful when single-unit isolation is poor or when studying population-level oscillatory dynamics. For the research questions addressed by this pipeline (attention, eye movements, subcortical circuits), single-unit activity is the appropriate level of analysis.

### Implications

- The `.ns5` file (30 kHz broadband) is read only for spike sorting; lower-sample-rate `.ns2`/`.ns3` files are not processed.
- No bandpass filtering for LFP extraction is performed.
- The `session_data.mat` output contains no LFP-related fields.

### Future Considerations

If a specific scientific need for LFP analysis arises, the pipeline could be extended. The raw data remains available, and a parallel LFP extraction path could be added without disrupting the existing spike-focused workflow.

---

## Manifest-Driven Workflow

### Decision

All pipeline operations are controlled by a single CSV file (`session_manifest.csv`) that defines jobs and tracks their completion status.

### Rationale

The manifest serves multiple purposes:

1. **Single source of truth**: All metadata about a recording session (subject, date, probe type, channels, etc.) lives in one place.
2. **Status tracking**: Each processing stage has its own status column (`pending`, `complete`, `error`), enabling:
   - Selective re-processing of failed jobs
   - Resume-after-interruption capability
   - Clear visibility into pipeline state
3. **Batch processing**: The `run_preprocessing.m` script iterates through all jobs, processing only those with `pending` status.
4. **Reproducibility**: The manifest documents exactly which parameters were used for each session.

### Status Column Semantics

| Status | Meaning |
|--------|---------|
| `pending` | Ready to be processed |
| `complete` | Successfully processed |
| `error` | Processing failed; requires investigation |

### Multi-Probe Sessions

The `session_group_id` column links multiple `unique_id` entries that belong to the same physical recording session (e.g., simultaneous SC and SNc recordings). Each probe gets its own row because:

- Different probes may have different channel orderings
- Kilosort processes each probe's data independently
- This maintains a clean one-to-one mapping between manifest rows and output files

---

## Channel Ordering and Probe Validation

### Decision

Channel ordering for the `.dat` file is determined by the `probe_type` field in the manifest. An automated validation step checks whether the data's correlation structure matches the expected probe geometry, but the manifest is always the authoritative source.

### Rationale

Different probe types have different physical channel layouts. For Kilosort to correctly estimate spike positions along the probe, channels must be written to the `.dat` file in an order that reflects their physical arrangement.

The validation step (computing inter-channel correlations and checking against known probe templates) serves as a safety check. Adjacent channels on a probe should show higher correlation due to shared local signals. If the best-matching template differs from the manifest, a warning is issued—but the manifest value is still used.

### Why Manifest Over Auto-Detection?

- **Explicitness**: The user explicitly declares the probe type, reducing silent errors.
- **Edge cases**: Correlation-based detection can fail with unusual data (e.g., very low signal, artifacts).
- **Auditability**: The manifest provides a clear record of what was intended.

### Known Probe Orderings

The following probe types are currently defined in `prepare_spikes_for_kilosort.m`:

| Probe Type | Channel Ordering |
|------------|------------------|
| `vProbe` | `[32:-2:2, 31:-2:1]` |
| `nnVector` | `[17:2:31 18:2:32 2:2:16 1:2:15]` |

Adding new probe types requires updating the `knownOrderings` map in `prepare_spikes_for_kilosort.m`.

---

## Trial Alignment: The Anchor-and-Step Algorithm

### The Problem

Neural data (`.nev` file) and behavioral data (PLDAPS `.mat` files) are recorded by separate systems. To analyze neural activity in relation to behavioral events, trials must be aligned between these two data streams.

Challenges include:
- PLDAPS may be restarted mid-session, producing multiple `.mat` files
- Trial numbering may reset between files
- Minor differences in strobe sequences can occur due to timing edge cases

### The Solution: Anchor-and-Step

The alignment algorithm operates in two phases:

**Phase 1: Find an Anchor**
- Iterate through the first ~100 NEV trials
- For each, search for an exact strobe sequence match in the PLDAPS data
- Accept only unambiguous (single) matches
- If exact matching fails, fall back to Longest Common Subsequence (LCS) matching
- The first unambiguous match becomes the "anchor"

**Phase 2: Step Forward**
- From the anchor, iterate through subsequent NEV trials
- Use the `trialCount` strobe value (which increments monotonically within a PLDAPS session) for matching
- Search only within the current monotonic block (handle `trialCount` resets from PLDAPS restarts)
- This is faster and more robust than repeated LCS computation

### Rationale

The two-phase approach balances robustness with efficiency:
- Phase 1 handles the "cold start" problem where we have no prior alignment information
- Phase 2 exploits the structure of the data (monotonic trial counts) for fast subsequent matching
- LCS fallback handles cases where exact matching fails due to minor strobe differences

### Diagnostic Output

A figure comparing NEV and PDS `trialCount` values is saved to the `diagnostics/` folder. This allows visual verification that alignment succeeded.

---

## Timestamp Correction for Legacy Tasks

### The Problem

Two behavioral tasks (`gSac_4factors` and `tokens`) contain a bug in their PLDAPS implementation: stimulus onset strobes were not sent immediately after the screen flip but were instead queued and sent at the end of each trial.

This means the strobe timestamps in the `.nev` file do not reflect the true time of stimulus events for these tasks.

### The Solution

PLDAPS maintains its own record of event times (relative to `trialStartPTB`). The `trialEnd` strobe, which is sent via a different mechanism, has correct timing in both systems.

The correction procedure:
1. Identify trials with valid `trialEnd` times in both NEV and PDS
2. Fit a linear model: `Ripple_time = f(PLDAPS_time)`
3. Remove outliers (>3 standard deviations) and refit
4. Apply this mapping to convert PLDAPS-recorded event times to Ripple clock times

### Rationale

This approach works because:
- Clock drift between systems is approximately linear over a session
- The `trialEnd` strobe provides reliable anchor points throughout the session
- PLDAPS accurately records event times relative to trial start, just fails to strobe them at the right moment

### Scope and Future

This correction is applied only to `gSac_4factors` and `tokens` tasks. The underlying bug has been fixed in the PLDAPS task code, so future recordings do not require this correction.

The correction model parameters and diagnostic plots are saved for transparency and verification.

---

## Dynamic Field Discovery from PLDAPS

### Decision

Rather than hardcoding expected field names from PLDAPS data, the pipeline dynamically discovers all fields present in `p.trVars` and `p.trData` structures.

### Rationale

Different behavioral tasks save different variables. A task studying saccades will have fields like `targDegX` and `SRT`, while a task studying attention might have `cueLoc` and `chgLoc`. Hardcoding field names would require updating the pipeline for each new task.

Dynamic discovery means:
- New tasks work automatically without pipeline changes
- All available data is preserved
- The pipeline is robust to task evolution over time

### Implementation Details

- Scalar numeric fields become columns in the `trialInfo` table
- Non-scalar or non-numeric fields are stored in cell arrays
- Fields that are empty (`[]`) in the raw data become `NaN` in numeric arrays
- Field names are prefixed appropriately (e.g., `pds` prefix for PLDAPS timing fields) to distinguish data sources

### Multiple PLDAPS Files

When multiple PLDAPS files exist for a session (due to restarts):
1. Files are sorted chronologically by timestamp in filename (or modification date as fallback)
2. The union of all field names across files is computed
3. Data is merged into a single structure with proper trial indexing

---

## Waveform Extraction Strategy

### Decision

Mean waveforms are extracted in a separate step (`extract_waveforms`) before final consolidation, and are computed by averaging raw snippets from the `.dat` file.

### Rationale

**Why a separate step?**
- Waveform extraction is computationally intensive (reading large `.dat` files)
- Separating it allows re-running consolidation without re-extracting waveforms
- The preliminary `session_data.mat` serves as a checkpoint

**Why extract from `.dat` rather than `.ns5`?**
- The `.dat` file has channels already reordered to match probe geometry
- This ensures waveform channel indices are consistent with Kilosort's cluster assignments

### Parameters

- Window size: `[-40, +41]` samples around each spike (configurable in `pipeline_config.m`)
- At 30 kHz, this is approximately ±1.3 ms
- Both mean and standard deviation are computed and stored

### Edge Case Handling

Spikes near the beginning or end of the recording may have incomplete windows. These are zero-padded to maintain consistent array dimensions.

---

## Error Handling Philosophy

### Decision

Error handlers include `keyboard` statements that pause execution for interactive debugging rather than logging errors and continuing.

### Rationale

This pipeline is designed for interactive, supervised use rather than unattended batch processing. When an error occurs:
- The researcher can immediately inspect the workspace
- Variables like `ME` (the exception), `job`, and `config` are available
- The cause can be diagnosed and potentially fixed on the spot
- Execution can resume with `dbcont` or abort with `dbquit`

### Trade-offs

This approach prioritizes:
- **Debuggability** over automation
- **Data integrity** over throughput (errors are not silently skipped)
- **Transparency** over convenience

For a pipeline processing dozens of sessions, this is appropriate. For hundreds or thousands, a logging-based approach might be preferable.

---

## Known Issues and Future Considerations

### `.dat` Filename Discrepancy

There is a potential inconsistency between how `prepare_spikes_for_kilosort.m` names the `.dat` file (using `unique_id`) and how `extract_waveforms.m` looks for it (using `raw_filename_base`). This has not caused failures in practice, possibly due to Kilosort/Phy file handling or manual intervention. If errors arise, this should be investigated and standardized.

**Recommendation**: Standardize on `unique_id` for all pipeline-generated files.

### Anti-Correlation Report

The `generate_anticorrelation_report` function was added to identify cluster pairs that might be candidates for merging (based on complementary firing patterns). In practice, the utility of this report has been uncertain.

**Recommendation**: Consider making this optional (via a config flag) or removing it if it does not inform curation decisions.

### Hardcoded Probe Orderings

New probe types require code changes to `prepare_spikes_for_kilosort.m`. 

**Recommendation**: Consider moving probe definitions to a configuration file or the manifest itself.

### Figure Cleanup

Multiple `close all force` calls are scattered through the code, suggesting figures were being left open. This works but is somewhat ad-hoc.

**Recommendation**: Consider more structured figure management (e.g., explicit handles, cleanup in `finally` blocks).

---

## Document History

| Date | Author | Description |
|------|--------|-------------|
| 2026-01-05 | Claude (with James Herman) | Initial creation based on codebase review and Q&A |

---

## References

- `README.md` - Pipeline overview and usage instructions
- `session_data_dictionary.md` - Complete field documentation for output files
- `analysis_info.md` - Guide for downstream analysis
- Task-specific data dictionaries in `docs/task_data_dictionaries/`
