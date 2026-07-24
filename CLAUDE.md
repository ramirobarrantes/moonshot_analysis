# moonshot_pipeline_analysis

## Overview

Single-sample somatic variant analysis pipeline written in R. Intended to become a module within a larger Nextflow pipeline. Takes variant call files from one or more callers for a tumor/normal pair, annotates them, and produces RDS outputs for downstream use.

## Entry Point

`analysis.R` — run via `Rscript`:

```bash
Rscript analysis.R \
  <samplesheet.csv> \
  <sequenceDB> \
  <uniprotDBFile> \
  <cosmicDB> \
  <patient_id> \
  <output_dir>
```

### Positional arguments

| # | Name | Description |
|---|------|-------------|
| 1 | `samplesheet` | CSV with columns: `variant_file`, `name`, `caller`, `normal_name`, `tumor_name` |
| 2 | `sequenceDB` | Reference sequence database path |
| 3 | `uniprotDBFile` | UniProt annotation database file |
| 4 | `cosmicDB` | COSMIC recurrence database |
| 5 | `patient` | Patient/sample identifier (used in output filenames) |
| 6 | `output_dir` | Directory for output files (created if absent) |

## Outputs

Both written to `<output_dir>/`:

- `differenceMatrix_<patient>.RDS` — pairwise variant set differences across callers
- `allVariants_<patient>.RDS` — annotated variant table with peptide sequences, UniProt data, and COSMIC recurrence counts

## Key Helper Functions (from BSRBioinformaticsTools)

Sourced at runtime from a sibling directory — **not bundled here**. Do not read or modify files outside this directory.

| Function | Purpose |
|----------|---------|
| `createDifferenceMatrix()` | Compares variant calls across callers |
| `processVariantFiles()` | Parses VCF/variant files, applies caller-specific logic, annotates with sequence and UniProt data |
| `addReferenceAndMutantPeptide()` | Derives peptide sequences around each variant |
| `annotate_cosmic_recurrence()` | Looks up COSMIC hit counts via `Existing_variation` field |

## Pipeline Integration Notes

- This script is designed to run **per patient** — one invocation per sample.
- In the Nextflow context, each sample in the samplesheet becomes a channel item; `patient` maps to the sample ID.
- Sequence and protein fields are truncated to 30,000 characters before serialization to avoid oversized RDS files.
- The `toolsDirectory` path at line 3 is currently hardcoded — this must be parameterized or bundled before the script can run portably in Nextflow.

## Known Issues / TODOs

- Hardcoded `toolsDirectory` path needs to be replaced with an env var or Nextflow `params` input.
- No input validation — missing files will produce cryptic R errors rather than clear messages.
- No logging beyond the final `message("Done: ", patient)`.
