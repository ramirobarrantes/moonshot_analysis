# moonshot_pipeline_analysis

Single-sample somatic variant analysis and reporting for the Moonshot project. This repo holds the **final analysis stage** — variant annotation, cross-caller comparison, and an HTML report — plus an orchestration layer that ties together the upstream pipelines needed to produce its inputs.

## Pipeline overview

A full run for one patient involves four pipelines, run in this order:

```
moonshot_pipeline_db (infrequent — refreshes reference DBs)
       │
       ▼ (ensembl_peptides, uniprot_database)
       │
seqinspector ──────────────┐   both consume raw FASTQs only,
sarek (mutect2, strelka)   │   so they run in parallel
       └───────┬───────────┘
               ▼ (recalibrated CRAMs from sarek)
       moonshot_pipeline (muse, tnscope, ichorCNA, purple; also
                          produces BAMs from sarek's CRAMs)
               ▼ (VCFs from sarek + moonshot_pipeline, DB outputs)
       moonshot_analysis  (this repo)
```

| Stage | Repo | Produces |
|---|---|---|
| DB pipeline | `moonshot_pipeline_db` | `ensembl_peptides/`, `uniprot_database/` (COSMIC DB is added manually, not generated) |
| seqinspector | `seqinspector` | Raw-read QC (FastQC, FastQ Screen, MultiQC) |
| sarek | `nf-core/sarek` | Alignment, mutect2 + strelka calls, VEP annotation, pipeline QC |
| moonshot_pipeline | `moonshot_pipeline` | muse2 + tnscope calls, ichorCNA tumor fraction, PURPLE purity/ploidy, BAMs |
| **moonshot_analysis** | **this repo** | Cross-caller variant comparison, annotated variant table, subcellular location summary, tumor purity, VAF plots, QC summary, HTML report |

## This repo

### Entry point

`analysis.R`, run via Rscript:

```bash
Rscript analysis.R \
  <samplesheet.csv> \
  <sequenceDB> \
  <uniprotDBFile> \
  <cosmicDB> \
  <patient_id> \
  <output_dir>
```

See [CLAUDE.md](CLAUDE.md) for argument details, helper function reference, and known issues.

### Docker

```bash
docker build -t moonshot-analysis:latest .

docker run --rm \
  -e ANTHROPIC_API_KEY="$ANTHROPIC_API_KEY" \
  -v /path/to/patient/output:/data \
  -v /path/to/DBs:/dbs \
  moonshot-analysis:latest \
  /data/analysisSamplesheet.csv \
  /dbs/ensembl_peptides/ensembl_peptides.csv \
  /dbs/uniprot_database/HUMAN_9606_Database.RDS \
  /dbs/Cosmic_MutantCensus_Tsv_v104_GRCh37/Cosmic_MutantCensus_v104_GRCh37.tsv \
  MS008 \
  /data/finalAnalysis
```

`ANTHROPIC_API_KEY` is optional — it enables the automated QC summary paragraph in the report (calls the Claude API on MultiQC's `llms-full.txt` output from both seqinspector and sarek). Without it, the report shows a placeholder in that section.

### Report

`template.qmd` renders one HTML report per patient, expecting `data_dir` to be the patient's base output directory containing `seqinspector/`, `sarek_output/`, and `output/` (the analysis output) as siblings — this is the layout `orchestration/run_pipeline.sh` produces.

Report sections, top to bottom: QC summary (LLM-generated) → QC report links → download links (CSV exports) → variant caller comparison → subcellular location summary (overall and per-caller) → tumor purity → VAF distribution.

## Orchestration

`orchestration/` runs the full four-pipeline chain for one patient from a single config file.

1. Copy the example conf and fill in patient/sample info and reference paths:

   ```bash
   cp orchestration/conf/run_MS008.conf.example orchestration/conf/run_MS008.conf
   # edit run_MS008.conf
   ```

2. Run the full pipeline:

   ```bash
   ./orchestration/run_pipeline.sh orchestration/conf/run_MS008.conf
   ```

   `generate_samplesheets.sh` is called automatically and derives every downstream samplesheet (`qcSamplesheet.csv`, sarek's `samplesheet.csv`, `sampleSheetMoonshot.csv`, `analysisSamplesheet.csv`) from the sample list declared once in the conf — you don't hand-maintain each CSV.

3. Re-run a subset of stages with skip flags, e.g. after fixing a failed sarek run:

   ```bash
   ./orchestration/run_pipeline.sh orchestration/conf/run_MS008.conf --skip-db --skip-seqinspector
   ```

   Available flags: `--skip-db`, `--skip-seqinspector`, `--skip-sarek`, `--skip-moonshot-pipeline`, `--skip-analysis`.

Logs for each stage are written to `<BASE_OUTDIR>/orchestration_logs/`.

### Config file

One conf file per patient (`orchestration/conf/run_<PATIENT>.conf`, gitignored — copy from the `.example`). It declares:

- Patient ID, normal sample + FASTQs, and a list of tumor samples + FASTQs (supports multiple tumor samples per patient)
- Output/work directories for every stage
- Paths to the four upstream pipeline repos
- Reference genome, VEP cache, and PURPLE/ichorCNA resource paths
- DB pipeline output locations (feeds moonshot_analysis directly)
- Tool versions (sarek, seqinspector) and the Docker image tag for this repo

### Assumed directory layouts

The samplesheet generator assumes:

- sarek: `<SAREK_OUTDIR>/preprocessing/recalibrated/<SAMPLE>/<SAMPLE>.recal.cram`, `<SAREK_OUTDIR>/annotation/{mutect2,strelka}/<TUMOR>_vs_<NORMAL>/...vcf.gz`
- moonshot_pipeline: `<MOONSHOT_PIPELINE_OUTDIR>/{muse,tnscope}/annotation/<TUMOR>.vep.vcf.gz`

If either pipeline's output naming changes, update `orchestration/generate_samplesheets.sh` accordingly.

## Requirements

- Docker (for the analysis stage)
- Nextflow + a configured executor profile (`slurm`, `apptainer`) for the upstream pipelines
- An Anthropic API key (optional, for the automated QC summary)
