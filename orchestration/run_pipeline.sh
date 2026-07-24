#!/usr/bin/env bash
# Full moonshot pipeline orchestration for a single patient.
#
# Usage: ./run_pipeline.sh conf/run_MS008.conf [--skip-db] [--skip-seqinspector] [--skip-sarek] [--skip-moonshot-pipeline] [--skip-analysis]
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <conf-file> [--skip-db] [--skip-seqinspector] [--skip-sarek] [--skip-moonshot-pipeline] [--skip-analysis]"
  exit 1
fi

CONF_FILE="$1"; shift
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SKIP_DB=false
SKIP_SEQINSPECTOR=false
SKIP_SAREK=false
SKIP_MOONSHOT_PIPELINE=false
SKIP_ANALYSIS=false

for arg in "$@"; do
  case "$arg" in
    --skip-db) SKIP_DB=true ;;
    --skip-seqinspector) SKIP_SEQINSPECTOR=true ;;
    --skip-sarek) SKIP_SAREK=true ;;
    --skip-moonshot-pipeline) SKIP_MOONSHOT_PIPELINE=true ;;
    --skip-analysis) SKIP_ANALYSIS=true ;;
    *) echo "Unknown flag: $arg"; exit 1 ;;
  esac
done

# shellcheck disable=SC1090
source "$CONF_FILE"
: "${PATIENT:?PATIENT not set in $CONF_FILE}"

LOG_DIR="${BASE_OUTDIR}/orchestration_logs"
mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Step 0: reference databases (infrequent — skip by default once populated) ─
if [ "$SKIP_DB" = false ]; then
  log "Running moonshot_pipeline_db (ensembl_version=${ENSEMBL_VERSION})"
  ( cd "$MOONSHOT_DB_REPO" && nextflow run main.nf \
      -profile slurm \
      --ensembl_version "$ENSEMBL_VERSION" \
      --outdir "$DB_OUTDIR" \
  ) 2>&1 | tee "${LOG_DIR}/00_db.log"
else
  log "Skipping DB pipeline"
fi

# ── Generate samplesheets ──────────────────────────────────────────────────────
log "Generating samplesheets"
"${SCRIPT_DIR}/generate_samplesheets.sh"

# ── Step 1: seqinspector + sarek in parallel (both only need raw FASTQs) ──────
PIDS=()

if [ "$SKIP_SEQINSPECTOR" = false ]; then
  log "Launching seqinspector"
  (
    nextflow run "${SEQINSPECTOR_REPO}/main.nf" \
      -profile apptainer \
      -r "$SEQINSPECTOR_VERSION" \
      --input "${BASE_OUTDIR}/qcSamplesheet.csv" \
      --fastq_screen_references "$FASTQ_SCREEN_REFERENCES" \
      --fasta "$FASTA" \
      --outdir "$SEQINSPECTOR_OUTDIR" \
      --work "$WORKDIR"
  ) > "${LOG_DIR}/01_seqinspector.log" 2>&1 &
  PIDS+=($!)
else
  log "Skipping seqinspector"
fi

if [ "$SKIP_SAREK" = false ]; then
  log "Launching sarek"
  (
    export NXF_APPTAINER_CACHEDIR="${NXF_APPTAINER_CACHEDIR:-/gpfs1/mbsr_tools/NXF_APPTAINER_CACHEDIR}"
    NXF_VER=25.10.4 nextflow run nf-core/sarek \
      -r "$SAREK_VERSION" \
      -profile apptainer \
      --input "${BASE_OUTDIR}/samplesheet.csv" \
      --outdir "$SAREK_OUTDIR" \
      --tools "$SAREK_TOOLS" \
      --genome "$SAREK_GENOME" \
      --wes \
      --work "$WORKDIR"
  ) > "${LOG_DIR}/02_sarek.log" 2>&1 &
  PIDS+=($!)
else
  log "Skipping sarek"
fi

if [ ${#PIDS[@]} -gt 0 ]; then
  log "Waiting for seqinspector/sarek to finish (PIDs: ${PIDS[*]})"
  FAIL=0
  for pid in "${PIDS[@]}"; do
    wait "$pid" || FAIL=1
  done
  if [ "$FAIL" -ne 0 ]; then
    log "seqinspector or sarek failed — check ${LOG_DIR}/01_seqinspector.log and ${LOG_DIR}/02_sarek.log"
    exit 1
  fi
  log "seqinspector and sarek completed"
fi

# ── Step 2: moonshot_pipeline (muse, tnscope, ichorCNA, purple) — needs sarek CRAMs
if [ "$SKIP_MOONSHOT_PIPELINE" = false ]; then
  log "Running moonshot_pipeline"
  ( cd "$MOONSHOT_PIPELINE_REPO" && nextflow run main.nf \
      --input "${BASE_OUTDIR}/sampleSheetMoonshot.csv" \
      --fasta "$FASTA" \
      --fai "$FAI" \
      --dbsnp "$DBSNP" \
      --dbsnp_tbi "$DBSNP_TBI" \
      --gc_wig "$GC_WIG" \
      --map_wig "$MAP_WIG" \
      --wgs false \
      --vep_cache "$VEP_CACHE" \
      --vep_cache_version "$VEP_CACHE_VERSION" \
      --germline_vcf "$GERMLINE_VCF" \
      --germline_vcf_tbi "$GERMLINE_VCF_TBI" \
      --contamination_vcf "$CONTAMINATION_VCF" \
      --contamination_vcf_tbi "$CONTAMINATION_VCF_TBI" \
      --het_pon "$HET_PON" \
      --het_pon_tbi "$HET_PON_TBI" \
      --gc_profile "$GC_PROFILE" \
      --outdir "$MOONSHOT_PIPELINE_OUTDIR" \
      --dict "$DICT" \
      --ensembl_data_dir "$ENSEMBL_DATA_DIR" \
      --work "$WORKDIR" \
  ) 2>&1 | tee "${LOG_DIR}/03_moonshot_pipeline.log"
else
  log "Skipping moonshot_pipeline"
fi

# ── Step 3: final analysis + report (this repo, dockerized) ────────────────────
if [ "$SKIP_ANALYSIS" = false ]; then
  log "Running moonshot_analysis"
  mkdir -p "$ANALYSIS_OUTDIR"
  docker run --rm \
    -e ANTHROPIC_API_KEY="${ANTHROPIC_API_KEY:-}" \
    -v "${BASE_OUTDIR}:/data" \
    -v "${DB_OUTDIR}:/dbs" \
    "$ANALYSIS_DOCKER_IMAGE" \
    "/data/analysisSamplesheet.csv" \
    "/dbs/$(realpath --relative-to="$DB_OUTDIR" "$ENSEMBL_PEPTIDES_DB")" \
    "/dbs/$(realpath --relative-to="$DB_OUTDIR" "$UNIPROT_DB")" \
    "/dbs/$(realpath --relative-to="$DB_OUTDIR" "$COSMIC_DB")" \
    "$PATIENT" \
    "/data/$(basename "$ANALYSIS_OUTDIR")" \
    2>&1 | tee "${LOG_DIR}/04_analysis.log"
else
  log "Skipping analysis"
fi

log "Pipeline complete for $PATIENT"
