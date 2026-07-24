#!/usr/bin/env bash
# Generates all derived samplesheets from a sourced run conf.
# Usage: source conf/run_<PATIENT>.conf && ./generate_samplesheets.sh
set -euo pipefail

: "${PATIENT:?PATIENT not set — source a run conf first}"

mkdir -p "$BASE_OUTDIR"

# ── qcSamplesheet.csv (seqinspector input) ────────────────────────────────────
QC_SAMPLESHEET="${BASE_OUTDIR}/qcSamplesheet.csv"
{
  echo "sample,fastq_1,fastq_2,rundir,tags"
  echo "${NORMAL_SAMPLE},${NORMAL_FASTQ_1},${NORMAL_FASTQ_2},,paired_sample"
  for i in "${!TUMOR_SAMPLES[@]}"; do
    echo "${TUMOR_SAMPLES[$i]},${TUMOR_FASTQ_1[$i]},${TUMOR_FASTQ_2[$i]},,paired_sample"
  done
} > "$QC_SAMPLESHEET"
echo "Wrote $QC_SAMPLESHEET"

# ── sarek samplesheet.csv ──────────────────────────────────────────────────────
SAREK_SAMPLESHEET="${BASE_OUTDIR}/samplesheet.csv"
{
  echo "patient,status,sample,lane,fastq_1,fastq_2"
  echo "${PATIENT},0,${NORMAL_SAMPLE},1,${NORMAL_FASTQ_1},${NORMAL_FASTQ_2}"
  for i in "${!TUMOR_SAMPLES[@]}"; do
    echo "${PATIENT},1,${TUMOR_SAMPLES[$i]},1,${TUMOR_FASTQ_1[$i]},${TUMOR_FASTQ_2[$i]}"
  done
} > "$SAREK_SAMPLESHEET"
echo "Wrote $SAREK_SAMPLESHEET"

# ── moonshot_pipeline sampleSheetMoonshot.csv (CRAMs from sarek recalibration) ─
MOONSHOT_SAMPLESHEET="${BASE_OUTDIR}/sampleSheetMoonshot.csv"
NORMAL_CRAM="${SAREK_OUTDIR}/preprocessing/recalibrated/${NORMAL_SAMPLE}/${NORMAL_SAMPLE}.recal.cram"
{
  echo "patient,tumor_cram,tumor_crai,normal_cram,normal_crai,tumor_sample,normal_sample"
  for tumor in "${TUMOR_SAMPLES[@]}"; do
    tumor_cram="${SAREK_OUTDIR}/preprocessing/recalibrated/${tumor}/${tumor}.recal.cram"
    echo "${tumor},${tumor_cram},${tumor_cram}.crai,${NORMAL_CRAM},${NORMAL_CRAM}.crai,${PATIENT}_${tumor},${PATIENT}_${NORMAL_SAMPLE}"
  done
} > "$MOONSHOT_SAMPLESHEET"
echo "Wrote $MOONSHOT_SAMPLESHEET"

# ── final analysis samplesheet.csv (variant_file,name,caller,normal_name,tumor_name)
# Paths follow sarek / moonshot_pipeline's deterministic output naming.
ANALYSIS_SAMPLESHEET="${BASE_OUTDIR}/analysisSamplesheet.csv"
{
  echo "variant_file,name,caller,normal_name,tumor_name"
  for tumor in "${TUMOR_SAMPLES[@]}"; do
    pair_name="${tumor}_vs_${NORMAL_SAMPLE}"

    mutect_vcf="${SAREK_OUTDIR}/annotation/mutect2/${pair_name}/${pair_name}.mutect2.filtered_VEP.ann.vcf.gz"
    echo "${mutect_vcf},${PATIENT}_mutect_${tumor},mutect,${PATIENT}_${NORMAL_SAMPLE},${PATIENT}_${tumor}"

    strelka_vcf="${SAREK_OUTDIR}/annotation/strelka/${pair_name}/${pair_name}.strelka.somatic_snvs_VEP.ann.vcf.gz"
    echo "${strelka_vcf},${PATIENT}_strelka_${tumor},strelka,${NORMAL_SAMPLE},${tumor}"

    muse_vcf="${MOONSHOT_PIPELINE_OUTDIR}/muse/annotation/${tumor}.vep.vcf.gz"
    echo "${muse_vcf},${PATIENT}_muse2_${tumor},muse2,${PATIENT}_${NORMAL_SAMPLE},${PATIENT}_${tumor}"

    tnscope_vcf="${MOONSHOT_PIPELINE_OUTDIR}/tnscope/annotation/${tumor}.vep.vcf.gz"
    echo "${tnscope_vcf},${PATIENT}_tnscope_${tumor},tnscope,${PATIENT}_${NORMAL_SAMPLE},${PATIENT}_${tumor}"
  done
} > "$ANALYSIS_SAMPLESHEET"
echo "Wrote $ANALYSIS_SAMPLESHEET"
