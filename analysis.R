#!/usr/bin/env Rscript

library(BSRBioinformaticsTools)
library(VariantAnnotation)
library(SummarizedExperiment)
source("/opt/helperFunctions.R")

args <- commandArgs(trailingOnly = TRUE)

samplesheet   <- args[1]
sequenceDB    <- args[2]
uniprotDBFile <- args[3]
cosmicDB      <- args[4]
patient       <- args[5]
output_dir    <- args[6]

# ── Load samplesheet ──────────────────────────────────────────────────────────
samples             <- read.csv(samplesheet, stringsAsFactors = FALSE)
variantFiles        <- samples$variant_file
variantListNames    <- samples$name
variantListCallers  <- samples$caller
normalNames         <- samples$normal_name
tumorNames          <- samples$tumor_name

# ── Create output directory ───────────────────────────────────────────────────
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Analysis ──────────────────────────────────────────────────────────────────
differenceMatrix <- createDifferenceMatrix(
  variantFiles, variantListNames, variantListCallers, normalNames, tumorNames
)
saveRDS(differenceMatrix,
        file = file.path(output_dir, paste0("differenceMatrix_", patient, ".RDS")))

allVariants <- processVariantFiles(
  variantFiles, variantListCallers, sequenceDB,
  uniprotDBFile = uniprotDBFile, normalNames, tumorNames
)

allVariantsWithPeptide <- addReferenceAndMutantPeptide(allVariants)

cosmic_db <- read.csv(cosmicDB, sep = "\t", check.names = FALSE)
cosmic_hits <- annotate_cosmic_recurrence(
  allVariantsWithPeptide, "Existing_variation", cosmic_db
)
allVariantsWithPeptide[, c("cosmic_sample_count",
                           "cosmic_study_count",
                           "confirmed_somatic_n")] <- cosmic_hits

allVariantsWithPeptide$sequence        <- substr(allVariantsWithPeptide$sequence, 0, 30000)
allVariantsWithPeptide$proteinSequence <- substr(allVariantsWithPeptide$proteinSequence, 0, 30000)

saveRDS(allVariantsWithPeptide,
        file = file.path(output_dir, paste0("allVariants_", patient, ".RDS")))

# ── Extract VAFs ──────────────────────────────────────────────────────────────
vaf_list <- lapply(seq_len(nrow(samples)), function(i) {
  vcf <- VariantAnnotation::readVcf(samples$variant_file[i], genome = "hg38")
  vcf <- vcf[rowRanges(vcf)$FILTER == "PASS"]
  tumor_col <- samples$tumor_name[i]
  g <- VariantAnnotation::geno(vcf)

  if (!is.null(g$AD) && tumor_col %in% colnames(g$AD)) {
    vaf <- sapply(seq_len(nrow(vcf)), function(j) {
      counts <- g$AD[j, tumor_col][[1]]
      if (length(counts) >= 2 && sum(counts) > 0) counts[2] / sum(counts) else NA
    })
  } else if (!is.null(g$DP) && tumor_col %in% colnames(g$DP)) {
    # Strelka: derive alt count from nucleotide-specific fields (AU/CU/GU/TU)
    alt_alleles <- as.character(sapply(alt(vcf), function(x) as.character(x[1])))
    dp <- as.numeric(g$DP[, tumor_col])
    vaf <- sapply(seq_len(nrow(vcf)), function(j) {
      field <- paste0(alt_alleles[j], "U")
      if (field %in% names(g) && tumor_col %in% colnames(g[[field]])) {
        alt_count <- g[[field]][j, tumor_col, 1]
        if (!is.na(dp[j]) && dp[j] > 0) alt_count / dp[j] else NA
      } else NA
    })
  } else {
    return(NULL)
  }
  data.frame(caller = samples$caller[i], vaf = vaf)
})
vaf_df <- do.call(rbind, Filter(Negate(is.null), vaf_list))
saveRDS(vaf_df, file = file.path(output_dir, paste0("VAF_", patient, ".RDS")))

message("Done: ", patient)

# ── Generate report ───────────────────────────────────────────────────────────
qmd_file <- file.path(output_dir, paste0(patient, ".qmd"))
html_file <- paste0(patient, ".html")
file.copy("/opt/template.qmd", qmd_file, overwrite = TRUE)
quarto::quarto_render(
  qmd_file,
  execute_params = list(patient = patient, data_dir = dirname(output_dir)),
  output_file = html_file
)
message("Report: ", file.path(output_dir, html_file))
