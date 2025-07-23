#!/usr/bin/env Rscript

# --- Load libraries ---
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(stringr)

# --- Setup directories ---
projectRoot <- "/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_190_viking_somalogic/"
inputDir <- file.path(projectRoot, "Scripts/twoSampleMr/3significantResultsMR/Split2/")
gwasDir <- file.path(projectRoot, "GWAS/Full/p05_comb_chr/")
outputDir <- file.path(projectRoot, "Scripts/twoSampleMr/4resultsReverseMR/")

# --- Load available outcomes once ---
ao <- available_outcomes()

# --- Define genomic regions to exclude (ABO and HLA) ---
exclude_regions <- list(
  list(chr = 9, start = 136131052, end = 136150605),   # ABO
  list(chr = 6, start = 29645000,  end = 33365000)     # HLA
)

# --- Loop through all protein files ---
for (proteinName in list.files(inputDir)) {
  proteinFilePath <- file.path(inputDir, proteinName)
  proteinDf <- fread(proteinFilePath)
  proteinFileLength <- nrow(proteinDf)

  # Derive file name for protein GWAS summary statistics
  splitName <- str_split(proteinName, '-')[[1]]
  if (length(splitName) < 2) next  # skip malformed names
  proteinGwasFilename <- paste0(splitName[1], "__", substr(splitName[2], 1, nchar(splitName[2]) - 4))
  proteinGwasFilePath <- file.path(gwasDir, paste0("str", proteinGwasFilename, ".tsv.gz"))

  # --- Loop through all outcome rows (exposures) in the file ---
  for (i in 1:proteinFileLength) {
    outcomeValue <- proteinDf[i, 2][[1]]
    outcomeId <- tail(str_split(outcomeValue, ':')[[1]], n = 1)

    # Try to extract instruments
    exp_df <- extract_instruments(outcomes = outcomeId)
    if (is.null(nrow(exp_df))) {
      cat("NOT FOUND:", outcomeValue, "\n")
      next
    }

    # Convert chromosome and position to numeric
    exp_df$chr.exposure <- as.numeric(exp_df$chr.exposure)
    exp_df$pos.exposure <- as.numeric(exp_df$pos.exposure)

    # Exclude SNPs in ABO and HLA regions
    for (region in exclude_regions) {
      exp_df <- exp_df[!(exp_df$chr.exposure == region$chr &
                         exp_df$pos.exposure > region$start &
                         exp_df$pos.exposure < region$end), ]
    }

    # Read outcome (protein) GWAS data
    out_df <- read_outcome_data(
      snps = exp_df$SNP,
      filename = proteinGwasFilePath,
      sep = "\t",
      snp_col = "rsid",
      beta_col = "beta1",
      se_col = "se",
      effect_allele_col = "a1",
      other_allele_col = "a0",
      eaf_col = "freq1",
      pval_col = "p",
      samplesize_col = "n"
    )

    out_df$outcome <- proteinName

    # Harmonise exposure and outcome
    harmonised <- harmonise_data(exposure_dat = exp_df, outcome_dat = out_df)

    # Run MR
    res <- mr(harmonised, method_list = c(
      "mr_wald_ratio",
      "mr_ivw",
      "mr_ivw_fe",
      "mr_weighted_mode",
      "mr_weighted_median",
      "mr_egger_regression",
      "mr_two_sample_ml"
    ))

    # Format and save results
    resFormatted <- data.table(res)[, .(
      exposure,
      outcome,
      method = gsub(".", "_", tolower(make.names(method)), fixed = TRUE),
      nsnp,
      beta = b,
      se,
      p = pval
    )]

    outputFile <- file.path(outputDir, paste0(proteinName, "_", i))
    fwrite(resFormatted, outputFile, sep = ",", quote = FALSE, na = "NA")
  }
}
