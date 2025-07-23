# --- Load libraries ---
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(TwoSampleMR)
library(stringr)

# --- Define directories ---
projectRoot <- '/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_190_viking_somalogic/'
analysisDir <- file.path(projectRoot, 'Scripts/twoSampleMr/somamerSNPlists/singleSNP/')
harmonisedOutDir <- file.path(projectRoot, 'Scripts/twoSampleMr/2resultsMR/harmonisedPreMRData/')
mrResultOutDir <- file.path(projectRoot, 'Scripts/twoSampleMr/2resultsMR/')

# --- Load available outcome datasets ---
ao <- available_outcomes()

# Optionally filter to EBI, IEU, UKB-A/B outcomes
ao2 <- ao %>% filter(str_detect(id, "ebi|ieu|ukb-b"))
ao2 <- ao %>% filter(str_detect(id, "ukb-a"))  # Only one active at a time

# --- Run MR for each file ---
for (fileName in list.files(analysisDir)) {
  exposureFilePath <- file.path(analysisDir, fileName)
  exposureName <- fileName

  exposure <- fread(exposureFilePath)

  # Format exposure data
  exp_df <- format_data(
    dat = data.table(exposure, phenotype = exposureName, id = exposureName),
    snp_col = "rsid",
    beta_col = "beta1",
    se_col = "se",
    effect_allele_col = "a1",
    other_allele_col = "a0",
    eaf_col = "freq1",
    pval_col = "p",
    samplesize_col = "n",
    phenotype_col = "phenotype",
    id_col = "id",
    type = "exposure"
  )

  # Extract outcomes for exposure SNPs
  outcome <- extract_outcome_data(snps = exp_df$SNP, outcomes = ao2$id)

  # Harmonise data
  harmonisedData <- harmonise_data(exposure_dat = exp_df, outcome_dat = outcome)
  harmonisedData <- harmonisedData %>% distinct()

  # Export harmonised data
  harmonisedFile <- file.path(harmonisedOutDir, paste0("harmonised_data_", exposureName, ".tsv"))
  fwrite(harmonisedData, harmonisedFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  # Load harmonised data back (you had hardcoded one file earlier — this generalises it)
  harmonisedData <- fread(harmonisedFile)

  # Perform Mendelian Randomisation
  res <- mr(harmonisedData, method_list = c("mr_wald_ratio", "mr_ivw_fe"))

  # Format results
  res1 <- data.table(res)[, .(
    exposure,
    outcome,
    method = gsub(".", "_", tolower(make.names(method)), fixed = TRUE),
    nsnp,
    beta = b,
    se,
    p = pval
  )]

  # Export MR results
  resultFile <- file.path(mrResultOutDir, paste0("MRresult_", exposureName))
  fwrite(res1, resultFile, sep = ",", quote = FALSE, na = "NA")
}
