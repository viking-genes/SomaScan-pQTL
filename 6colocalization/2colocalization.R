# Refactored version of Arianna's s18_cisPQTLs_metabolomics_FIRST_COLOC.R

# Load packages
Add_Directory_To_Rlibs <- function(directory) {
  .libPaths(unique(c(directory, .libPaths())))
}
Add_Directory_To_Rlibs("./R/3.6.3")

library(data.table)
library(dplyr)
library(tidyr)
library(coloc)
library(ieugwasr)

# Manueal sample size assignment where it's missing in the database
sampleSizeOverrides <- list(
  "ieu-b-4870" = 214989, "ieu-b-4871" = 185221, "ieu-b-111" = 441016,
  "ieu-b-4808" = 441291, "ieu-b-4812" = 441291, "ieu-b-4865" = 199569,
  "ieu-b-102" = 500199, "ieu-b-4861" = 244207, "ieu-b-4862" = 205527,
  "ieu-b-4864" = 199569, "ieu-b-107" = 393193, "ieu-b-108" = 439214,
  "ieu-b-110" = 440546, "ieu-b-109" = 403943, "ieu-b-104" = 46368,
  "ieu-b-105" = 353315, "ieu-b-4869" = 180386, "ebi-a-GCST005647" = 80610,
  "ieu-b-72" = 326938
)

# Incorrect assignment of base pair positions. Skipping for consistency
skipDatasets <- c("ebi-a-GCST005071", "ebi-a-GCST005647")

# Helper function to fix sample size issues
fixSampleSizes <- function(gwasData, id) {
  if (id %in% names(sampleSizeOverrides)) {
    gwasData$N <- sampleSizeOverrides[[id]]
  }
  # ieu-b-4810 has each row duplicated. This removes duplicate rows
  if (id == "ieu-b-4810") {
    gwasData <- gwasData[!duplicated(gwasData), ]
  }
  return(gwasData)
}

# Load and filter protein file
loadProteinFile <- function(filepath, chr, start, end, exposure) {
  fread(filepath, sep = "\t", data.table = FALSE,
        select = c("rsid", "chr", "pos", "a1", "a0", "freq1", "beta1", "se", "p", "n")) %>%
    filter(chr == chr, pos >= start, pos <= end) %>%
    mutate(maf = ifelse(freq1 < 0.5, freq1, 1 - freq1),
           id = exposure) %>%
    rename(beta = beta1, pvalue = p, ea = a1, oa = a0, eaf = freq1)
}

# Load and filter IEU GWAS data
loadIEUGWAS <- function(chr, start, end, id) {
  variants <- paste0(chr, ":", start, "-", end)
  ieuData <- as.data.frame(associations(variants = variants, id = id, proxies = 0)) %>%
    mutate(maf = ifelse(eaf < 0.5, eaf, 1 - eaf)) %>%
    rename(N = n, pvalue = p, pos = position, oa = nea) %>%
    select(rsid, chr, pos, ea, oa, eaf, beta, se, pvalue, N, maf, trait)

  fixSampleSizes(ieuData, id)
}

# Merge GWAS datasets by rsID
mergeDatasets <- function(proteinData, metaanalysisData) {
  proteinData <- proteinData %>%
    rename_at(vars(-rsid, -chr, -pos, -id), ~ paste0(., "Protein"))
  metaanalysisData <- metaanalysisData %>%
    rename_at(vars(-rsid, -chr, -pos, -trait), ~ paste0(., "Metaanalysis"))
  
  merge(proteinData, metaanalysisData, by = c("rsid", "chr", "pos"))
}

# Run coloc analysis
runColoc <- function(mergedData) {
  coloc.abf(
    dataset1 = list(
      beta = mergedData$betaProtein,
      varbeta = mergedData$seProtein^2,
      N = mean(mergedData$nProtein),
      type = "quant",
      snp = mergedData$rsid
    ),
    dataset2 = list(
      beta = mergedData$betaMetaanalysis,
      varbeta = mergedData$seMetaanalysis^2,
      N = mean(mergedData$NMetaanalysis),
      type = "quant",
      snp = mergedData$rsid
    ),
    MAF = mergedData$mafProtein
  )
}

# Main pipeline
runPipeline <- function(inputFilePath, resultsDir, gwasBasePath) {
  inputDf <- fread(inputFilePath)
  
  for (i in seq_len(nrow(inputDf))) {
    exposure <- inputDf[i, "exposure"]
    somamerID <- inputDf[i, "somamerID"]
    chr <- inputDf[i, "chr"]
    bpStart <- inputDf[i, "bpStart"]
    bpEnd <- inputDf[i, "bpEnd"]
    metaanalysisID <- tail(strsplit(inputDf[i, "outcome"], ":")[[1]], 1)

    if (metaanalysisID %in% skipDatasets) {
      next
    }

    # Paths and file loading
    gwasFile <- file.path(gwasBasePath, paste0(somamerID, ".tsv.gz"))
    proteinDf <- loadProteinFile(gwasFile, chr, bpStart, bpEnd, exposure)
    
    # Replace missing rsIDs
    proteinDf$rsid[proteinDf$rsid == "."] <- paste0(proteinDf$chr[proteinDf$rsid == "."], ":", proteinDf$pos[proteinDf$rsid == "."])
    
    ieuDf <- loadIEUGWAS(chr, bpStart, bpEnd, metaanalysisID)
    
    # Remove duplicated rsIDs
    proteinDf <- proteinDf[!duplicated(proteinDf$rsid) & !duplicated(proteinDf$rsid, fromLast = TRUE), ]
    ieuDf <- ieuDf[!duplicated(ieuDf$rsid) & !duplicated(ieuDf$rsid, fromLast = TRUE), ]
    
    # Merge and run coloc
    mergedDf <- mergeDatasets(proteinDf, ieuDf)
    colocResults <- runColoc(mergedDf)
    
    # Save results
    outputFile <- file.path(resultsDir, paste0(exposure, "_", metaanalysisID, ".tsv.gz"))
    write.table(colocResults$summary, outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# Run the pipeline
runPipeline(
  inputFilePath = "./prj_190_viking_somalogic/Scripts/colocalization/1mergedTwoSampleMRResults.tsv",
  resultsDir = "./prj_190_viking_somalogic/Scripts/colocalization/results",
  gwasBasePath = "./prj_190_viking_somalogic/GWAS/Full/p05_comb_chr"
)
