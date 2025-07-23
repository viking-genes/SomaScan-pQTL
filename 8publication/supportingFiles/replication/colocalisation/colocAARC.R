#module load igmm/apps/R/3.6.3

# --- Load Required Libraries ---
library(data.table)
library(coloc)
library(ieugwasr)
library(tidyr)
library(dplyr)
library(biomaRt)

# --- Define Project Directories ---
projectRoot <- "./prj_190_viking_somalogic"

gwasVikingDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/GWASCatalogUpload/convertedGWASFiles")
gwasAaskDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/Surapaneni2022Data/fullStats")
supplementaryTable1Path <- file.path(projectRoot, "Scripts/Publication/supplementary1.xlsx")


#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwasAASKdir and perform colocalisation
for (filename in list.files(gwasAASKdir)) {

  #Print progress
  print(paste("Processing", filename))

  #Extract the somamer ID from the file name
  somamerID <- gsub(".gz", "", filename)

  #Check if the file has already been processed
  resultFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", somamerID, ".tsv"))
  if (file.exists(resultFilepath)) {
    next
  }

  #Find the row corresponding to the somamer ID and retrieve the sentinel SNP
  somamerRow <- supplementary1[supplementary1$somamerID == somamerID,]
  sentinelSNP <- somamerRow$SNP

  # Load the GWAS summary statistics
  gwasVIKING <- fread(cmd = paste("gunzip -c", paste0(gwasVIKINGdir, "/", somamerID, ".tsv.gz")), sep = "\t", header = TRUE)

  # Window selection
  windowSize <- 300000

  SNPchromosome <- as.numeric(gwasVIKING[gwasVIKING$rs_id == sentinelSNP, "chromosome"])
  SNPposition <- as.numeric(gwasVIKING[gwasVIKING$rs_id == sentinelSNP, "base_pair_location"])

  # Construct VIKING table
  startPosition <- SNPposition - windowSize
  endPosition <- SNPposition + windowSize

  gwasVIKINGFiltered <- subset(gwasVIKING, chromosome == SNPchromosome & base_pair_location >= startPosition & base_pair_location <= endPosition)

  # Construct AASK table around the queried sentinel SNP
  convert_AASK_to_GRCh37 <- function(somamerID, sentinelSNP) {

    outputFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/AARC_GRCh37/", somamerID, ".tsv"))

    #Check if the file already exists
    if (file.exists(outputFilepath)) {
      return(fread(outputFilepath, sep = "\t", header = TRUE))
    }

    gwasAASKFull <- fread(cmd = paste("gunzip -c", paste0(gwasAASKdir, "/", somamerID, ".gz")), sep = "\t", header = TRUE)

    # Extract larger region for liftover to be reduced later, just in case the region changed between GRCh37 and GRCh38
    # Query ensembl for GRCh38 position of the sentinel SNP
    ensemblGrCh38 <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
    martOutputGrCh38 <- getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),
      filters = "snp_filter",
      values = sentinelSNP,
      mart = ensemblGrCh38,
      useCache = FALSE)

    SNPpositionGrCh38 <- as.numeric(martOutputGrCh38$chrom_start)
    startPositionGrCh38 <- SNPpositionGrCh38 - 3*windowSize
    endPositionGrCh38 <- SNPpositionGrCh38 + 3*windowSize

    gwasAASK <- subset(gwasAASKFull, chromosome == SNPchromosome & base_pair_location >= startPositionGrCh38 & base_pair_location <= endPositionGrCh38)

    #Convert SNP positions to GRCh37
    variantIds <- gwasAASK$variant_id
    #Filter out empty, non-rsid SNPs
    variantIds <- variantIds[variantIds != ""]

    # Query ensembl for GRCh37 positions for the variant IDs
    ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org")
    martOutput <- getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),
      filters = "snp_filter",
      values = variantIds,
      mart = ensembl,
      useCache = FALSE)

    colnames(martOutput) <- c("variant_id", "GRCh37_chr", "GRCh37_position")
    # Merge the GRCh37 positions with your original dataframe
    annotated_gwasAASK <- gwasAASK %>%
      left_join(martOutput, by = "variant_id")

    #Filter out SNPs that did not have a GRCh37 position
    annotated_gwasAASK <- annotated_gwasAASK[!is.na(annotated_gwasAASK$GRCh37_chr),]

    #Further filter out the region to only include the original window now that the liftover is complete
    startPosition <- SNPposition - windowSize
    endPosition <- SNPposition + windowSize
    annotated_gwasAASK <- subset(annotated_gwasAASK, GRCh37_position >= startPosition & GRCh37_position <= endPosition)

    #Remove rare variants as the AASK dataset has their allele frequencies coded as 0 or 1
    annotated_gwasAASK <- annotated_gwasAASK[annotated_gwasAASK$effect_allele_frequency > 0 & annotated_gwasAASK$effect_allele_frequency < 1,]

    #Remove ambiguous duplicate SNP as coloc doesn't support them
    annotated_gwasAASK <- annotated_gwasAASK[!duplicated(annotated_gwasAASK$variant_id),]

    #Export the GRCh37 annotated GWAS
    write.table(annotated_gwasAASK, file = outputFilepath, sep = "\t", quote = FALSE, row.names = FALSE)

    return(annotated_gwasAASK)
  }

  gwasAASK <- convert_AASK_to_GRCh37(somamerID, sentinelSNP)

  # Prepare the data for colocalisation
  gwasVIKINGcoloc <- list(
      beta = gwasVIKING$beta,
      varbeta = (gwasVIKING$standard_error)^2,
      type = "quant",
      snp = gwasVIKING$rs_id,
      position = gwasVIKING$base_pair_location,
      MAF = gwasVIKING$effect_allele_frequency,
      N = gwasVIKING$n
  )

  # Annotate AASK data with static sample size of 466
  gwasAASKColoc <- list(
      beta = gwasAASK$beta,
      varbeta = (gwasAASK$standard_error)^2,
      type = "quant",
      snp = gwasAASK$variant_id,
      position = gwasAASK$GRCh37_position,
      MAF = gwasAASK$effect_allele_frequency,
      N = 466
  )

  # Perform colocalisation
  colocResults <- coloc.abf(
      gwasVIKINGcoloc,
      gwasAASKColoc
  )

  #Export the results
  outputFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", somamerID, ".tsv"))
  output <- capture.output(print(colocResults))
  writeLines(output, con = outputFilepath)
}