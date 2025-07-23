#module load igmm/apps/R/3.6.3

# --- Load Required Libraries ---
library(data.table)
library(coloc)
library(ieugwasr)
library(tidyr)
library(dplyr)
library(biomaRt)

# --- Define Project Root ---
projectRoot <- "./prj_190_viking_somalogic"

# --- Define Directory Paths ---
gwasVikingDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/GWASCatalogUpload/convertedGWASFiles")
gwasEqtlGenDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract")
supplementaryTable1Path <- file.path(projectRoot, "Scripts/Publication/supplementary1.xlsx")

#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwasAACSdir and perform colocalisation
for (filename in list.files(gwaseQTLGendir)) {

  #Print progress
  print(paste("Processing", filename))

  #Extract the Gene name from the file name
  geneName <- strsplit(filename, "\\.")[[1]][1]

  #Find the row corresponding to the somamer ID and retrieve the sentinel SNP
  somamerRow <- supplementary1[supplementary1$HUGO == geneName,]
  sentinelSNP <- somamerRow$SNP
  somamerID <- somamerRow$somamerID

  #Check if the file has already been processed
  resultFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", somamerID, "_eQTLGen.tsv"))
  if (file.exists(resultFilepath)) {
    next
  }

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

  #Read the eQTLGen summary statistics stored in .tsv files
  eQTLGen <- fread(paste0(gwaseQTLGendir, "/", filename), sep = "\t", header = TRUE)

  #Annotate eQTLGen with beta and standard error according to https://www.nature.com/articles/ng.3538, supplementary note 2
  #beta = z-score/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
  #se = 1/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
  eQTLGen$beta <- eQTLGen$Zscore/sqrt(2*eQTLGen$AssessedAllele_freq*(1-eQTLGen$AssessedAllele_freq)*(eQTLGen$NrSamples+eQTLGen$Zscore^2))
  eQTLGen$standard_error <- 1/sqrt(2*eQTLGen$AssessedAllele_freq*(1-eQTLGen$AssessedAllele_freq)*(eQTLGen$NrSamples+eQTLGen$Zscore^2))

  gwaseQTLGenColoc <- list(
      beta = eQTLGen$beta,
      varbeta = (eQTLGen$standard_error)^2,
      type = "quant",
      snp = eQTLGen$SNP,
      position = eQTLGen$SNPPos,
      MAF = eQTLGen$AssessedAllele_freq,
      N = eQTLGen$NrSamples
  )

  # Perform colocalisation
  colocResults <- coloc.abf(
      gwasVIKINGcoloc,
      gwaseQTLGenColoc
  )

  #Export the results
  outputFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", somamerID, "_eQTLGen.tsv"))
  output <- capture.output(print(colocResults))
  writeLines(output, con = outputFilepath)
}