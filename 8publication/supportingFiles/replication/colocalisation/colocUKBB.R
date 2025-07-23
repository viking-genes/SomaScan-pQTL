# Load the required libraries
library(data.table)
library(coloc)z
library(ieugwasr)
library(tidyr)
library(dplyr)

# --- Define Project Root ---
projectRoot <- "./prj_190_viking_somalogic"

#Define directories
gwasVIKINGdir <- file.path(projectRoot, "/Scripts/Publication/supportingFiles/GWASCatalogUpload/convertedGWASFiles/")
gwasUKBBdir <- file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/Sun2023SummaryStats/")
supplementary1Df <- file.path(projectRoot, "/Scripts/Publication/supplementary1.xlsx")

#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwasAACSdir and perform colocalisation
for (filename in list.files(gwasUKBBdir)) {

  #Only process .tar files
  if (!grepl(".tar", filename)) {
    next
  }

  #Print progress
  print(paste("Processing", filename))

  #Extract the Uniprot ID from the file name
  uniprotID <- strsplit(filename, "_")[[1]][1]

  #Check if the file has already been processed
  resultFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", uniprotID, ".tsv"))
  if (file.exists(resultFilepath)) {
    next
  }

  #Find the row corresponding to the somamer ID and retrieve the sentinel SNP
  somamerRow <- supplementary1[supplementary1$Uniprot == uniprotID,]
  sentinelSNP <- somamerRow$SNP
  somamerID <- somamerRow$somamerID

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

  # Construct UKBB table around the queried sentinel SNP
  convert_UKBB_to_GRCh37 <- function(uniprotID, sentinelSNP) {
    
    #Create output directory if it doesn't exist already
    dir.create(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/UKBB_GRCh37", showWarnings = FALSE))

    outputFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/UKBB_GRCh37/", uniprotID, ".tsv"))

    #Check if the file already exists
    if (file.exists(outputFilepath)) {
      return(fread(outputFilepath, sep = "\t", header = TRUE))
    }

    #Unpack the UKBB .tar file into ./unzipped
    dir.create(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/Sun2023SummaryStats/unzipped", showWarnings = FALSE))
    untar(paste0(gwasUKBBdir, filename), exdir = paste0(gwasUKBBdir, "unzipped/"))

    #Read the .gz file of the sentinel SNP chromosome
    folderName <- sub("^[^_]*_", "", filename)
    folderName <- paste0(gwasUKBBdir, "unzipped/", gsub(".tar", "", folderName))
    
    #Retrieve the filename matching the chromosome of the sentinel SNP
    unzippedFiles <- list.files(folderName)
    unzippedFiles <- unzippedFiles[grep(paste0("chr", SNPchromosome, "_"), unzippedFiles)]
    unzippedFiles <- unzippedFiles[1]

    gwasUKBBFull <- fread(cmd = paste("gunzip -c", paste0(folderName, "/", unzippedFiles)), sep = " ", header = TRUE)

    #Create a column for the GrCh37 positions by splitting the ID column
    gwasUKBBFull$GRCh37_position <- as.numeric(sapply(strsplit(gwasUKBBFull$ID, ":"), "[[", 2))

    #Filter out the region to only include the coloc window
    startPosition <- SNPposition - windowSize
    endPosition <- SNPposition + windowSize
    annotated_gwasUKBB <- subset(gwasUKBBFull, GRCh37_position >= startPosition & GRCh37_position <= endPosition)

    #Annotate the UKBB data with rsid by mapping overlapping SNP with the VIKING dataset
    annotated_gwasUKBB$rs_id <- as.character(NA)

    for (i in 1:nrow(annotated_gwasUKBB)) {
      SNP <- annotated_gwasUKBB[i,]
      overlappingSNP <- gwasVIKINGFiltered[gwasVIKINGFiltered$base_pair_location == SNP$GRCh37_position,]
      if (nrow(overlappingSNP) > 1) {
        print(paste("Multiple overlapping SNPs found for", SNP$GRCh37_position, "in", uniprotID))
      }
      if (nrow(overlappingSNP) > 0) {
        annotated_gwasUKBB$rs_id[i] <- as.character(overlappingSNP$rs_id)
      }
    }

    #Drop rows with missing rsid
    annotated_gwasUKBB <- annotated_gwasUKBB[!is.na(annotated_gwasUKBB$rs_id),]

    #Drop rows with duplicate rsid as some have been separated by effect allele in the UKBB data
    annotated_gwasUKBB <- annotated_gwasUKBB[!duplicated(annotated_gwasUKBB$rs_id),]

    #Export the GRCh37 annotated GWAS
    write.table(annotated_gwasUKBB, file = outputFilepath, sep = "\t", quote = FALSE, row.names = FALSE)

    return(annotated_gwasUKBB)
  }

  gwasUKBB <- convert_UKBB_to_GRCh37(uniprotID, sentinelSNP)

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

  gwasUKBBColoc <- list(
      beta = gwasUKBB$BETA,
      varbeta = (gwasUKBB$SE)^2,
      type = "quant",
      snp = gwasUKBB$rs_id,
      position = gwasUKBB$GRCh37_position,
      MAF = gwasUKBB$A1FREQ,
      N = gwasUKBB$N
  )

  # Perform colocalisation
  colocResults <- coloc.abf(
      gwasVIKINGcoloc,
      gwasUKBBColoc
  )

  #Export the results
  outputFilepath <- paste0(file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/colocalisation/colocResults/", somamerID, "_UKBB.tsv"))
  output <- capture.output(print(colocResults))
  writeLines(output, con = outputFilepath)
}