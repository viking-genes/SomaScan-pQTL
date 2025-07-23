# --- Load required libraries ---
library(biomaRt)
library(data.table)

# --- Define project root and paths ---
projectRoot <- "./prj_190_viking_somalogic/"

gwasUKBBdir <- file.path(
  projectRoot,
  "Scripts/Publication/supportingFiles/replication/Sun2023SummaryStats"
)

supplementary1Df <- file.path(
  projectRoot,
  "Scripts/Publication/supplementary1.xlsx"
)

outputDir <- file.path(
  projectRoot,
  "Scripts/Publication/supportingFiles/replication/MR/0regionData"
)

UKBBMetadataDir <- file.path(
  projectRoot,
  "Scripts/Publication/supportingFiles/replication/Sun2023SummaryStats/metadata"
)


windowSize <- 1000000


# Construct UKBB table around the queried sentinel SNP
convert_UKBB_to_GRCh37 <- function(uniprotID, VIKINGSentinelRow) {

    #Extract the TSS and chromosome from the VIKING data
    TSSchromosome <- VIKINGSentinelRow$Chromosome
    TSSpositionGrCh37 <- as.numeric(VIKINGSentinelRow$'Somamer targeted protein TSS')

    #Read the .gz file of the sentinel SNP chromosome
    folderName <- sub("^[^_]*_", "", filename)
    folderName <- paste0(gwasUKBBdir, "unzipped/")
    
    #Unzip the file
    untar(tarfile = paste0(gwasUKBBdir, filename), exdir = folderName)

    #Update folderName to the unzipped folder
    # Remove the first segment before the underscore to align with the unzipped folder name
    modifiedFilename <- sub("^[^_]+_", "", filename)
    #Remove the .tar extension
    modifiedFilename <- sub("\\.tar$", "", modifiedFilename)
    folderName <- paste0(folderName, modifiedFilename)

    #Retrieve the filename matching the chromosome of the sentinel SNP
    unzippedFiles <- list.files(folderName)
    unzippedFiles <- unzippedFiles[grep(paste0("chr", TSSchromosome, "_"), unzippedFiles)]
    unzippedFiles <- unzippedFiles[1]

    gwasUKBBFull <- fread(cmd = paste("gunzip -c", paste0(folderName, "/", unzippedFiles)), sep = " ", header = TRUE)

    #Create a column of GrCh37 positions taking the 2nd element of the ID column
    gwasUKBBFull$GrCh37Position <- as.numeric(sapply(strsplit(gwasUKBBFull$ID, ":"), "[[", 2))

    #Filter the UKBB data to the expanded region around the TSS of the measured protein
    gwasUKBB <- subset(gwasUKBBFull, CHROM == TSSchromosome & GrCh37Position >= TSSpositionGrCh37 - windowSize & GrCh37Position <= TSSpositionGrCh37 + windowSize)

    #Create a p-value column, converting from -logp
    gwasUKBB$p_value <- 10^(-gwasUKBB$LOG10P)
    
    #Annotate the UKBB data with rsid using UKBB metadata file
    #Read the metadata file with SNP mapping
    metadataFile <- paste0(UKBBMetadataDir, "/", "olink_rsid_map_mac5_info03_b0_7_chr", TSSchromosome, "_patched_v2.tsv.gz")
    metadataFile <- fread(cmd = paste("gunzip -c", metadataFile), sep = "\t", header = TRUE)
    metadataFile <- metadataFile[, c("ID", "rsid")]

    #Map the UKBB data by the ID column in both dataframes, extracting only the rsid column from the metadata
    gwasUKBB <- merge(gwasUKBB, metadataFile, by.x = "ID", by.y = "ID", all.x = TRUE)

    #Export the annotated GWAS data
    outputFilepath <- paste0(outputDir, "Olink/", uniprotID, ".tsv")
    #Create output directory if it does not exist
    dir.create(dirname(outputFilepath), showWarnings = FALSE, recursive = TRUE)
    write.table(gwasUKBB, file = outputFilepath, sep = "\t", quote = FALSE, row.names = FALSE)

    return(gwasUKBB)
}


#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwasUKBBdir
for (filename in list.files(gwasUKBBdir)) {

  #Skip metadata and unzipped folders
  if (filename == "metadata" | filename == "unzipped") {
    next
  }

  #Extract the uniprot ID from the file name
  uniprotID <- strsplit(filename, "_")[[1]][1]

  #Check if the file has already been processed
  resultFilepath <- paste0(outputDir, 'Olink/', uniprotID, ".tsv")
  if (file.exists(resultFilepath)) {
    next
  }

  #Print progress
  print(paste("Processing", filename))
  
  #Extract the VIKING sentinel SNP for the protein
  somamerRow <- supplementary1[supplementary1$Uniprot == uniprotID,]

  #Find the UKBB sentinel SNP, convert the data to GRCh37 and export the window size region
  gwasUKBB <- convert_UKBB_to_GRCh37(uniprotID, somamerRow)
}