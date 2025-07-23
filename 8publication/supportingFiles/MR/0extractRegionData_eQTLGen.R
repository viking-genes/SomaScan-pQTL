# --- Load required libraries ---
library(biomaRt)
library(data.table)
library(httr)
library(readr)
library(dplyr)

# --- Define project root and paths ---
projectRoot <- "./prj_190_viking_somalogic/"

gwaseQTLGendir <- file.path(
  projectRoot,
  "Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract"
)

supplementary1Df <- file.path(
  projectRoot,
  "Scripts/Publication/supplementary1.xlsx"
)

outputDir <- file.path(
  projectRoot,
  "Scripts/Publication/supportingFiles/replication/MR/0regionData"
)

# --- Parameters ---
windowSize <- 1000000


#Only process the proteins that colocalised with VIKING
colocalisingSignals <- c("AAMDC", "B3GAT1", "LRFN4", "GALNT4", "NIF3L1", "NTAQ1")

#Add non-colocalising proteins. CYP2C19 seems to be not present in the eQTLGen data
colocalisingSignals <- c(colocalisingSignals, "PROCR", "BCL7A", "COMMD10", "GIMAP4", "ITGA4", "LTK", "SCPEP1", "ITGA4|ITGB1")

#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwaseQTLGendir
for (filename in list.files(gwaseQTLGendir)) {

  #Extract the uniprot ID from the file name
  geneName <- strsplit(filename, "\\.")[[1]][1]

  #Check if the file has already been processed
  resultFilepath <- paste0(outputDir, 'eQTLGen/', geneName, ".tsv")
  if (file.exists(resultFilepath)) {
    next
  }

  #Check if the protein colocalised with VIKING
  if (!(geneName %in% colocalisingSignals)) {
    next
  }
  
  #Print progress
  print(paste("Processing", filename))
  
  #Extract the TSS and Chr for the protein
  somamerRow <- supplementary1[supplementary1$HUGO == geneName,]
  TSSpositionGrCh37 <- as.numeric(somamerRow$'Somamer targeted protein TSS')
  SNPchromosome <- somamerRow$Chromosome

  #Process ITGA4|ITGB1 separately
  if (geneName == "ITGA4|ITGB1") {
    TSSpositionGrCh37 <- 182321929
    SNPchromosome <- 2
  }

  #Get the window size region
  gwaseQTLGen <- fread(paste0(gwaseQTLGendir, "/", filename))
  gwaseQTLGen <- gwaseQTLGen %>% filter(SNPChr == SNPchromosome & SNPPos >= TSSpositionGrCh37 - windowSize & SNPPos <= TSSpositionGrCh37 + windowSize)

  #Extract genome-wide significant SNPs
  gwaseQTLGen <- gwaseQTLGen %>% filter(Pvalue < 5e-8)

  #Create directory if it doesn't exist
  dir.create(file.path(outputDir, 'eQTLGen'), showWarnings = FALSE)

  #Export the data if there are any SNPs
  if (nrow(gwaseQTLGen) > 0) {
    write.table(gwaseQTLGen, resultFilepath, sep = "\t", row.names = FALSE, quote = FALSE)
  }

}