#Extracts all genome-wide significant pQTL within a defined window size around the transcription start site of the measured protein

# --- Load required libraries ---
library(biomaRt)
library(data.table)
library(dplyr)
library(httr)
library(readr)

# --- Define project root and directories ---
projectRoot <- './prj_190_viking_somalogic/'

gwasAASKdir <- file.path(
  projectRoot,
  'Scripts/Publication/supportingFiles/replication/Surapaneni2022Data/fullStats'
)

supplementary1Df <- file.path(
  projectRoot,
  'Scripts/Publication/supplementary1.xlsx'
)

outputDir <- file.path(
  projectRoot,
  'Scripts/Publication/supportingFiles/replication/MR/0regionData'
)

windowSize <- 1000000

#Only process the proteins that colocalised with VIKING
colocalisingSignals <- c("21688-50", "21696-80", "21722-21", "23257-14", "21331-19", "24216-30", "21770-18")

#Add non-colocalising proteins too
colocalisingSignals <- c(colocalisingSignals, "24023-35", "25464-1", "24684-7", "23181-2", "23371-5", "24649-11", "23203-3")

# Function to get the TSS of the gene in GRCh38
get_GRCh38_TSS <- function(somamerRowVIKING) {

    geneName <- somamerRowVIKING$HUGO

    #For ITGA4|ITGB1, assign ITGA4 as the gene name as the cis pQTL is in that region
    if (geneName == "ITGA4|ITGB1") {
        geneName <- "ITGA4"
    }

    #Query ensembl for GRCh38 position of the gene's TSS
    ensemblGrCh38 <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    martOutputGrCh38 <- getBM(
        attributes = c("chromosome_name", "start_position"),
        filters = "external_gene_name",
        values = geneName,
        mart = ensemblGrCh38,
        useCache = FALSE)

    TSSpositionGrCh38 <- as.numeric(martOutputGrCh38$start_position)
    return(TSSpositionGrCh38)
}

# Checks if the liftover region extends to the limits of the window size
check_liftover_region <- function(df, startPositionGrCh37, windowSize) {

    #Check if the liftovered region extends past the GRCh37 TSS region in both directions
    if (min(df$GRCh37_position) < startPositionGrCh37 - windowSize | max(df$GRCh37_position) > startPositionGrCh37 + windowSize) {
        return(TRUE)
    }

    print(paste("Liftover region does not extend to the limits of the window size for", somamerID))
    print(paste("GRCh38 TSS position:", TSSpositionGrCh38))
    print(paste("GRCh37 TSS position:", startPositionGrCh37))
    print(paste("GRCh37 min position:", min(martOutput$GRCh37_position)))
    print(paste("GRCh37 max position:", max(martOutput$GRCh37_position)))
    return(FALSE)
}


# Construct AASK table around the queried sentinel SNP
convert_AASK_to_GRCh37 <- function(somamerID, somamerRowVIKING) {

    gwasAASKFull <- fread(cmd = paste("gunzip -c", paste0(gwasAASKdir, "/", somamerID, ".gz")), sep = "\t", header = TRUE)
    
    #Extract only the chromosome of interest
    SNPchromosome <- somamerRowVIKING$Chromosome
    gwasAASKFull <- gwasAASKFull[gwasAASKFull$chromosome == SNPchromosome,]

    #Query ensembl for GRCh38 position of the gene's TSS
    TSSpositionGrCh38 <- get_GRCh38_TSS(somamerRowVIKING)

    #Get the TSS in GRCh37
    TSSpositionGrCh37 <- as.numeric(somamerRowVIKING$'Somamer targeted protein TSS')

    #Replace ITGA4/ITGB1 with ITGA4
    if (somamerRowVIKING$HUGO == "ITGA4|ITGB1") {
        TSSpositionGrCh37 <- 182321929
    }

    # Extract larger region for liftover to be reduced later, just in case the region changed between GRCh37 and GRCh38
    startPositionGrCh38 <- TSSpositionGrCh38 - 1.5*windowSize
    endPositionGrCh38 <- TSSpositionGrCh38 + 1.5*windowSize

    gwasAASKQuery <- subset(gwasAASKFull, base_pair_location >= startPositionGrCh38 & base_pair_location <= endPositionGrCh38)

    #Filter out empty, non-rsid SNPs
    variantIds <- gwasAASKQuery$variant_id
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
    
    #Filter out rows with patch chromosome annotations as they are just duplicates
    martOutput <- martOutput[!grepl("PATCH", martOutput$GRCh37_chr),]
    
    #Check if the liftovered region extends to the GRCh37 TSS region
    if (!check_liftover_region(martOutput, TSSpositionGrCh37, windowSize)) {
        return(NULL)
    }

    # Merge the GRCh37 positions with the original dataframe
    annotated_gwasAASK <- gwasAASKQuery %>%
        left_join(martOutput, by = "variant_id")

    #Filter out SNPs that did not have a GRCh37 position
    annotated_gwasAASK <- annotated_gwasAASK[!is.na(annotated_gwasAASK$GRCh37_chr),]

    #Remove rare variants as the AASK dataset has their allele frequencies coded as 0 or 1
    annotated_gwasAASK <- annotated_gwasAASK[annotated_gwasAASK$effect_allele_frequency > 0 & annotated_gwasAASK$effect_allele_frequency < 1,]

    #Further filter out the region to only include the original window now that the liftover is complete
    startPosition <- TSSpositionGrCh37 - windowSize
    endPosition <- TSSpositionGrCh37 + windowSize
    annotated_gwasAASK <- subset(annotated_gwasAASK, GRCh37_position >= startPosition & GRCh37_position <= endPosition)

    #Find the genome-wide significant SNP in the AASK data
    annotated_gwasAASK <- annotated_gwasAASK[annotated_gwasAASK$p_value < 5e-8,]

    #Export the GRCh37 annotated GWAS
    outputFilepath <- paste0(outputDir, "Somalogic/", somamerID, ".tsv")
    #Create output directory if it does not exist
    dir.create(dirname(outputFilepath), showWarnings = FALSE, recursive = TRUE)
    write.table(annotated_gwasAASK, file = outputFilepath, sep = "\t", quote = FALSE, row.names = FALSE)

    return(annotated_gwasAASK)
}


#Read in the supplementary file
supplementary1 <- readxl::read_excel(supplementary1Df)

#Iterate over all the files in gwasAASKdir
for (filename in list.files(gwasAASKdir)) {

  #Extract the somamer ID from the file name
  somamerID <- gsub(".gz", "", filename)

  #Check if the protein colocalised with VIKING
  if (!(somamerID %in% colocalisingSignals)) {
    next
  }

  #Print progress
  print(paste("Processing", filename))

  #Check if the file has already been processed
  resultFilepath <- paste0(outputDir, 'Somalogic/', somamerID, ".tsv")
  if (file.exists(resultFilepath)) {
    next
  }

  #Find the row corresponding to the somamer ID and retrieve the sentinel SNP
  somamerRowVIKING <- supplementary1[supplementary1$somamerID == somamerID,]

  #Find the AASK sentinel SNP, convert the AASK data to GRCh37 and export the window size region
  gwasAASK <- convert_AASK_to_GRCh37(somamerID, somamerRowVIKING)

}