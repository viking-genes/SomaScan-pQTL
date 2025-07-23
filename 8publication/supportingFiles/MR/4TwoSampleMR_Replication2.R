#!/usr/bin/env R

# --- API Key ---
openGWASApiKey <- "yourAPIKeyHere"

# --- Load libraries ---
library(data.table)
library(dplyr)
library(TwoSampleMR)

# --- Define paths ---
projectRoot <- './prj_190_viking_somalogic/'

exposureDir <- file.path(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/1somamerSNPlists/multiSNP')
colocFile <- file.path(projectRoot, 'Scripts/Publication/supplementaryColoc.xlsx')
outputDir <- file.path(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/4resultsMR')


#Use the IEU GWAS API to get the available outcomes
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

#Get the available outcomes. Sometimes fails due to the API server being unresponsive
ao <- available_outcomes(opengwas_jwt = openGWASApiKey)

#Clumping reference populations
clumpingPopulations <- list(
  "Somalogic" = "AFR",
  "Olink" = "EUR",
  "eQTLGen" = "EUR"
)

#Read in the coloc file
colocDf <- readxl::read_excel(colocFile)
#Filter out the non-colocalising results
colocDf <- colocDf[colocDf$'PP.H4.abf' >= 0.8,]

clump_and_perform_MR <- function(data, population){
    #Clump the data
    independentSignals <- clump_data(data, clump_kb = 1000000, pop = population, clump_r2 = 0.001)

    #Create a directory to store the harmonised data
    dir.create(paste0(outputDir), showWarnings = FALSE)
    dir.create(paste0(outputDir, "/harmonisedPreMRData"), showWarnings = FALSE)

    write.table(independentSignals, paste0(outputDir, "/harmonisedPreMRData/MRHarmonised_", exposureName, '_', data$id.outcome[[1]], ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t")

    #Run the MR analysis
    result <- mr(independentSignals, method_list = c("mr_ivw", "mr_wald_ratio"))

    #Add the outcome sample size (samplesize.outcome) to the results by matching on the outcome ID
    result <- merge(result, ao, by.x = "id.outcome", by.y = "id", all.x = T)

    return (result)
}

#Loop through the exposures and perform the MR analysis
exposureFiles <- list.files(exposureDir, pattern = "*.tsv")

for (exposureFile in exposureFiles) {

    #Read in the exposure
    exposure <- read.table(paste0(exposureDir, "/", exposureFile), header = T, sep = "\t")

    #Get the assay name
    exposureName <- sub(".tsv", "", exposureFile)
    assayName = sub(".*_", "", exposureName)

    #Add the phenotype and id column with the filename
    exposure$phenotype <- exposureName
    exposure$id <- exposureName
    
    #Format the exposure data
    exposure <- format_data(
        dat=exposure,
        snp_col = "variant_id",
        beta_col = "beta",
        se_col = "standard_error",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        eaf_col = "effect_allele_frequency",
        pval_col = "p_value",
        samplesize_col = "samplesize",
        phenotype_col="phenotype",
        id_col="id",
        type="exposure"
    )

    #Get the studies that are to be used for the MR analysis from the coloc file
    geneName <- strsplit(exposureName, "_")[[1]][1]
    studies <- colocDf[colocDf$'HUGO' == geneName,]$'Open GWAS Dataset ID'

    #Loop through the studies and perform the MR analysis
    for (study in studies){

        #Skip the study if the results have already been generated
        if (file.exists(paste0(outputDir, "/MRresult_", exposureName, '_', study, ".tsv"))) {
            next
        }

        #Get outcome data
        outcome <- extract_outcome_data(snps = exposure$SNP, outcomes = study, opengwas_jwt = openGWASApiKey)

        #Harmonise the data
        exposureAndOutcomeDf <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
        #Filter out mr_keep == FALSE as these SNPs cannot be used in the MR analysis
        exposureAndOutcomeDf <- exposureAndOutcomeDf[exposureAndOutcomeDf$mr_keep == TRUE,]
        
        #Skip the study if there are no SNPs to perform MR on
        if (nrow(exposureAndOutcomeDf) == 0) {
            print(paste0("No SNPs to perform MR on for ", exposureName, " and ", study))
            next
        }

        #Clump the data and perform MR
        result <- clump_and_perform_MR(exposureAndOutcomeDf, clumpingPopulations[[assayName]])

        #Export the results
        fwrite(result, paste0(outputDir, "/MRresult_", exposureName, '_', study, ".tsv"), sep = "\t", quote = F, na = "NA")
    }
}

