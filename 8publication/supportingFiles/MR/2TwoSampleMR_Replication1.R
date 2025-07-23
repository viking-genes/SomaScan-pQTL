#!/usr/bin/env R

# --- Setup ---
openGWASApiKey <- "yourAPIKeyHere"

library(data.table)
library(dplyr)
library(TwoSampleMR)

# --- Directories ---
projectRoot <- "./prj_190_viking_somalogic/"

analysisDir <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/1somamerSNPlists/multiSNP")
outputDir <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/2resultsMR")
harmonisedDataOuptutDir <- file.path(outputDir, "harmonisedLocalOutcomes")

# --- Setup OpenGWAS API ---
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
ao <- available_outcomes(opengwas_jwt = openGWASApiKey)

# --- Manually Annotated List of Replication Datasets ---
replicationDatasetsInternal <- list(
  "GALNT4" = c(file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz")),
  "COMMD10" = c("Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/st005_03_cleaned_moth_alldr.tsv"), # External project
  "AAMDC" = c(
    file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90244093_buildGRCh37.tsv.gz"),
    file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310295.tsv.gz")
  ),
  "LRFN4" = c(
    file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/EA4_additive_p1e-5_clumped.txt"),
    file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")
  ),
  "BCL7A" = c(file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310294.tsv.gz")),
  "PROCR" = c(file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/SNP_gwas_mc_merge_nogc.tbl.uniq.gz")),
  "LTK" = c(file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip"))
)

# --- Study Name Mapping ---
studyNames <- list(
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz") = "GIANT_HEIGHT",
  "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/st005_03_cleaned_moth_alldr.tsv" = "eLife_Maternal_longevity",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90244093_buildGRCh37.tsv.gz") = "Forced_Vital_Capacity",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/EA4_additive_p1e-5_clumped.txt") = "Educational_Attainment",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310294.tsv.gz") = "Systolic_Blood_Pressure",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310295.tsv.gz") = "Diastolic_Blood_Pressure",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/SNP_gwas_mc_merge_nogc.tbl.uniq.gz") = "GIANT_BMI",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip") = "type_2_diabetes",
  file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz") = "smoking_cessation"
)


replicationDatasetsExternal <- list(
  "AAMDC" = c("ieu-a-1006"), #Mean platelet volume 22139419
  "LRFN4" = c("ebi-a-GCST90093322"), #Physical activity 34753499
  "B3GAT1" = c("ieu-b-4809"), #Prostate cancer https://data.bris.ac.uk/datasets/aed0u12w0ede20olb0m77p4b9/Genome-wide%20Association%20Study%20of%20Cancer%20Risk%20in%20UK%20Biobank.pdf
  "GIMAP4" = c("ieu-a-302") #Triglyceride levels 24097068
)

#Clumping reference populations
clumpingPopulations <- list(
  "Somalogic" = "AFR",
  "Olink" = "EUR",
  "eQTLGen" = "EUR"
)

# Extracts the outcome data from a locally downloaded GWAS summary statistics file
extract_local_outcome_data <- function(studyPath, exposureDf) {

    #Read in the outcome data
    result <- fread(studyPath)
    
    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$RSID %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "RSID" = "SNP", 
        "CHR" = "chr", 
        "POS" = "pos", 
        "BETA" = "beta.outcome", 
        "SE" = "se.outcome", 
        "N" = "samplesize.outcome", 
        "P" = "pval.outcome", 
        "EFFECT_ALLELE" = "effect_allele.outcome", 
        "OTHER_ALLELE" = "other_allele.outcome", 
        "EFFECT_ALLELE_FREQ" = "eaf.outcome"
        )

        #Add outcome and id.outcome columns
        result$outcome <- "GIANT_HEIGHT"
        result$id.outcome <- "GIANT_HEIGHT"
    }
    
    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/st005_03_cleaned_moth_alldr.tsv")) {
        #Extract SNP present in the exposure data
        result <- result[result$rsid %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "rsid" = "SNP", 
        "chr" = "chr", 
        "pos" = "pos", 
        "a1" = "effect_allele.outcome", 
        "a0" = "other_allele.outcome", 
        "n" = "samplesize.outcome", 
        "beta1" = "beta.outcome", 
        "se" = "se.outcome", 
        "p" = "pval.outcome", 
        "freq1" = "eaf.outcome"
        )

        #Add outcome and id.outcome columns
        result$outcome <- "eLife_Maternal_longevity"
        result$id.outcome <- "eLife_Maternal_longevity"
    }

    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90244093_buildGRCh37.tsv.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$variant_id %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "variant_id" = "SNP", 
        "chromosome" = "chr", 
        "base_pair_location" = "pos", 
        "effect_allele" = "effect_allele.outcome", 
        "other_allele" = "other_allele.outcome", 
        "N" = "samplesize.outcome", 
        "p_value" = "pval.outcome", 
        "effect_allele_frequency" = "eaf.outcome"
        )

        result$outcome <- "Forced_Vital_Capacity"
        result$id.outcome <- "Forced_Vital_Capacity"

        #Annotate with beta and standard error according to https://www.nature.com/articles/ng.3538, supplementary note 2
        #beta = z-score/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
        #se = 1/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
        result$beta.outcome <- result$Zscore/sqrt(2*result$effect_allele_frequency*(1-result$effect_allele_frequency)*(result$N+result$Zscore**2))
        result$se.outcome <- 1/sqrt(2*result$effect_allele_frequency*(1-result$effect_allele_frequency)*(result$N+result$Zscore**2))
    }

    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/EA4_additive_p1e-5_clumped.txt")) {
        #Extract SNP present in the exposure data
        result <- result[result$rsID %in% exposureDf$SNP,]

        #EAF is not provided for this study
        nameMapping <- c(
        "rsID" = "SNP", 
        "Chr" = "chr", 
        "BP" = "pos", 
        "Effect_allele" = "effect_allele.outcome", 
        "Other_allele" = "other_allele.outcome", 
        "Beta" = "beta.outcome", 
        "SE" = "se.outcome", 
        "P" = "pval.outcome"
        )

        result$outcome <- "Educational_Attainment"
        result$id.outcome <- "Educational_Attainment"
    }

    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310294.tsv.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$rs_id %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "chromosome" = "chr", 
        "base_pair_location" = "pos", 
        "effect_allele" = "effect_allele.outcome", 
        "other_allele" = "other_allele.outcome", 
        "beta" = "beta.outcome", 
        "standard_error" = "se.outcome", 
        "effect_allele_frequency" = "eaf.outcome", 
        "p_value" = "pval.outcome", 
        "rs_id" = "SNP", 
        "n" = "samplesize.outcome"
        )

        result$outcome <- "Systolic_Blood_Pressure"
        result$id.outcome <- "Systolic_Blood_Pressure"
    }
    
    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/downloadedOutcomes/GCST90310295.tsv.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$rs_id %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "chromosome" = "chr", 
        "base_pair_location" = "pos", 
        "effect_allele" = "effect_allele.outcome", 
        "other_allele" = "other_allele.outcome", 
        "beta" = "beta.outcome", 
        "standard_error" = "se.outcome", 
        "effect_allele_frequency" = "eaf.outcome", 
        "p_value" = "pval.outcome", 
        "rs_id" = "SNP", 
        "n" = "samplesize.outcome"
        )

        result$outcome <- "Diastolic_Blood_Pressure"
        result$id.outcome <- "Diastolic_Blood_Pressure"
    }
    
    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/SNP_gwas_mc_merge_nogc.tbl.uniq.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$SNP %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        # "chromosome" = "chr", #Missing column
        # "base_pair_location" = "pos", #Missing column
        "A1" = "effect_allele.outcome", 
        "A2" = "other_allele.outcome", 
        "b" = "beta.outcome", 
        "se" = "se.outcome", 
        "Freq1.Hapmap" = "eaf.outcome", 
        "p" = "pval.outcome", 
        "SNP" = "SNP", 
        "N" = "samplesize.outcome"
        )

        result$outcome <- "GIANT_BMI"
        result$id.outcome <- "GIANT_BMI"
    }

    if (studyPath == file.path(projectRoot, "/Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz")) {
        #Extract SNP present in the exposure data
        result <- result[result$RSID %in% exposureDf$SNP,]

        #Define the name mapping
        nameMapping <- c(
        "CHR" = "chr",
        "POS" = "pos",
        "EFFECT_ALLELE" = "effect_allele.outcome", 
        "OTHER_ALLELE" = "other_allele.outcome", 
        "BETA" = "beta.outcome", 
        "SE" = "se.outcome", 
        "AF_1000G" = "eaf.outcome", 
        "P" = "pval.outcome", 
        "RSID" = "SNP", 
        "N" = "samplesize.outcome"
        )

        result$outcome <- "smoking_cessation"
        result$id.outcome <- "smoking_cessation"

        #Remove the chr prefix from the chromosome column
        result$CHR <- gsub("chr", "", result$CHR)
    }

    #Rename the columns
    result <- result %>%
    rename_at(vars(names(nameMapping)), ~ nameMapping[.])

    #Export the harmonised data
    fwrite(result, paste0(harmonisedDataOuptutDir, "/harmonised_", result$id.outcome[[1]], ".tsv"), sep = "\t", quote = F, na = "NA")

    return (result)
}

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

for (filename in list.files(analysisDir)) {

    #Remove the .csv from the filename
    exposureName = sub("\\.tsv$", "", filename)

    #Get the gene name
    geneName = sub("_.*$", "", exposureName)

    #Get the assay name
    assayName = sub(".*_", "", exposureName)

    print(paste0("Processing ", exposureName))

    #Skip the protein if it had not been selected for replication in the external or internal datasets
    if (!geneName %in% c(names(replicationDatasetsInternal), names(replicationDatasetsExternal))) { 
        next
    }

    #Read in the exposure data
    exposure <- fread(paste0(analysisDir, "/", filename))
    
    #Add the sample size column with 466 as the sample size if the samplesize column does not exist
    if (!"samplesize" %in% colnames(exposure)) {
        exposure$samplesize <- 466
    }

    #Add the phenotype and id column with the filename
    exposure$phenotype <- geneName
    exposure$id <- exposureName

    exposure <- data.frame(exposure)
    
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

    #Get the studies that are to be used for the MR analysis
    replicationOutcomes <- c(replicationDatasetsExternal[[geneName]], replicationDatasetsInternal[[geneName]])
    
    for (studyPath in replicationOutcomes) {
          
        #Check if the file has already been processed
        if (studyPath %in% names(studyNames)) {
            study <- studyNames[studyPath]
        } else {
            study <- studyPath
        }
         
        #Skip the study if the results have already been generated
        if (file.exists(paste0(outputDir, "/MRresult_", exposureName, '_', study, ".tsv"))) {
            next
        }
        
        #Check if the study is stored in a file or is an OpenGWAS study
        if (grepl("/", studyPath)) {
            #Read in and harmonise locally stored outcome data
            outcome <- extract_local_outcome_data(studyPath, exposure)
            
        } else {

            #Extract the outcome data for the selected studies if it is in OpenGWAS
            outcome <- extract_outcome_data(snps = exposure$SNP, outcomes = study, opengwas_jwt = openGWASApiKey)
        }

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
