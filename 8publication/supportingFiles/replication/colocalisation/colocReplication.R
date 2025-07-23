# --- Load Required Libraries ---
library(data.table)
library(coloc)
library(ieugwasr)
library(tidyr)
library(dplyr)
library(biomaRt)

# --- Define Project Root ---
projectRoot <- "./prj_190_viking_somalogic"

# --- Define Directories and Filepaths ---
aaskDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/colocalisation/AARC_GRCh37")
ukbbDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/colocalisation/UKBB_GRCh37")
eqtlGenDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract")
gcstDatasetDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/GCSTdatasets")
replicationDatasetDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication")
colocOutputDirectory <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/colocalisation/colocResultsReplication")
supplementaryTable1Path <- file.path(projectRoot, "Scripts/Publication/supplementary1.xlsx")
mrReplicationTablePath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/5replicationDf.xlsx")

#Read in the df files
supplementary1Df <- readxl::read_excel(supplementary1Df)
MRReplicationDf <- readxl::read_excel(MRReplicationDf)

# Extracts the outcome data from a locally downloaded GWAS summary statistics file
extract_local_outcome_data <- function(outcomeID, geneName) {

    #Match the outcomeID to the correct file
    if (outcomeID == "GIANT_HEIGHT") {
        studyPath <- paste(replicationDir, "/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz", sep = "")
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
        
    } else if (outcomeID == "eLife_Maternal_longevity") {
        studyPath <- paste(replicationDir, "/st005_03_cleaned_moth_alldr.tsv", sep = "")

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

    } else if (outcomeID == "Forced_Vital_Capacity") {
        studyPath <- paste(replicationDir, "/GCST90244093_buildGRCh37.tsv.gz", sep = "")
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

    } else if (outcomeID == "Educational_Attainment") {
        studyPath <- paste(replicationDir, "/EA4_additive_p1e-5_clumped.txt", sep = "")
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
    } else if (outcomeID == "Systolic_Blood_Pressure") {
        studyPath <- paste(replicationDir, "/GCST90310294.tsv.gz", sep = "")
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
    } else if (outcomeID == "Diastolic_Blood_Pressure") {
        studyPath <- paste(replicationDir, "/GCST90310295.tsv.gz", sep = "")
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
    } else if (outcomeID == "GIANT_BMI") {
        studyPath <- paste(replicationDir, "/SNP_gwas_mc_merge_nogc.tbl.uniq.gz", sep = "")
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
    } else if (outcomeID == "smoking_cessation") {
        studyPath <- paste(replicationDir, "/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz", sep = "")
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

    } else {
        
        #Take only the last item split by -
        studyID <- strsplit(outcomeID, "-")[[1]]

        # Take the last item
        studyID <- tail(studyID, n = 1)

        #Check if there is a file in the replication directory with outcomeID somewhere in the name
        studyPath <- paste(GCSTdatasetDir, "/", grep(studyID, list.files(GCSTdatasetDir), value = TRUE), sep = "")
        
        #The study may be in the discovery outcome directory. If it is not, query Open GWAS
        if (file.exists(studyPath[1]) & startsWith(studyID, "GCST")) {

            #Sample size is not present in the dataframe so adding it manually
            sampleSizeList <- list(
                "ebi-a-GCST004608" = 169545,
                "ebi-a-GCST004599" = 164454,
                "ebi-a-GCST90000618" = 496946,
                "ebi-a-GCST90002386" = 408112,
                "ebi-a-GCST90012009" = 21758,
                "ebi-a-GCST90014290" = 34461,
                "ebi-a-GCST004609" = 170494,
                "ebi-a-GCST004611" = 170761,
                "ebi-a-GCST004627" = 171643,
                "ebi-a-GCST004632" = 171748,
                "ebi-a-GCST004633" = 171542,
                "ebi-a-GCST006097" = 377234,
                "ebi-a-GCST006696" = 412937,
                "ebi-a-GCST006867" = 655666,
                "ebi-a-GCST010723" = 105248,
                "ebi-a-GCST90093322" = 89683
            )
            
            #Read in the outcome data
            result <- fread(studyPath)

            #Remove rows where beta.outcome is not a number
            result <- result[!is.na(result$beta),]

            #Get the chr, pos of the SNP from the supplementary1
            snp_data <- supplementary1Df[supplementary1Df$HUGO == geneName, c("Chromosome", "Position")]
            chromosome <- snp_data$Chromosome[[1]]
            position <- snp_data$Position[[1]]

            #Filter the data to only include SNPs within 1Mb of the pQTL
            result <- result[chromosome == result$chromosome & abs(position - result$base_pair_location) < 1000000,]

            #Remove MAF = 0 and MAF = 1 rows
            result <- result[result$effect_allele_frequency > 0 & result$effect_allele_frequency < 1,]

            #Remove duplicate SNP, leaving only the first instance
            result <- result[!duplicated(result$variant_id),]

            #Remove missing SNP, MAF values
            result <- result[!is.na(result$variant_id),]
            result <- result[!is.na(result$effect_allele_frequency),]
            
            #Set the type of the outcome data
            result <- list(
                beta = result$beta,
                varbeta = (result$standard_error)^2,
                type = "quant",
                snp = result$variant_id,
                MAF = result$effect_allele_frequency,
                pvalues = result$p_value,
                N = sampleSizeList[outcomeID][[1]]
            )


            return (result)

        } else {

            #Get the chr, pos of the SNP from the supplementary1
            snp_data <- supplementary1Df[supplementary1Df$HUGO == geneName, c("Chromosome", "Position")]
            chromosome <- snp_data$Chromosome[[1]]
            position <- snp_data$Position[[1]]
            
            #Query Open GWAS for the outcomeID
            snp_data <- as.data.frame(ieugwasr::associations(
                id = outcomeID,
                variants = paste0(chromosome, ":", position - 1000000, "-", position + 1000000),
                proxies = 0)
            )

            #Convert the MAF to be between 0 and 0.5
            snp_data$eaf <- ifelse(snp_data$eaf > 0.5, 1 - snp_data$eaf, snp_data$eaf)
            
            #Add the varbeta column
            snp_data$varbeta <- (snp_data$se)^2

            #Remove duplicate SNP, leaving only the first instance
            snp_data <- snp_data[!duplicated(snp_data$rsid),]
            
            #Remove rows with no MAF data
            snp_data <- snp_data[!is.na(snp_data$eaf),]
            #Remove rows with MAF = 0 or MAF = 1
            snp_data <- snp_data[snp_data$eaf > 0 & snp_data$eaf < 1, ]
            # Check for non-positive values in varbeta and drop them (ieu-b-85)
            snp_data <- snp_data[snp_data$varbeta > 0,]

            #Apply fixes to the data
            if (outcomeID == "ieu-b-105") {
                #N is missing from the data so adding it manually
                snp_data$n <- 353315
            } else if (outcomeID == "ieu-b-111") {
                #N is missing from the data so adding it manually
                snp_data$n <- 441016
            } else if (outcomeID == "ieu-b-4865") {
                #N is missing from the data so adding it manually
                snp_data$n <- 199569
            }
         
            #Rename the columns for coloc
            result <- list(
                beta = snp_data$beta,
                varbeta = snp_data$varbeta,
                type = "quant",
                snp = snp_data$rsid,
                position = snp_data$position,
                MAF = snp_data$eaf,
                N = snp_data$n,
                pvalues = snp_data$p
            )

            return (result)
        }
    }
    
    #Read in the outcome data
    result <- fread(studyPath)

    #Apply fixes to the data
    if (outcomeID == "Forced_Vital_Capacity") {
        
        #Annotate with beta and standard error according to https://www.nature.com/articles/ng.3538, supplementary note 2
        #beta = z-score/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
        #se = 1/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
        result$beta.outcome <- result$Zscore/sqrt(2*result$effect_allele_frequency*(1-result$effect_allele_frequency)*(result$N+result$Zscore**2))
        result$se.outcome <- 1/sqrt(2*result$effect_allele_frequency*(1-result$effect_allele_frequency)*(result$N+result$Zscore**2))
    }

    if (outcomeID == "smoking_cessation") {

        #Remove the chr prefix from the chromosome column
        result$CHR <- gsub("chr", "", result$CHR)
    }

    #Rename the columns with the standardized names
    result <- result %>%
    rename_at(vars(names(nameMapping)), ~ nameMapping[.])

    #Remove rows where beta.outcome is not a number
    result <- result[!is.na(result$beta.outcome),]

    #Remove MAF = 0 and MAF = 1 rows
    result <- result[result$eaf.outcome > 0 & result$eaf.outcome < 1,]

    #Remove duplicate SNP, leaving only the first instance
    result <- result[!duplicated(result$SNP),]

    #Set the type of the outcome data
    result <- list(
        beta = result$beta.outcome,
        varbeta = (result$se.outcome)^2,
        type = "quant",
        snp = result$SNP,
        MAF = result$eaf.outcome,
        N = result$samplesize.outcome,
        pvalues = result$pval.outcome
    )

    #Change type to binary for smoking cessation
    if (outcomeID == "smoking_cessation") {
        result$type <- "cc"
        #Assign proportion of cases in the outcome data for binary outcomes
        result$s <- 0.53 # Calculated from supplementary table 1 PMID 36477530
    }

    return (result)
}

perform_colocalisation <- function(assayName, geneName, outcomeID, uniprotID = NULL, somamerID = NULL) {
    print(paste(geneName, assayName, outcomeID))
  # Load the exposure data depending on the assay
  if (assayName == "eQTLGen") {
    exposureData <- fread(paste0(eQTLGendir, "/", geneName, ".tsv"))
    
    # Annotate eQTLGen with beta and standard error
    exposureData$beta <- exposureData$Zscore / sqrt(2 * exposureData$AssessedAllele_freq * 
                                                     (1 - exposureData$AssessedAllele_freq) * 
                                                     (exposureData$NrSamples + exposureData$Zscore^2))
    exposureData$standard_error <- 1 / sqrt(2 * exposureData$AssessedAllele_freq * 
                                            (1 - exposureData$AssessedAllele_freq) * 
                                            (exposureData$NrSamples + exposureData$Zscore^2))
    
    exposureData <- list(
      beta = exposureData$beta,
      varbeta = (exposureData$standard_error)^2,
      type = "quant",
      snp = exposureData$SNP,
      position = exposureData$SNPPos,
      MAF = exposureData$AssessedAllele_freq,
      N = exposureData$NrSamples
    )
    
  } else if (assayName == "Olink") {
    exposureData <- fread(paste0(UKBBdir, "/", uniprotID, ".tsv"))
    
    exposureData <- list(
      beta = exposureData$BETA,
      varbeta = (exposureData$SE)^2,
      type = "quant",
      snp = exposureData$rs_id,
      position = exposureData$GRCh37_position,
      MAF = exposureData$A1FREQ,
      N = exposureData$N
    )
    
  } else if (assayName == "Somalogic") {
    exposureData <- fread(paste0(AASKdir, "/", somamerID, ".tsv"))
    
    exposureData <- list(
      beta = exposureData$beta,
      varbeta = (exposureData$standard_error)^2,
      type = "quant",
      snp = exposureData$variant_id,
      position = exposureData$GRCh37_position,
      MAF = exposureData$effect_allele_frequency,
      N = rep(466, length(exposureData$variant_id))
    )
  } else {
    stop("Invalid assayName provided.")
  }
  
  # Read the outcome data
  outcomeData <- extract_local_outcome_data(outcomeID, geneName)

  # Perform colocalisation
  colocResults <- coloc.abf(exposureData, outcomeData)
  
  # Save the results as .txt
  output <- capture.output(print(colocResults))
  writeLines(output, con = paste0(outputDir, "/", geneName, "_", assayName, "_", outcomeID, ".txt"))
  
  return(colocResults)
}



#Create the output directory if it doesn't exist
if (!dir.exists(outputDir)) {
    dir.create(outputDir)
}

#Perform colocalisation for each pair in the MRReplicationDf
for (rowIndex in 1:nrow(MRReplicationDf)) {

    print(rowIndex)

    exposureName <- MRReplicationDf[rowIndex, "Replication MR Exposure"][[1]]
    geneName <- strsplit(exposureName, "_")[[1]][1]
    assayName <- strsplit(exposureName, "_")[[1]][2]
    somamerID <- supplementary1Df[supplementary1Df$HUGO == geneName, "somamerID"][[1]]
    uniprotID <- supplementary1Df[supplementary1Df$HUGO == geneName, "Uniprot"][[1]]

    outcomeID <- MRReplicationDf[rowIndex, "Replication MR Stage 1 Study ID"][[1]]
    #If outcomeID is NA, that means Stage 2 was performed with Discovery Outcome
    if (outcomeID == "NA") {
        outcomeID <- MRReplicationDf[rowIndex, "Discovery MR Study ID"][[1]]
    #If outcomeID is not NA, also perform colocalisation with the Discovery Outcomes
    } else {
        outcomeIDDiscovery <- MRReplicationDf[rowIndex, "Discovery MR Study ID"][[1]]
        if (!is.na(outcomeIDDiscovery)) {

            #Split the outcomeID by | in case there are multiple outcomes
            outcomeIDDiscovery <- strsplit(outcomeIDDiscovery, "\\|")[[1]]
            for (outcome in outcomeIDDiscovery) {

                #Check if this has already been processed
                if (file.exists(paste0(outputDir, "/", geneName, "_", assayName, "_", outcome, ".txt"))) {
                    next
                }

                perform_colocalisation(assayName, geneName, outcome, uniprotID, somamerID)
            }
        }
    }

    #If outcomeID is educational attainment, skip it as there is no allele frequency provided for the dataset
    if (outcomeID == "Educational_Attainment") {
        next
    }

    #If outcomeID is ieu-a-1006, skip it as the data has no allele frequency
    if (outcomeID == "ieu-a-1006") {
        next
    }

    #If outcomeID is ebi-a-GCST90012009, skip it as it is galectin-3 proteomics study and is not MR significant
    if (outcomeID == "ebi-a-GCST90012009") {
        next
    }

    #Check if this has already been processed
    if (file.exists(paste0(outputDir, "/", geneName, "_", assayName, "_", outcomeID, ".txt"))) {
        next
    }

    perform_colocalisation(assayName, geneName, outcomeID, uniprotID, somamerID)
}
