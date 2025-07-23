#!/usr/bin/env R

openGWASApiKey = "yourAPIKey"

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(openxlsx)

# --- Directories ---
projectRoot <- './prj_190_viking_somalogic/'

analysisDir <- file.path(
  projectRoot,
  'Scripts/Publication/supportingFiles/replication/MR/1somamerSNPlists/multiSNP'
)

outputDir <- file.path(
  projectRoot,
  'Scripts/Publication/supportingFiles/replication/MR/outcomeSelection/outcomeSelection'
)

# --- Files ---
supplementary1 <- file.path(
  projectRoot,
  'Scripts/Publication/supplementary1.xlsx'
)

VIKINGResultsDf <- file.path(
  projectRoot,
  'Scripts/Publication/supplementaryColoc.xlsx'
)

#Read in the dataframes
supplementary1 <- readxl::read_excel(supplementary1)
VIKINGResultsDf <- readxl::read_excel(VIKINGResultsDf)

#Use the IEU GWAS API to get the available outcomes
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

ao <- available_outcomes(opengwas_jwt = openGWASApiKey)

#Remove all proteomics studies from the available outcome ID since we are only interested in the medical outcomes
ao <- ao %>% filter(!grepl("prot", id, ignore.case = TRUE))

#Lists of greedy keywords to use for the MR outcome search
keywordList <- list(
  "AAMDC" = list(
    "ebi-a-GCST90014290" = c("epigenetic age", "DNA methylation", "age acceleration"),
    "ieu-b-39" = c("diastolic", "DBP", "hypertension"),
    "ukb-b-7953" = c("forced vital capacity", "FVC", "lung function"),
    "ebi-a-GCST004599" = c("MPV", "platelet")
  ),
  "B3GAT1" = list(
    "ieu-b-85" = c("prostate"),
    "ieu-b-32" = c("lymphocyte")
  ),
  "BCL7A" = list(
    "ieu-b-38" = c("systolic", "SBP", "hypertension")
  ),
  "COMMD10" = list(
    "ebi-a-GCST006696" = c("longevity", "lifespan", "life expectancy")
  ),
  "GALNT4" = list(
    "ieu-b-31" = c("monocyte", "leukocyte", "WBC", "white blood", "white cell", "granulocyte"),
    "ieu-b-33" = c("eosinophil", "acidophil", "WBC", "white blood", "white cell", "granulocyte"),
    "ieu-a-89" = c("height", "stature"),
    "ukb-b-6235" = c("cholecystectomy", "gallbladder", "gall bladder", "gallstone", "gall stone", "cholelithiasis")
  ),
  "ITGA4|ITGB1" = list(
    "ebi-a-GCST90002386" = c("reticulocyte", "retic", "erythrocyte", "red blood", "red cell")
  ),
  "LRFN4" = list(
    "ukb-b-19842" = c("smoking", "cigarette", "nicotine", "tobacco", "smoke"),
    "ieu-a-1239" = c("schooling", "education", "education", "years of school", "college", "university", "degree", "qualification"),
    "ebi-a-GCST006097" = c("physical activity", "exercise", "sport", "fitness", "workout", "training"),
    "ukb-b-18275" = c("hearing", "deafness")
  ),
  "NIF3L1" = list(
    "ebi-a-GCST010723" = c("macula", "retina", "eye", "vision", "sight", "blind")
  ),
  "NTAQ1" = list(
    "ieu-b-4865" = c("testosterone", "androgen")
  )
)

#Some studies used in VIKING MR were measuring the same trait 
sameStudyList <- list(
  "ieu-b-39" = c("ukb-b-7992"),
  "ukb-b-7953" = c("ieu-b-105"),
  "ieu-b-31" = c("ebi-a-GCST004609", "ebi-a-GCST004608"),
  "ieu-b-33" = c("ebi-a-GCST004608"),
  "ebi-a-GCST90002386" = c("ebi-a-GCST004611")
)

#non-UKBB study list
nonUKBBList <- c("ieu-a-89", "ieu-b-85", "ebi-a-GCST90014290")

#List of UKBB studies with non-ukb identifiers or proteomics studies
ukbbPMIDList <- c("29892013", "33067605", "32888493", "34226706", "32888494", "27863252", "34594039", "33303764", "28957414", "33230300", "32042192", "30804560", "33728380", "32986727", "30643251")

for (filename in list.files(analysisDir)) {

    #Skip the Old directory
    if (filename == "Old") {
        next
    }

    #Read in the exposure data
    exposure <- fread(paste0(analysisDir, "/", filename))
    
    #Add the sample size column with 466 as the sample size if the samplesize column does not exist
    if (!"samplesize" %in% colnames(exposure)) {
        exposure$samplesize <- 466
    }

    #Remove the .csv from the filename
    exposureName = sub("\\.tsv$", "", filename)

    #Also remove everything after _ in the filename
    geneName = sub("_.*$", "", exposureName)

    #Add the phenotype and id column with the filename
    exposure$phenotype <- geneName
    exposure$id <- exposureName

    exposure <- data.frame(exposure)
    
    #Format the exposure data
    exp_df <- format_data(
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

    #Retrieve all colocalising and MR GWAS studies for this gene
    get_studies_by_HUGO <- function(data, geneName) {
      result <- data %>%
        filter(HUGO == geneName) %>%
        #Filter coloc significant results
        filter(`PP.H4.abf` > 0.8) %>%
        select(`Open GWAS Dataset ID`) %>%
      
      return(result)
    }

    #Get the studies that passed the colocalisation test for VIKING
    studies <- get_studies_by_HUGO(VIKINGResultsDf, geneName)$`Open GWAS Dataset ID`

    #Filter studies for specified keywords for the gene and study pair ignoring case
    # Define the function to filter based on keyword list
    filter_by_keywords <- function(data, keywordList, study) {
      # Initialize an empty data frame to store the filtered results
      result <- data.frame()

      # Loop through each keyword
      for (keyword in keywordList) {
        # Filter the data where the 'trait' column matches the keyword, but not as part of a larger word
        matchingRows <- data %>% filter(grepl(paste0("\\b", keyword, "\\b"), trait, ignore.case = TRUE))
        # Combine the matching rows with the filtered_data
        result <- rbind(result, matchingRows)
      }

      # Remove duplicate rows from the filtered data
      result <- result %>% distinct()

      #Remove the studies measured in VIKING from the filtered data
      result <- result %>% filter(!id %in% study)

      # Return the filtered data
      return(result)
    }

    #Create dir for Jim to select the most appropriate outcome
    dir.create(paste0(outputDir), showWarnings = FALSE)

    # Loop through each key for this protein in the keyword list and filter the available outcomes
    for (study in names(keywordList[[geneName]])) {
      keywords <- keywordList[[geneName]][[study]]
      
      # Filter the available outcomes for the specified keywords leaving out the originally linked MR studies
      VIKINGstudyList <- c(study, sameStudyList[[study]])
      filteredOutcomes <- filter_by_keywords(ao, keywords, VIKINGstudyList)

      #Filter out studies from the same cohort
      #Remove UKBB studies from the filtered outcomes if the original MR study was UKBB
      if (!(study %in% nonUKBBList)) {
        #Remove UKBB studies from the filtered outcomes by removing the ukb- prefix studies
        filteredOutcomes <- filteredOutcomes %>% filter(!grepl("ukb-", id, ignore.case = TRUE))

        #Search the entire dataframe for keywords related to UK Biobank and remove those studies
        UKBiobankKeywordList <- c("ukb", "UK Biobank")
        filteredOutcomes <- filteredOutcomes %>% filter(!apply(filteredOutcomes, 1, function(row) any(grepl(paste(UKBiobankKeywordList, collapse = "|"), row, ignore.case = TRUE))))

        #Remove the studies that were checked to include UKBB via PMID reference 
        filteredOutcomes <- filteredOutcomes %>% filter(!pmid %in% ukbbPMIDList)
      }

      #Sort by sample size then by number of cases so that finngen studies are at the top
      filteredOutcomes <- filteredOutcomes[order(filteredOutcomes$sample_size, decreasing = T),]
      filteredOutcomes <- filteredOutcomes[order(filteredOutcomes$ncase, decreasing = T),]

      #Add the row with the VIKING MR studies as the first rows and create an empty row below it
      VIKINGstudies <- ao %>% filter(id %in% VIKINGstudyList)

      #Create empty row
      emptyRow <- as.data.frame(matrix("", nrow = 1, ncol = ncol(VIKINGstudies)))
      colnames(emptyRow) <- colnames(VIKINGstudies)

      #Construct the final filtered outcomes
      filteredOutcomes <- rbind(VIKINGstudies, emptyRow, filteredOutcomes)
      
      # Write the filtered outcomes to an excel file
      write.xlsx(filteredOutcomes, paste0(outputDir, '/', study, ".xlsx"))
    }

}
