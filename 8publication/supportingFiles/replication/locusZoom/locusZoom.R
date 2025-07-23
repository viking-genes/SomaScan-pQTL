# Load necessary modules
#python/2.7.10
#bcftools/1.20

# --- Define Project Root ---
projectRoot <- "./prj_190_viking_somalogic"

# --- Define Explicit Paths for Studies ---
vikingGwasPath <- file.path(projectRoot, "GWAS/Full/p05_comb_chr")
eqtlgenPath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract")
somalogicPath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/Surapaneni2022Data/fullStats")
olinkPath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/MR/0regionData/Olink")
discoveryMRPath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/discovery")
replicationMRPath <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/datasets/replication")

outputDir <- file.path(projectRoot, "Scripts/Publication/supportingFiles/replication/locusZoom/output")


# Define a list with protein names and corresponding rsids
#studies = c('VIKING', 'eQTLGen', 'Somalogic', 'Olink')
proteinData <- list(

  AAMDC = list(
    studies = c("str24216__30.tsv.gz", "AAMDC.tsv", "24216-30.gz", "Q9H7C9.tsv"),
    discovery = c('ieu-b-39.vcf.gz', 'ukb-b-7992.vcf.gz', 'ieu-b-105.vcf.gz', 'ukb-b-7953.vcf.gz', 'ebi-a-GCST004599.vcf.gz', '34187551-GCST90014290-GO_0007568.h.tsv.gz'),
    replication = c('ieu-a-1006.vcf.gz', 'GCST90244093_buildGRCh37.tsv.gz', 'GCST90310295.tsv.gz'),
    rsid = 'rs72941336',
    chromosome = 11,
    windowSize = 1000
  ),

  B3GAT1 = list(
    studies = c("str21770__18.tsv.gz", "B3GAT1.tsv", "21770-18.gz", ""),
    discovery = c('ieu-b-85.vcf.gz', 'ieu-b-32.vcf.gz'),
    replication = c('ieu-b-4809.vcf.gz'),
    rsid = 'rs78760579',
    chromosome = 11,
    windowSize = 100
  ),

  BCL7A = list(
    studies = c("str21331__19.tsv.gz", "BCL7A.tsv", "21331-19.gz", "Q4VC05.tsv"),
    discovery = c('ieu-b-38.vcf.gz'),
    replication = c('GCST90310294.tsv.gz'),
    rsid = 'rs1169084',
    chromosome = 12,
    windowSize = 1000
  ),

  COMMD10 = list(
    studies = c("str23257__14.tsv.gz", "COMMD10.tsv", "23257-14.gz", ""),
    discovery = c('ebi-a-GCST006696.vcf.gz'),
    replication = c('st005_03_cleaned_moth_alldr.tsv'),
    rsid = 'rs56953556',
    chromosome = 5,
    windowSize = 500
  ),

  CYP2C19 = list(
    studies = c("str25464__1.tsv.gz", "", "25464-1.tsv", ""),
    discovery = c('32242144-GCST90000618-EFO_0004631.h.tsv.gz'),
    replication = c(''),
    rsid = 'rs80159642',
    chromosome = 10,
    windowSize = 500
  ),

  GALNT4 = list(
    studies = c("str21722__21.tsv.gz", "GALNT4.tsv", "21722-21.tsv", ""),
    discovery = c('ieu-b-31.vcf.gz', 'ieu-b-33.vcf.gz', 'ebi-a-GCST004609.vcf.gz', 'ebi-a-GCST004608.vcf.gz', 'ieu-a-89.vcf.gz', 'ukb-b-6235.vcf.gz'),
    replication = c('GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz'),
    rsid = 'rs7960289',
    chromosome = 12,
    windowSize = 1000
  ),

  GIMAP4 = list(
    studies = c("str24684__7.tsv.gz", "GIMAP4.tsv", "24684-7.tsv", ""),
    discovery = c('ebi-a-GCST004627.vcf.gz', 'ebi-a-GCST004632.vcf.gz', 'ieu-b-111.vcf.gz', 'ieu-b-32.vcf.gz', 'ebi-a-GCST004633.vcf.gz', '33067605-GCST90012009-EFO_0008137.h.tsv.gz'),
    replication = c('ieu-a-302.vcf.gz', 'GCST90239661.h.tsv.gz'),
    rsid = 'rs62491812',
    chromosome = 7,
    windowSize = 200
  ),

  ITGA4ITGB1 = list(
    studies = c("str21688__50.tsv.gz", "ITGA4|ITGB1.tsv", "21688-50.tsv", ""),
    discovery = c('ebi-a-GCST004611.vcf.gz', '32888494-GCST90002386-EFO_0010700.h.tsv.gz'),
    replication = c(''),
    rsid = 'rs767500',
    chromosome = 2,
    windowSize = 200
  ),

  LRFN4 = list(
    studies = c("str21696__80.tsv.gz", "LRFN4.tsv", "21696-80.tsv", ""),
    discovery = c('ukb-b-19842.vcf.gz', 'ieu-a-1239.vcf.gz', 'ebi-a-GCST006097.vcf.gz', 'ukb-b-18275.vcf.gz'),
    replication = c('GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz', 'EA4_additive_p1e-5_clumped.txt', '34753499-GCST90093322-EFO_0008002.h.tsv.gz'),
    rsid = 'rs2077432',
    chromosome = 11,
    windowSize = 1000
  ),

  LTK = list(
    studies = c("str23181__2.tsv.gz", "LTK.tsv", "23181-2.gz", ""),
    discovery = c('ukb-b-14609.vcf.gz', 'ukb-b-10753.vcf.gz', 'ebi-a-GCST006867.vcf.gz'),
    replication = c('Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip'),
    rsid = 'rs1473781',
    chromosome = 15,
    windowSize = 500
  ),
  
  NIF3L1 = list(
    studies = c("str23371__5.tsv.gz", "NIF3L1.tsv", "23371-5.gz", ""),
    discovery = c('ebi-a-GCST010723.vcf.gz'),
    replication = c(''),
    rsid = 'rs10931931',
    chromosome = 2,
    windowSize = 500
  ),
  
  NTAQ1 = list(
    studies = c("str24649__11.tsv.gz", "NTAQ1.tsv", "24649-11.gz", ""),
    discovery = c('ieu-b-4865.vcf.gz'),
    replication = c(''),
    rsid = 'rs13258747',
    chromosome = 8,
    windowSize = 300
  ),

  PROCR = list(
    studies = c("str24023__35.tsv.gz", "PROCR.tsv", "24023-35.gz", "Q9UNN8.tsv"),
    discovery = c('ukb-b-20531.vcf.gz', 'ukb-b-12854.vcf.gz', 'ukb-b-19953.vcf.gz', 'ukb-b-8909.vcf.gz'),
    replication = c('SNP_gwas_mc_merge_nogc.tbl.uniq.gz'),
    rsid = 'rs6060241',
    chromosome = 20,
    windowSize = 1300
  ),

  SCPEP1 = list(
    studies = c("str23203__3.tsv.gz", "SCPEP1.tsv", "", "Q9HB40.tsv"),
    discovery = c('ieu-b-30.vcf.gz', 'ebi-a-GCST004627.vcf.gz', 'ieu-b-32.vcf.gz', 'ieu-b-34.vcf.gz', 'ieu-b-31.vcf.gz'),
    replication = c(''),
    rsid = 'rs59721290',
    chromosome = 17,
    windowSize = 1300
  )
)

nonStandardStudies <- list(
  '34187551-GCST90014290-GO_0007568.h.tsv.gz' = list(
    outputName = "GCST90014290",
    rsidColumnName = "hm_rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "hm_chrom",
    islogTransformed = FALSE
  ),
  '32242144-GCST90000618-EFO_0004631.h.tsv.gz' = list(
    outputName = "GCST90000618",
    rsidColumnName = "hm_rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "hm_chrom",
    islogTransformed = FALSE
  ),
  '32888494-GCST90002386-EFO_0010700.h.tsv.gz' = list(
    outputName = "GCST90002386",
    rsidColumnName = "hm_rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "hm_chrom",
    islogTransformed = FALSE
  ),
  '33067605-GCST90012009-EFO_0008137.h.tsv.gz' = list(
    outputName = "GCST90012009",
    rsidColumnName = "hm_rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "hm_chrom",
    islogTransformed = FALSE
  ),
  'GCST90244093_buildGRCh37.tsv.gz' = list(
    outputName = "GCST90244093",
    rsidColumnName = "variant_id",
    pValueColumnName = "p_value",
    chromosomeColumnName = "chromosome",
    islogTransformed = FALSE
  ),
  'GCST90310295.tsv.gz' = list(
    outputName = "GCST90310295",
    rsidColumnName = "rs_id",
    pValueColumnName = "p_value",
    chromosomeColumnName = "chromosome",
    islogTransformed = FALSE
  ),
  'GCST90310294.tsv.gz' = list(
    outputName = "GCST90310294",
    rsidColumnName = "rs_id",
    pValueColumnName = "p_value",
    chromosomeColumnName = "chromosome",
    islogTransformed = FALSE
  ),
  'st005_03_cleaned_moth_alldr.tsv' = list(
    outputName = "maternal_longevity_eLife_Timmers",
    rsidColumnName = "rsid",
    pValueColumnName = "p",
    chromosomeColumnName = "chr",
    islogTransformed = FALSE
  ),
  'GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz' = list(
    outputName = "height_GIANT_Yengo",
    rsidColumnName = "RSID",
    pValueColumnName = "P",
    chromosomeColumnName = "CHR",
    islogTransformed = FALSE
  ),
  'GCST90239661.h.tsv.gz' = list(
    outputName = "GCST90239661",
    rsidColumnName = "rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "chromosome",
    islogTransformed = FALSE
  ),
  'GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR_without_UKB.txt.gz' = list(
    outputName = "smoking_cessation_36477530",
    rsidColumnName = "RSID",
    pValueColumnName = "P",
    chromosomeColumnName = "CHR",
    islogTransformed = FALSE
  ),
  'EA4_additive_p1e-5_clumped.txt' = list(
    outputName = "educational_attainment_35361970",
    rsidColumnName = "rsID",
    pValueColumnName = "P",
    chromosomeColumnName = "Chr",
    islogTransformed = FALSE
  ),
  '34753499-GCST90093322-EFO_0008002.h.tsv.gz' = list(
    outputName = "GCST90093322",
    rsidColumnName = "hm_rsid",
    pValueColumnName = "p_value",
    chromosomeColumnName = "hm_chrom",
    islogTransformed = FALSE
  ),
  'Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip' = list(
    outputName = "type2_diabetes",
    rsidColumnName = "None",
    pValueColumnName = "Pval",
    chromosomeColumnName = "Chromsome",
    islogTransformed = FALSE,
    positionColumnName = "Position",
    positionArchitecture = "GRCh37",
    allele1ColumnName = "EffectAllele",
    allele2ColumnName = "NonEffectAllele",
    targetRsidPosition = 41818917
  ),
  'SNP_gwas_mc_merge_nogc.tbl.uniq.gz' = list(
    outputName = "BMI_GIANT",
    rsidColumnName = "SNP",
    pValueColumnName = "p",
    chromosomeColumnName = "None",
    islogTransformed = FALSE
  )
)

extract_data_from_downloaded_file <- function(studyName, tsvFilepath, protein, studyFilename, pathToFiles) {
  # Change the working directory to the path of the files
  setwd(pathToFiles)

  # Create .tsv file from .vcf file or extract it from .gz for column renaming before locuszoom
  if (grepl(".vcf", studyFilename)) {
    
    #Convert the .vcf file to .tsv, keeping only the necessary data
    cmdConvertToTsv <- sprintf("bcftools query -f '%%CHROM\t%%POS\t%%ID\t%%REF\t%%ALT\t[%%LP]\n' %s/%s > %s.tsv", pathToFiles, studyFilename, studyName)
    system(cmdConvertToTsv)

    #Read in the created .tsv file
    df <- read.table(tsvFilepath, header = FALSE, sep = "\t")
    #Rename the columns
    colnames(df) <- c("CHR", "BP", "RSID", "REF", "ALT", "LOGP")

    #Extract the chromosome that is cis to the protein
    df <- df[df$CHR == proteinData[[protein]]$chromosome, ]

  } else {

    #Extract the file to studyName.tsv if it is zipped
    if (grepl(".gz", studyFilename)) {
      cmdExtractTsv <- sprintf("gunzip -c %s/%s > %s.tsv", pathToFiles, studyFilename, studyName)
      system(cmdExtractTsv)
    } else if (grepl(".zip", studyFilename)) {
      cmdExtractZip <- sprintf("unzip -p %s/%s > %s.tsv", pathToFiles, studyFilename, studyName)
      system(cmdExtractZip)

      # Otherwise copy the .tsv file to the output directory, even if it is a symbolic link for processing
    } else {
      cmdCopyTsv <- sprintf("cp %s/%s %s.tsv", pathToFiles, studyFilename, studyName)
      system(cmdCopyTsv)
    }

    #Read in the .tsv file
    df <- read.table(tsvFilepath, header = TRUE, sep = "\t")
    #Rename the columns according to the nonStandardStudies list
    if (studyFilename %in% names(nonStandardStudies)) {
      rsidColumnName <- nonStandardStudies[[studyFilename]]$rsidColumnName
      pValueColumnName <- nonStandardStudies[[studyFilename]]$pValueColumnName
      colnames(df)[colnames(df) == rsidColumnName] <- "RSID"
      colnames(df)[colnames(df) == pValueColumnName] <- "LOGP"
      colnames(df)[colnames(df) == nonStandardStudies[[studyFilename]]$chromosomeColumnName] <- "CHR"

      #Extract the chromosome that is cis to the protein if annotation is present
      #Check if CHR column is present in the dataframe
      if ("CHR" %in% colnames(df)) {
        #Remove rows with NA values in the chromosome column as they interfere with filtering
        df <- df[!is.na(df$CHR), ]

        #Remove the chr prefix if present (smoking_cessation)
        df$CHR <- gsub("chr", "", df$CHR)

        #Match the data type
        df$CHR <- as.character(df$CHR)
        chrToExtract <- as.character(proteinData[[protein]]$chromosome)
        df <- df[df$CHR == chrToExtract, ]
      }

      #Transform p-values to -log10 if necessary
      if (!nonStandardStudies[[studyFilename]]$islogTransformed) {
        df$LOGP <- -log10(df$LOGP)
      }

      #Annotate with rsid if it is not present in the dataframe
      if (nonStandardStudies[[studyFilename]]$rsidColumnName == "None") {
        
        # Main function to annotate rsIDs by position and allele information using Ensembl
        annotate_rsid_by_position <- function(df, chromosomeColumn, positionColumn, alleleColumn1, alleleColumn2, assembly = "GRCh37") {

            library(biomaRt)
            library(dplyr)

            # Specify the host based on the assembly
            if (assembly == "GRCh37") {
                hostName = "https://grch37.ensembl.org"
            } else if (assembly == "GRCh38") {
                hostName = "https://www.ensembl.org"
            } else {
                stop("Invalid assembly specified. Please use 'GRCh37' or 'GRCh38'")
            }

            # Connect to the Ensembl database
            ensemblConnect <- useMart(host = hostName, 
                                    biomart = "ENSEMBL_MART_SNP",
                                    dataset = "hsapiens_snp")

            # Create chr:pos:pos format for the query
            df <- within(df, ensemblQuery <- paste(df[[chromosomeColumn]], ":", df[[positionColumn]], ":", df[[positionColumn]], sep = ""))

            # Initialize the RSID column
            df$RSID <- NA

            #Set batch size to 50 for querying Ensembl to avoid timeouts
            batchSize <- 50
            numberOfBatches <- ceiling(nrow(df) / batchSize)

            for (batch in 1:numberOfBatches) {
              
              print(paste("Querying Ensembl API, ", batch, "/", numberOfBatches))

              startIndex <- (batch - 1) * batchSize + 1
              endIndex <- min(nrow(df), batch * batchSize)

              currentBatch <- df[startIndex:endIndex, ]

              # Query Ensembl for the rsIDs
              retrieveRsid <- getBM(attributes = c('refsnp_id', 'allele', 'chrom_start'), 
                                  filters = 'chromosomal_region',
                                  values = currentBatch$ensemblQuery, 
                                  mart = ensemblConnect)

              #Merge the results with the original dataframe
              for (i in 1:nrow(currentBatch)) {
                  allele1 <- currentBatch[[alleleColumn1]][i]
                  allele2 <- currentBatch[[alleleColumn2]][i]

                  #Subset matching position rows from the retrieved query
                  matchingRows <- retrieveRsid[retrieveRsid$chrom_start == currentBatch[[positionColumn]][i], ]

                  #Subset matching allele rows by matching allele info
                  matchingRowsAllele <- matchingRows[
                    sapply(strsplit(matchingRows$allele, "/"), function(alleles) {
                      all(c(allele1, allele2) %in% alleles) & length(intersect(alleles, c(allele1, allele2))) == 2
                    }),
                  ]
                  
                  # If there are multiple rsIDs for the same position and allele, print warning and select the first one
                  if (nrow(matchingRowsAllele) > 1) {
                      warning(paste("Multiple rsIDs found for position", currentBatch[[positionColumn]][i], "and alleles", allele1, ",", allele2))

                      print(matchingRows)

                      print("Selecting the first matched row")

                      # Select the first row
                      matchingRowsAllele <- matchingRowsAllele[1, ]

                  } else if (nrow(matchingRowsAllele) == 0) {
                      warning(paste("No RSID found for position", currentBatch[[positionColumn]][i], "and alleles", allele1, ",", allele2))
                      
                      print(matchingRows)
                      # Skip to the next row
                      next
                  }

                  # Annotate the RSID with the original dataframe
                  df$RSID[startIndex + i - 1] <- matchingRowsAllele$refsnp_id
              }
            }

            # Remove the ensemblQuery column
            df <- subset(df, select = -ensemblQuery)

            return (df)
        }

        #Reduce the dataframe to the window size
        positionColumnName <- nonStandardStudies[[studyFilename]]$positionColumnName
        targetPosition <- nonStandardStudies[[studyFilename]]$targetRsidPosition

        df[[positionColumnName]] <- as.numeric(df[[positionColumnName]])
        df <- df[abs(df[[positionColumnName]] - targetPosition) < proteinData[[protein]]$windowSize * 1000, ]

        #Filter out non-significant p-values
        df <- df[df$LOGP > 1, ]

        #Annotate the rsid
        df <- annotate_rsid_by_position(df, "CHR", positionColumnName, nonStandardStudies[[studyFilename]]$allele1ColumnName, nonStandardStudies[[studyFilename]]$allele2ColumnName, nonStandardStudies[[studyFilename]]$positionArchitecture)

        #Drop rows without rsid as they cannot be plotted
        df <- df[!is.na(df$RSID), ]
        }
      }
    }

    #Remove all rows that have rsid start with TMP (GCST90239661)
    df <- df[!grepl("^TMP", df$RSID), ]

    #Leave only necessary columns
    df <- df[, c("RSID", "LOGP")]

    #Save the dataframe for locuszoom
    write.table(df, tsvFilepath, sep = "\t", quote = FALSE, row.names = FALSE)

}

# Function to run locuszoom on discovery and replication MR studies
locus_zoom_on_study_list <- function(studyList, protein, rsid, pathToFiles) {
  for (studyFilename in studyList) {
    if (studyFilename != "") {

      # Define output file name
      studyName <- gsub(".vcf.gz", "", studyFilename)
      if (studyName %in% names(nonStandardStudies)) {
        studyName <- nonStandardStudies[[studyName]]$outputName
      }
      outputFilepath <- file.path(outputDir, sprintf("%s_%s_%s.pdf", protein, studyName, rsid))

      #Skip if the output file already exists
      if (file.exists(outputFilepath)) {
        next
      }

      print(paste0('Processing ', studyName, ' for ', protein))

      tsvFilepath <- paste0(pathToFiles, '/', studyName, ".tsv")
      #Check if the tsv file has already been created
      if (!file.exists(tsvFilepath)) {
        #Process the raw data file into locuszoom usable format
        extract_data_from_downloaded_file(studyName, tsvFilepath, protein, studyFilename, pathToFiles)
      }

    #Get the pQTL-specific manually annotated window size for the locuszoom command
    windowSize <- proteinData[[protein]]$windowSize

    # Run locuszoom command for the discovery MR study
    setwd(outputDir)
    cmdLocuszoom <- sprintf("locuszoom --build hg19 --metal '%s' --pop EUR --source 1000G_Nov2014 --refsnp %s --markercol RSID --pvalcol LOGP --flank %skb --plotonly --no-transform --no-date --prefix %s_%s",
                            tsvFilepath, rsid, windowSize, protein, studyName)
    system(cmdLocuszoom)

    }
  }
}


# Function to generate and run locuszoom commands
run_locus_zoom <- function(protein, studies, rsid) {

  # Define output file names
  vikingOutput <- file.path(outputDir, sprintf("%s_VIKING_%s.pdf", protein, rsid))
  eqtlgenOutput <- file.path(outputDir, sprintf("%s_eQTLGen_%s.pdf", protein, rsid))
  somalogicOutput <- file.path(outputDir, sprintf("%s_Somalogic_%s.pdf", protein, rsid))
  olinkOutput <- file.path(outputDir, sprintf("%s_Olink_%s.pdf", protein, rsid))

  # Join paths with study files
  vikingGwasFile <- ifelse(studies[1] != "", file.path(vikingGwasPath, studies[1]), NA)
  eqtlgenFile <- ifelse(studies[2] != "", file.path(eqtlgenPath, studies[2]), NA)
  somalogicFile <- ifelse(studies[3] != "", file.path(somalogicPath, studies[3]), NA)
  olinkFile <- ifelse(studies[4] != "", file.path(olinkPath, studies[4]), NA)
  
  # Get window size for locuszoom
  windowSize <- proteinData[[protein]]$windowSize

  # VIKING GWAS study command
  if (!file.exists(vikingOutput) && !is.na(vikingGwasFile) && file.exists(vikingGwasFile)) {
    cmdViking <- sprintf("locuszoom --build hg19 --metal '%s' --pop EUR --source 1000G_Nov2014 --refsnp %s --markercol rsid --pvalcol p --flank %skb --plotonly --no-date --prefix %s_VIKING",
                          vikingGwasFile, rsid, windowSize, protein)
    cat("Running:", cmdViking, "\n")
    system(cmdViking)
  }
  
  # eQTLGen study command
  if (!file.exists(eqtlgenOutput) && !is.na(eqtlgenFile) && file.exists(eqtlgenFile)) {
    cmdEqtlgen <- sprintf("locuszoom --build hg19 --metal '%s' --pop EUR --source 1000G_Nov2014 --refsnp %s --markercol SNP --pvalcol Pvalue --flank %skb --plotonly --no-date --prefix %s_eQTLGen",
                           eqtlgenFile, rsid, windowSize, protein)
    cat("Running:", cmdEqtlgen, "\n")
    system(cmdEqtlgen)
  }
  
  # Somalogic study command
  if (!file.exists(somalogicOutput) && !is.na(somalogicFile) && file.exists(somalogicFile)) {
    cmdSomalogic <- sprintf("locuszoom --build hg19 --metal '%s' --pop AFR --source 1000G_Nov2014 --refsnp %s --markercol variant_id --pvalcol p_value --flank %skb --plotonly --no-date --prefix %s_Somalogic",
                             somalogicFile, rsid, windowSize, protein)
    cat("Running:", cmdSomalogic, "\n")
    system(cmdSomalogic)
  }

  # Olink study command
  if (!file.exists(olinkOutput) && !is.na(olinkFile) && file.exists(olinkFile)) {
    cmdOlink <- sprintf("locuszoom --build hg19 --metal '%s' --pop EUR --source 1000G_Nov2014 --refsnp %s --markercol rsid --pvalcol p_value --flank %skb --plotonly --no-date --prefix %s_Olink",
                         olinkFile, rsid, windowSize, protein)
    cat("Running:", cmdOlink, "\n")
    system(cmdOlink)
  }

  # Discovery MR study command
  discoveryStudies <- proteinData[[protein]]$discovery
  locus_zoom_on_study_list(discoveryStudies, protein, rsid, discoveryMRPath)

  # # Replication MR study command
  replicationStudies <- proteinData[[protein]]$replication
  locus_zoom_on_study_list(replicationStudies, protein, rsid, replicationMRPath)
}

# Loop through each protein and run the locuszoom commands
for (protein in names(proteinData)) {
  data <- proteinData[[protein]]
  run_locus_zoom(protein, data$studies, data$rsid)
}

cat("All locuszoom commands have been completed.\n")
