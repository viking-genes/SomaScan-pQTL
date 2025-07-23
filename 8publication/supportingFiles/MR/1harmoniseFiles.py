import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

# --- Define input/output directories ---
inputDir = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/MR/0regionData'
)

outputDir = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/MR/1somamerSNPlists/multiSNP'
)

# --- Load data ---
colocResultPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/colocalisation/colocResultDf.tsv'
)
dfColoc = pd.read_csv(colocResultPath, sep='\t')

vikingSummaryPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supplementary1.xlsx'
)
dfVIKING = pd.read_excel(vikingSummaryPath)


for QTLSource in ["eQTLGen", "Somalogic", "Olink"]:

    for filename in os.listdir(os.path.join(inputDir, QTLSource)):
        
        #Standardise the files and their column names
        if QTLSource == "eQTLGen":

            geneName = filename.split('.')[0]
            print('Processing', geneName, "in", QTLSource)

            #Read the summary statistics file
            df = pd.read_csv(os.path.join(inputDir, QTLSource, filename), sep="\t")

        if QTLSource == "Somalogic":

            #Crossreference the filename for the colocalising protein in the colocalisation results table
            somamerID = filename.split('.')[0]
            #Find somamerID index in the colocalisation results table
            index = dfColoc[dfColoc["somamerID"] == somamerID].index[0]
            #Check if the somamerID is a colocalising signal
            geneName = dfColoc.loc[index, "Gene name"]

            print('Processing', geneName, "in", QTLSource)

            #Read the summary statistics file
            df = pd.read_csv(os.path.join(inputDir, QTLSource, filename), sep="\t")

            #Rename columns to standardise
            df.rename(columns = {'variant_id': 'SNP', 'effect_allele': 'AssessedAllele', 'other_allele': 'OtherAllele', 'effect_allele_frequency': 'AssessedAllele_freq', 'p_value': 'Pvalue'}, inplace = True)

            #Add missing columns
            df['NrSamples'] = 466

        if QTLSource == "Olink":
            #Crossreference the filename for the colocalising protein in the colocalisation results table
            uniprotID = filename.split('.')[0]
            #Find somamerID index in the colocalisation results table
            index = dfVIKING[dfVIKING["Uniprot"] == uniprotID].index[0]
            #Check if the somamerID is a colocalising signal
            geneName = dfVIKING.loc[index, "HUGO"]

            print('Processing', geneName, "in", QTLSource)

            #Read the summary statistics file
            df = pd.read_csv(os.path.join(inputDir, QTLSource, filename), sep="\t")

            #Rename columns to standardise
            df.rename(columns = {'rsid': 'SNP', 'BETA': 'beta', 'SE': 'standard_error', 'ALLELE0': 'AssessedAllele', 'ALLELE1': 'OtherAllele', 'A1FREQ': 'AssessedAllele_freq', 'N': 'NrSamples', 'p_value': 'Pvalue'}, inplace = True)

            #Filter out the non-genome-wide significant SNPs
            df = df[df['Pvalue'] < 5e-8]

            
            
        #Check if already processed
        if os.path.exists(os.path.join(outputDir, geneName + '_' + QTLSource + '.tsv')):
            continue

        #Annotate with missing summary stats
        if QTLSource == "eQTLGen":
            df.loc[:, "beta"] = ''
            df.loc[:, "standard_error"] = ''
            
            #Annotate eQTLGen with beta and standard error according to https://www.nature.com/articles/ng.3538, supplementary note 2
            #beta = z-score/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
            #se = 1/SQRT(2*MAF*(1-MAF)*(n+z-score**2))
            for index in df.index:
                Zscore = df.at[index, "Zscore"]
                alleleFreq = df.at[index, "AssessedAllele_freq"]
                sampleSize = df.at[index, "NrSamples"]
                df.at[index, "beta"] = Zscore / (2 * alleleFreq * (1 - alleleFreq) * (sampleSize + Zscore**2))**0.5
                df.at[index, "standard_error"] = 1 / (2 * alleleFreq * (1 - alleleFreq) * (sampleSize + Zscore**2))**0.5

        #Change the column names to match
        df = df.rename(columns = {'SNP': 'variant_id', 'beta': 'beta', 'standard_error': 'standard_error', 'AssessedAllele': 'effect_allele', 'OtherAllele': 'other_allele', 'AssessedAllele_freq': 'effect_allele_frequency', 'Pvalue': 'p_value', 'NrSamples': 'samplesize'})

        #Export the independent SNP list for MR if there is at least one genome-wide significant SNP
        if len(df) > 0:
            df.to_csv(os.path.join(outputDir, geneName + '_' + QTLSource + '.tsv'), sep="\t", index=False)
