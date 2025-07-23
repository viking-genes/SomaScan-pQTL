#VEP command:
#./vep --af --buffer_size 500 --check_existing --distance 5000 --most_severe --plugin MPC,[path_to]/fordist_constraint_official_mpc_values.txt.gz --pubmed --regulatory --species homo_sapiens --transcript_version --cache --input_file [input_data] --output_file [output_file] --port 3337

import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File Paths ---
df = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/retrieveWithinLDForPublishedpQTL/VEPoutput.txt')
outputDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/retrieveWithinLDForPublishedpQTL')

# --- Load Data ---
df = pd.read_csv(df, sep='\t')

# --- Extract Rows Matching Specific Consequences ---
def extractConsequence(df, consequenceList, columnsToExtract):
    result = []
    for i in range(len(df)):
        if df.loc[i, 'Consequence'] in consequenceList:
            result.append(df.loc[i, :])

    result = pd.DataFrame(result)
    result.reset_index(drop=True, inplace=True)
    result = result[columnsToExtract]
    
    return result

# --- Keep Only Unique rsIDs ---
def uniqueRsid(df, columnName):
    uniqueValues = set()
    extractedRows = []

    for _, row in df.iterrows():
        value = row[columnName]
        if value not in uniqueValues:
            uniqueValues.add(value)
            extractedRows.append(row)

    return pd.DataFrame(extractedRows)

# --- Filter for Missense Variants ---
df = extractConsequence(df, consequenceList=['missense_variant'], columnsToExtract=['#Uploaded_variation'])

# --- Deduplicate by rsID ---
df = uniqueRsid(df, '#Uploaded_variation')

# --- Export Result ---
outputFile = os.path.join(outputDir, '2extractMissense.txt')
df.to_csv(outputFile, header=False, index=False)
