import pandas as pd
import numpy as np

# Define parameters
basePairWindow = 300_000

# Load input data
reverseMrResultsPath = './prj_190_viking_somalogic/Scripts/twoSampleMr/5combineReverseMRResults.xlsx'
supplementaryDataPath = './prj_190_viking_somalogic/Scripts/Publication/supplementary1.xlsx'

# Define output path
outputPath = './prj_190_viking_somalogic/Scripts/colocalization/1mergedTwoSampleMRResults.tsv'

reverseMrResultsDf = pd.read_excel(reverseMrResultsPath)
supplementaryDf = pd.read_excel(supplementaryDataPath)

# Initialize new columns in reverse MR dataframe
for column in ['bpStart', 'bpEnd', 'chr', 'somamerID']:
    reverseMrResultsDf[column] = ''

def getAnnotatedRow(targetRow):
    """
    Extracts chromosome and position window info from a supplementary row.
    Also constructs a string version of somamer ID.
    """
    minPosition = targetRow['Position'].min()
    maxPosition = targetRow['Position'].max()
    chromosome = targetRow.iloc[0]['Chromosome']
    somamerId = 'str' + '__'.join(targetRow.iloc[0]['somamerID'].split('-'))

    return minPosition, maxPosition, chromosome, somamerId

# Iterate through each MR result and annotate with supplementary info
for index, row in reverseMrResultsDf.iterrows():
    exposureGene = row['exposure']
    matchingRows = supplementaryDf[supplementaryDf['HUGO'] == exposureGene]

    if matchingRows.empty:
        continue  # No matching supplementary entry found

    if len(matchingRows) > 1:
        # Check for consistent annotations if there are multiple entries
        for column in ['Chromosome', 'cisTrans', 'somamerID']:
            if matchingRows[column].nunique() != 1:
                print('⚠️ Warning: inconsistent entries found for exposure:', exposureGene)
                print(matchingRows)
                break
        else:
            # All entries are consistent; proceed to annotate
            bpStart, bpEnd, chromosome, somamerId = getAnnotatedRow(matchingRows)
            reverseMrResultsDf.at[index, 'bpStart'] = bpStart - basePairWindow
            reverseMrResultsDf.at[index, 'bpEnd'] = bpEnd + basePairWindow
            reverseMrResultsDf.at[index, 'chr'] = chromosome
            reverseMrResultsDf.at[index, 'somamerID'] = somamerId
    else:
        # Single matching entry
        bpStart, bpEnd, chromosome, somamerId = getAnnotatedRow(matchingRows)
        reverseMrResultsDf.at[index, 'bpStart'] = bpStart - basePairWindow
        reverseMrResultsDf.at[index, 'bpEnd'] = bpEnd + basePairWindow
        reverseMrResultsDf.at[index, 'chr'] = chromosome
        reverseMrResultsDf.at[index, 'somamerID'] = somamerId

# Save annotated output
reverseMrResultsDf.to_csv(outputPath, sep='\t', index=False)
