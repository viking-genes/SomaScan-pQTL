import os
import pandas as pd

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

# Filepaths
measuredProteinPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/measuredProteins/measuredProteinsIntermediateFile.tsv'
)
uniprotMappingPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/measuredProteins/idmapping_2024_06_27.tsv'
)
exportPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/measuredProteins/previouslyMeasuredProteins.tsv'
)

# --- Load Data ---
measuredProteinDf = pd.read_csv(measuredProteinPath, sep='\t')
uniprotMappingDf = pd.read_csv(uniprotMappingPath, sep='\t')

# --- Setup Output Columns ---
measuredProteinDf['Updated UniprotID'] = ''
measuredProteinDf['Gene Name'] = ''
measuredProteinDf['Alternate Gene Names'] = ''
measuredProteinDf['Full Protein Name'] = ''

# --- Mapping Corrections for Merged/Deprecated Uniprot IDs ---
uniprotCorrectionDict = {
    'Q8IXS6': 'Q9Y2D5',
    'Q8NF90': 'P12034',
    'Q8WWJ7': 'P30203',
    'P04745': 'P0DUB6',
    'P30042': 'P0DPI2',
    'P01233': 'P0DN86',
    'Q8NFS9': 'Q8N0V5',
    'P01916': 'P04440',
    'P04232': 'P04440',
    'P13763': 'P04440',
    'P08107': 'P0DMV8',
    'J3QR46': 'Q15053',
    'P50224': 'P0DMM9',
    'Q9Y4X1': 'P0DTE4'
}

# --- Annotate Each Row ---
for index in measuredProteinDf.index:
    uniprotList = measuredProteinDf.at[index, 'UniprotID'].split('~')

    updatedUniprotIds = []
    proteinNames = []
    geneNames = []
    geneSynonyms = []

    for uniprotId in uniprotList:
        # Apply correction if necessary
        correctedId = uniprotCorrectionDict.get(uniprotId, uniprotId)

        # Match against Uniprot mapping table
        matchIndexList = uniprotMappingDf[uniprotMappingDf['Entry'].str.contains(correctedId)].index.tolist()

        if len(matchIndexList) == 0:
            print('Error: UniprotID', uniprotId, 'not found in Uniprot mapping file.')
            exit()

        matchIndex = matchIndexList[0]
        proteinName = uniprotMappingDf.loc[matchIndex, 'Protein names']
        geneName = uniprotMappingDf.loc[matchIndex, 'Gene Names (primary)']
        geneSynonym = uniprotMappingDf.loc[matchIndex, 'Gene Names (synonym)']

        # Fill NAs
        proteinNames.append(proteinName if pd.notna(proteinName) else 'NA')
        geneNames.append(geneName if pd.notna(geneName) else 'NA')
        geneSynonyms.append(geneSynonym if pd.notna(geneSynonym) else 'NA')
        updatedUniprotIds.append(correctedId)

    # Store combined annotation
    measuredProteinDf.at[index, 'Updated UniprotID'] = '~'.join(updatedUniprotIds)
    measuredProteinDf.at[index, 'Full Protein Name'] = '~'.join(proteinNames)
    measuredProteinDf.at[index, 'Gene Name'] = '~'.join(geneNames)
    measuredProteinDf.at[index, 'Alternate Gene Names'] = '~'.join(geneSynonyms)

# --- Identify Non-Unique IDs (for manual inspection) ---
nonUniqueRows = measuredProteinDf[measuredProteinDf.duplicated(subset='Updated UniprotID')]
print(nonUniqueRows)

# --- Export Annotated Table ---
measuredProteinDf.to_csv(exportPath, sep='\t', index=False)
