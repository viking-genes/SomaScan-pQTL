import pandas as pd
import numpy as np
import os

# --- Define project-relative paths ---
projectRoot = './prj_190_viking_somalogic/'
locationPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportUniprot.csv')
gwasPath = os.path.join(projectRoot, 'GWAS/GWAShitsv2.csv')
rawProteinPath = os.path.join(projectRoot, 'GWAS/viking1_proteomics_somalogic7k_2021_flagged.csv')
uniprotMappingPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExporting/uniprotMappingTable.txt')

# --- Load data ---
locationDf = pd.read_csv(locationPath)
gwasDf = pd.read_csv(gwasPath)
rawDf = pd.read_csv(rawProteinPath, low_memory=False)
uniprotToGeneMapping = pd.read_csv(uniprotMappingPath, sep='\t')


# --- Helper functions ---
def dropUnnamedAndFixProteinCol(df):
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df.rename(columns={'protein': 'SomalogicSeqID'}, inplace=True)
    return df

def extractUniprotDf(rawDf):
    # Extract row 6 (index 5) and row 1 (index 0) for SeqID and UniProt mapping
    tempDf = pd.DataFrame([rawDf.iloc[5], rawDf.iloc[0]])
    tempDf = tempDf.transpose().dropna()
    tempDf.columns = tempDf.iloc[0]
    tempDf = tempDf.iloc[1:].reset_index(drop=True)
    return tempDf

def mapUniprot(df, rawDf):
    uniprotDf = extractUniprotDf(rawDf)
    df['Uniprot'] = df['SomalogicSeqID'].map(
        lambda x: uniprotDf.loc[uniprotDf['SeqId'] == x, 'UniProt'].values[0] 
        if x in uniprotDf['SeqId'].values else np.nan
    )
    return df

def addGeneStableID(df, mappingDf):
    df['GeneStableID'] = ''

    noMatch = []

    for idx, uniprot in df['Uniprot'].items():
        if pd.isna(uniprot):
            noMatch.append(df.loc[idx])
            continue

        targets = uniprot.split('|')
        geneIds = []

        for t in targets:
            matches = mappingDf[mappingDf['From'] == t]['To'].tolist()
            if matches:
                geneIds.extend(matches)
        
        if geneIds:
            df.at[idx, 'GeneStableID'] = '|'.join(sorted(set(geneIds)))
        else:
            noMatch.append(df.loc[idx])

    with open(os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/addGeneStableIDNoMatch.txt'), 'w') as f:
        for row in noMatch:
            f.write(str(row.to_dict()) + '\n')

    return df

def addTSSInfo(df, locationDf):
    df['Transcription start site (TSS)'] = ''
    df['TSS Chromosome'] = ''

    noMatch = []

    for idx, geneId in df['GeneStableID'].items():
        if pd.isna(geneId):
            noMatch.append(df.loc[idx])
            continue

        geneIds = geneId.split('|')
        tssList = []
        chrList = []

        for gid in geneIds:
            subset = locationDf[locationDf['Gene stable ID'] == gid]
            if not subset.empty:
                tss = subset['Transcription start site (TSS)'].values[0]
                chrom = subset['Chromosome/scaffold name'].values[0]
                tssList.append(str(tss))
                chrList.append(str(chrom))

        # Remove 'PATCH' entries if other chromosomes exist
        if any(chrVal != 'PATCH' for chrVal in chrList):
            tssList = [tss for tss, chrVal in zip(tssList, chrList) if chrVal != 'PATCH']
            chrList = [chrVal for chrVal in chrList if chrVal != 'PATCH']

        if tssList and chrList:
            df.at[idx, 'Transcription start site (TSS)'] = '|'.join(tssList)
            df.at[idx, 'TSS Chromosome'] = '|'.join(chrList)
        else:
            noMatch.append(df.loc[idx])

    with open(os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/noMatch.txt'), 'w') as f:
        for row in noMatch:
            f.write(str(row.to_dict()) + '\n')

    return df

# --- Apply pipeline ---
gwasDf = dropUnnamedAndFixProteinCol(gwasDf)
gwasDf = mapUniprot(gwasDf, rawDf)
gwasDf = addGeneStableID(gwasDf, uniprotToGeneMapping)
gwasDf = addTSSInfo(gwasDf, locationDf)

# --- Save ---
outputPath = os.path.join(projectRoot, 'GWAS/GWAShitsv2_annotated.csv')
gwasDf.to_csv(outputPath, index=False)
print(f"✅ Annotated GWAS exported to: {outputPath}")
