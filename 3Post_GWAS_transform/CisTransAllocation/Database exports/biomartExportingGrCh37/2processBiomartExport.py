import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

mainExportPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/mart_export.txt')
uniprotMappingPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/uniprotMappingTable.txt')
hugoMappingPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/uniprotToHugo.xlsx')
supplementPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/mart_export_supplement.txt')

outputPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/TSSmartExport.csv')

# --- Define Patch/Scaffold Fixes ---
patchToChrMap = {
    'HG999_1_PATCH': 1, 'HG544_PATCH': 10, 'HG183_PATCH': 17, 'HG305_PATCH': 11, 'HG256_PATCH': 11,
    'HG1257_PATCH': 7, 'HG1211_PATCH': 10, 'HG1293_PATCH': 1, 'HG1287_PATCH': 1, 'HG311_PATCH': 10,
    'HG1436_HG1432_PATCH': 'X', 'HG1463_PATCH': 'X', 'HG1497_PATCH': 'X', 'HG1443_HG1444_PATCH': 'X',
    'HG1459_PATCH': 'X', 'HG385_PATCH': 17, 'HG957_PATCH': 3, 'HG75_PATCH': 17, 'HG1079_PATCH': 19
}

# --- Step 1: Clean chromosome names and deduplicate ---
def cleanChromosomes(df):
    def normalizeChromosome(entry):
        entry = str(entry)
        if entry.startswith('HSCHR'):
            chrPart = entry.split('_')[0][5:]
            return chrPart[:2] if chrPart[:2].isdigit() else chrPart
        elif entry.startswith('CHR_HSCHR19'):
            return '19'
        elif entry in patchToChrMap:
            return patchToChrMap[entry]
        print(f'⚠️ Excluded patch-only entry: {entry}')
        return 'PATCH'

    df['Chromosome/scaffold name'] = df['Chromosome/scaffold name'].apply(normalizeChromosome)
    df.drop_duplicates(subset='Gene stable ID', inplace=True, ignore_index=True)
    df.rename(columns={'Gene start (bp)': 'Transcription start site (TSS)'}, inplace=True)
    return df

# --- Step 2: Add UniProt from Ensembl-to-UniProt mapping ---
def mapUniprot(df, mappingDf):
    df['Uniprot'] = ''
    mappingDf.dropna(subset=['From', 'To'], inplace=True)

    for idx, geneId in enumerate(df['Gene stable ID']):
        matched = mappingDf[mappingDf['To'].astype(str).str.contains(geneId)]
        if not matched.empty:
            uniprotIds = matched['From'].dropna().unique()
            df.at[idx, 'Uniprot'] = '|'.join(uniprotIds)
        else:
            print(f'⚠️ No UniProt ID found for gene ID: {geneId}')
    return df

# --- Step 3: Add supplemental Hugo gene name mapping (Ensembl fallback) ---
def addSupplemental(df, supplementDf, hugoDf):
    supplementDf.rename(columns={'Gene start (bp)': 'Transcription start site (TSS)'}, inplace=True)
    supplementDf['Uniprot'] = ''

    hugoDf.dropna(subset=['From', 'To'], inplace=True)

    for idx, geneName in enumerate(supplementDf['Gene name']):
        matches = hugoDf[hugoDf['To'].astype(str).str.contains(geneName)]
        if not matches.empty:
            uniprots = matches['From'].dropna().unique()
            supplementDf.at[idx, 'Uniprot'] = '|'.join(uniprots)
        else:
            print(f'⚠️ No UniProt found for gene symbol: {geneName}')

    supplementDf.drop(columns='Gene name', inplace=True)
    return pd.concat([df, supplementDf], ignore_index=True)


# --- Run pipeline ---

# Load inputs
biomartDf = pd.read_csv(mainExportPath, sep='\t')
supplementDf = pd.read_csv(supplementPath, sep='\t')
uniprotMapDf = pd.read_csv(uniprotMappingPath, sep='\t')
hugoMapDf = pd.read_excel(hugoMappingPath)

# Clean and annotate
biomartDf = cleanChromosomes(biomartDf)
supplementDf = cleanChromosomes(supplementDf)
biomartDf = mapUniprot(biomartDf, uniprotMapDf)
mergedDf = addSupplemental(biomartDf, supplementDf, hugoMapDf)

# Export
mergedDf.to_csv(outputPath, index=False)
print(f'✅ Exported annotated TSS mapping to: {outputPath}')
