import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
inputPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/uniprotMappingTable.txt')
outputDir = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/partialBiomart/')
chunkSize = 499

# --- Create output directory if needed ---
os.makedirs(outputDir, exist_ok=True)

# --- Load UniProt mapping table ---
df = pd.read_csv(inputPath, sep='\t')

# --- Extract column of interest ---
uniprotIds = df['To'].dropna().tolist()

# --- Split and export chunks ---
totalChunks = (len(uniprotIds) + chunkSize - 1) // chunkSize

for i in range(totalChunks):
    chunk = uniprotIds[i * chunkSize : (i + 1) * chunkSize]
    chunkDf = pd.DataFrame(chunk)
    outputPath = os.path.join(outputDir, f'partial_{i + 1:02d}.txt')
    chunkDf.to_csv(outputPath, index=False, header=False)

print(f'✅ Successfully exported {totalChunks} chunks of UniProt IDs to: {outputDir}')
