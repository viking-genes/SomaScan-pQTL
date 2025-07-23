import os
import pandas as pd

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

# Directories and paths
fullGwasDir = os.path.join(projectRoot, 'GWAS/Full/p05_comb_chr/')
gwasHitDir = os.path.join(projectRoot, 'GWAS/GWAShits/')
rawCsvDirs = [
    os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/'),
    os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/processed/')
]
significantFullExportDir = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/significantHitFullExport/')
finalSubsetDir = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/onlySignficantHits/')
chromosomeTargetPath = os.path.join(projectRoot, 'Post_GWAS_transform/chromosomeTargets.txt')

# Parameters
pvalueThresholdForHitExport = 1e-5
pvalueThresholdGenomeWide = 5e-8

# --- Step 1: Extract hits below p-value threshold and annotate with protein name ---
def extractSuggestiveHits():
    os.makedirs(gwasHitDir, exist_ok=True)
    fileCount = 0
    for fileName in os.listdir(fullGwasDir):
        if not fileName.endswith('.tsv.gz'):
            continue

        filePath = os.path.join(fullGwasDir, fileName)
        df = pd.read_csv(filePath, sep='\t', compression='gzip', low_memory=False)
        fileCount += 1
        print(f'Processing {fileName} ({fileCount}/{len(os.listdir(fullGwasDir))})')

        significantHits = df[df['p'] < pvalueThresholdForHitExport]
        if not significantHits.empty:
            # Annotate with protein name from file name
            proteinName = fileName[3:-7].replace('__', '-')
            significantHits['protein'] = proteinName

            outputPath = os.path.join(gwasHitDir, proteinName + '8e-8.csv')
            significantHits.to_csv(outputPath, index=False)

# --- Step 2: Identify genome-wide significant hits and extract full summaries ---
def extractFullSummariesFromHits():
    os.makedirs(significantFullExportDir, exist_ok=True)
    filesToExtract = []

    # Identify genome-wide significant hits
    for directory in rawCsvDirs:
        for fileName in os.listdir(directory):
            if not fileName.endswith('.csv'):
                continue
            df = pd.read_csv(os.path.join(directory, fileName))
            if (df['p'] < pvalueThresholdGenomeWide).any():
                filesToExtract.append(fileName)

    # Export full .tsv.gz summaries for those hits
    for fileName in filesToExtract:
        if fileName in os.listdir(significantFullExportDir):
            continue
        print(f'Exporting full summary for: {fileName}')
        placeholderPath = os.path.join(significantFullExportDir, fileName)
        with open(placeholderPath, 'w') as tempFile:
            tempFile.write('temp')

        formattedName = 'str' + fileName.replace('-', '__').replace('.csv', '')
        fullPath = os.path.join(fullGwasDir, formattedName + '.tsv.gz')
        df = pd.read_csv(fullPath, sep='\t', compression='gzip', low_memory=False)
        df.to_csv(placeholderPath, sep=' ', index=False)

# --- Step 3: Identify chromosomes with genome-wide significant hits ---
def generateChromosomeTargetList():
    chromosomeTargets = []
    for fileName in os.listdir(significantFullExportDir):
        fullPath = os.path.join(significantFullExportDir, fileName)
        df = pd.read_csv(fullPath, sep=' ', low_memory=False)
        targetChromosomes = df.loc[df['p'] < pvalueThresholdGenomeWide, 'chr'].dropna().unique()
        if len(targetChromosomes) > 0:
            chromosomeTargets.append([fileName, ','.join(map(str, targetChromosomes))])

    os.makedirs(os.path.dirname(chromosomeTargetPath), exist_ok=True)
    with open(chromosomeTargetPath, 'w') as f:
        for entry in chromosomeTargets:
            f.write(' '.join(entry) + '\n')

# --- Step 4: Export only significant hits from specified chromosomes ---
def exportSignificantChromosomes():
    os.makedirs(finalSubsetDir, exist_ok=True)
    targetDf = pd.read_csv(chromosomeTargetPath, sep=' ', header=None)

    for index, row in targetDf.iterrows():
        fileName, chromosomeStr = row[0], row[1]
        if fileName in os.listdir(finalSubsetDir):
            continue

        print(f'Exporting subset by chromosome: {fileName}')
        fullDf = pd.read_csv(os.path.join(significantFullExportDir, fileName), sep=' ', low_memory=False)
        chromosomeList = chromosomeStr.split(',')
        subsetDf = fullDf[fullDf['chr'].astype(str).isin(chromosomeList)]
        subsetDf.to_csv(os.path.join(finalSubsetDir, fileName), index=False)


# --- Run full pipeline ---
extractSuggestiveHits()
extractFullSummariesFromHits()
generateChromosomeTargetList()
exportSignificantChromosomes()
