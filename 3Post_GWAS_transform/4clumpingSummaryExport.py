import os
import pandas as pd
from pathlib import Path
import shutil

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
inputDir = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/onlySignficantHits/')
combinedDir = os.path.join(inputDir, 'combinedResults/')
rawProteinDfPath = os.path.join(projectRoot, 'GWAS/viking1_proteomics_somalogic7k_2021_flagged.csv')
gwasResultDir = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/significantHitFullExport/')
significanceThreshold = 5e-8

# Output directory for final result files
exportDir = os.path.join(os.path.dirname(inputDir), f'clumpResults{significanceThreshold}')
os.makedirs(exportDir, exist_ok=True)

# Files to ignore in directory scans
directoryIgnoreList = ['processed', 'combinedResults', '3extractSentinelSNP.py']

# --- Function: Combine clumped and not-clumped SNPs by chromosome ---
def combineClumpedChromosomes(parentDir, folderName):
    resultsDir = os.path.join(parentDir, folderName, 'clumpResults')
    chrDir = os.path.join(parentDir, folderName, 'separatedByChromosome')
    finalDir = os.path.join(parentDir, 'combinedResults')
    os.makedirs(finalDir, exist_ok=True)

    clumpedFrames = []
    notClumpedFrames = []

    for chrNum in range(1, 23):
        clumpFile = f'{chrNum}.csv.clumped'
        logFile = f'{chrNum}.csv.log'
        chrFile = f'{chrNum}.csv'

        clumpPath = os.path.join(resultsDir, clumpFile)
        logPath = os.path.join(resultsDir, logFile)
        chrPath = os.path.join(chrDir, chrFile)

        if os.path.exists(clumpPath):
            clumpedFrames.append(pd.read_csv(clumpPath, sep=r'\s+'))
        elif os.path.exists(logPath) and os.path.exists(chrPath):
            notClumpedFrames.append(pd.read_csv(chrPath, sep=' '))

    if clumpedFrames:
        combinedClumped = pd.concat(clumpedFrames, ignore_index=True) if len(clumpedFrames) > 1 else clumpedFrames[0]
        combinedClumped.to_csv(os.path.join(finalDir, f'{folderName}.csv'), index=False)

    if notClumpedFrames:
        combinedNotClumped = pd.concat(notClumpedFrames, ignore_index=True) if len(notClumpedFrames) > 1 else notClumpedFrames[0]
        combinedNotClumped.to_csv(os.path.join(finalDir, f'{folderName}notclumped.csv'), index=False)

# --- Function: Extract SNPs below threshold from clumped or not-clumped files ---
def extractBelowThreshold(filePath, fileName, threshold):
    df = pd.read_csv(filePath)
    outFileName = 'unclumpedSNP.csv' if 'notclumped' in fileName else 'clumpedSNP.csv'

    if 'notclumped' in fileName:
        if df.loc[0, 'p'] < threshold:
            writeMode = 'a' if os.path.exists(os.path.join(exportDir, outFileName)) else 'w'
            df.to_csv(os.path.join(exportDir, outFileName), mode=writeMode, index=False, header=(writeMode == 'w'))
    else:
        result = df[df['P'] < threshold].copy()
        result['somamerID'] = fileName.replace('.csv', '')
        if not result.empty:
            writeMode = 'a' if os.path.exists(os.path.join(exportDir, outFileName)) else 'w'
            result.to_csv(os.path.join(exportDir, outFileName), mode=writeMode, index=False, header=(writeMode == 'w'))

# --- Function: Add UniProt column using raw proteomics reference ---
def addUniprotColumn(resultDf, rawDf):
    resultDf['Uniprot'] = ''
    for index, row in resultDf.iterrows():
        somamerID = row['somamerID']
        matchColumn = rawDf.columns[rawDf.isin([somamerID]).any()]
        if not matchColumn.empty:
            resultDf.at[index, 'Uniprot'] = rawDf.loc[5, matchColumn[0]]
    return resultDf

# --- Step 1: Combine clumped and not-clumped SNPs per file ---
for folderName in os.listdir(inputDir):
    folderPath = os.path.join(inputDir, folderName)
    if os.path.isdir(folderPath) and folderName not in directoryIgnoreList:
        combineClumpedChromosomes(inputDir, folderName)

# --- Step 2: Extract significant SNPs below threshold ---
for fileName in os.listdir(combinedDir):
    filePath = os.path.join(combinedDir, fileName)
    extractBelowThreshold(filePath, fileName, significanceThreshold)

# --- Step 3: Add Uniprot column to final outputs ---
rawProteinDf = pd.read_csv(rawProteinDfPath, low_memory=False)

for fileName in os.listdir(exportDir):
    filePath = os.path.join(exportDir, fileName)
    df = pd.read_csv(filePath)
    updatedDf = addUniprotColumn(df, rawProteinDf)
    updatedDf.to_csv(filePath, index=False)

# --- Step 4: Annotate SNPs with effect size and frequency ---
def annotateSnpsWithEffectSize():
    snpDbPath = os.path.join(exportDir, 'clumpedSNP.csv')
    if not os.path.exists(snpDbPath):
        print("No clumped SNP database found.")
        return

    snpDb = pd.read_csv(snpDbPath)
    snpDb['Effect Size'] = ''
    snpDb['Frequency'] = ''

    gwasFiles = {f for f in os.listdir(gwasResultDir) if f.endswith('.csv')}
    grouped = snpDb.groupby('somamerID')

    for idx, (somamerID, group) in enumerate(grouped):
        print(f'Annotating {somamerID} ({idx+1}/{len(grouped)})')
        gwasFileName = f'{somamerID}.csv'
        if gwasFileName in gwasFiles:
            gwasDf = pd.read_csv(os.path.join(gwasResultDir, gwasFileName), sep=' ', low_memory=False)
            for snp, rowIndex in zip(group['SNP'], group.index):
                match = gwasDf[gwasDf['rsid'] == snp]
                if not match.empty:
                    snpDb.at[rowIndex, 'Effect Size'] = match['beta1'].values[0]
                    snpDb.at[rowIndex, 'Frequency'] = match['freq1'].values[0]
        else:
            print(f'Warning: {gwasFileName} not found in GWAS results.')

    snpDb.to_excel(os.path.join(projectRoot, 'Post_GWAS_transform/3-1export.xlsx'), index=False)


annotateSnpsWithEffectSize()
