# --- Setup ---
import pandas as pd
import os

# --- Paths ---
projectRoot = './prj_190_viking_somalogic/'

inputSpreadsheetPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
resultDir = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/significantHitFullExport')

exportDir = os.path.join(resultDir, 'MR/1somamerSNPlists/singleSNP')
multipleSNPexportDir = os.path.join(resultDir, 'MR/1somamerSNPlists/multipleSNP')
ABOHLADir = os.path.join(resultDir, 'MR/1somamerSNPlists/nonanalyzedSNP/HBOHLA')

# --- Load input ---
df = pd.read_excel(inputSpreadsheetPath)

# --- Directory creation ---
def createDirs():
    requiredDirs = [
        os.path.join(resultDir, 'MR/1somamerSNPlists'),
        os.path.join(resultDir, 'MR/1somamerSNPlists/nonanalyzedSNP/'),
        exportDir,
        multipleSNPexportDir,
        ABOHLADir
    ]
    for directory in requiredDirs:
        if not os.path.isdir(directory):
            os.makedirs(directory, exist_ok=True)

# --- SNP Extraction ---
def extractSNP(somamerID, SNP, num):
    filePath = os.path.join(resultDir, somamerID + '.csv')
    if not os.path.isfile(filePath):
        print(somamerID, 'not found!')
        return

    somamerDf = pd.read_csv(filePath, sep=' ')
    row = somamerDf.loc[somamerDf['rsid'] == SNP]
    row.to_csv(os.path.join(exportDir, f'{num}_{somamerID}.csv'), index=False)

# --- Combine Somamers with multiple SNPs ---
def combineSomamers():
    dirFiles = os.listdir(exportDir)
    somamers = []
    for fileName in dirFiles:
        if fileName != 'multipleSNP':
            somamers.append(fileName.split('_')[1])

    seen = set()
    duplicateSomamers = [x for x in somamers if x in seen or seen.add(x)]

    for somamerID in duplicateSomamers:
        result = ''
        for fileName in dirFiles:
            if fileName.endswith('.csv') and fileName.split('_')[1] == somamerID:
                filePath = os.path.join(exportDir, fileName)
                file = pd.read_csv(filePath)
                if len(result) == 0:
                    result = file
                else:
                    result = pd.concat([result, file], ignore_index=True)
        outputPath = os.path.join(multipleSNPexportDir, somamerID)
        result.to_csv(outputPath, index=False)

# --- Remove now-merged single SNP files ---
def removeMultipleSNPFiles():
    toRemove = os.listdir(multipleSNPexportDir)
    for fileName in os.listdir(exportDir):
        if fileName.endswith('.csv') and fileName.split('_')[1] in toRemove:
            os.remove(os.path.join(exportDir, fileName))

# --- Separate ABO and HLA regions ---
def separateABOHLA():
    def checkABOHLA(df):
        keep = []
        remove = []
        for num, chrom in enumerate(df['chr']):
            pos = df.loc[num, 'pos']
            if chrom == 6 and 29645000 < pos < 33365000:
                remove.append(df.loc[num, :])
                continue
            elif chrom == 9 and 136131052 < pos < 136150605:
                remove.append(df.loc[num, :])
                continue
            keep.append(df.loc[num, :])
        return pd.DataFrame(keep, index=None), pd.DataFrame(remove, index=None)

    def moveFiles(directory):
        for fileName in os.listdir(directory):
            filePath = os.path.join(directory, fileName)
            file = pd.read_csv(filePath)
            keepDf, removeDf = checkABOHLA(file)
            if len(keepDf) > 0:
                keepDf.to_csv(filePath, index=False)
            else:
                os.remove(filePath)
            if len(removeDf) > 0:
                removeDf.to_csv(os.path.join(ABOHLADir, fileName), index=False)

    for directory in [exportDir, multipleSNPexportDir]:
        moveFiles(directory)

# --- Main Execution ---
createDirs()

for num, somamerID in enumerate(df['somamerID']):
    print(df.loc[num, 'cisTrans'], df.loc[num, 'novel'])
    if df.loc[num, 'cisTrans'] == 'Cis' and df.loc[num, 'novel'] is True:
        extractSNP(somamerID, df.loc[num, 'SNP'], num)

combineSomamers()
removeMultipleSNPFiles()
separateABOHLA()
