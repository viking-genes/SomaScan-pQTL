# --- Setup ---
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Paths and parameters ---
projectRoot = './prj_190_viking_somalogic/'

resultDir = os.path.join(projectRoot, 'Scripts/twoSampleMr/2resultsMR/')
outputDir = os.path.join(projectRoot, 'Scripts/twoSampleMr/3significantResultsMR/')
annotationPath = os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/openGWASChatGPTAnnotated.xlsx')
fdrPlotPath = os.path.join(projectRoot, 'Scripts/twoSampleMr/3FDRplot.png')

medicallyRelevantAnnotation = pd.read_excel(annotationPath)

FDRthreshold = 0.01
methodList = ['wald_ratio', 'inverse_variance_weighte=d__fixed_effects_']

# --- Utilities ---
def getColumnNumber(path):
    with open(path, 'r') as file:
        col_count = [len(line.split(',')) for line in file.readlines()]
    return list(range(max(col_count)))

def fixDfDelimiter(df, columnNames):
    expectedMethods = [
        'wald_ratio',
        'inverse_variance_weighted',
        'inverse_variance_weighted__fixed_effects_',
        'maximum_likelihood',
        'weighted_mode',
        'weighted_median',
        'mr_egger'
    ]
    for num, method in enumerate(df[2][1:]):
        while method not in expectedMethods:
            df.at[num + 1, 1] = str(df.loc[num + 1, 1]) + '<><>' + str(df.loc[num + 1, 2])
            for j in range(2, len(columnNames) - 1):
                df.at[num + 1, j] = df.loc[num + 1, j + 1]
            method = df.loc[num + 1, 2]
    df.drop(columns=[x for x in range(7, len(columnNames))], inplace=True)
    df = df.rename(columns=df.iloc[0]).drop(index=0).reset_index(drop=True)
    return df

def removeCohorts(df):
    exclude = [
        'Asthma (childhood onset) || id:ebi-a-GCST007800',
        'Asthma (adult onset) || id:ebi-a-GCST007799',
        'Breast cancer || id:ebi-a-GCST007236',
        'Schizophrenia || id:ieu-b-5070'
    ]
    return df[~df['outcome'].isin(exclude)].reset_index(drop=True)

def filterMedicallyRelevantResults(df):
    outcomeIds = [out.split(' || id:')[1] for out in df['outcome']]
    keepIndices = []
    for idx, outcomeId in enumerate(outcomeIds):
        annotationIdx = medicallyRelevantAnnotation[medicallyRelevantAnnotation['id.outcome'] == outcomeId].index.tolist()
        if not annotationIdx:
            continue
        if len(annotationIdx) == 1 and medicallyRelevantAnnotation.loc[annotationIdx[0], 'Clinically relevant'] == 'Yes':
            keepIndices.append(idx)
        elif len(annotationIdx) > 1:
            print('Error in the annotation database! Multiple matches found for', outcomeId)
    return df.iloc[keepIndices].reset_index(drop=True)

def removeUnwantedMethods(df, methods):
    return df[df['method'].isin(methods)].reset_index(drop=True)

def extractSignificant(df, pThreshold):
    return df[df['p'].astype(float) < pThreshold].reset_index(drop=True)

# --- FDR Calculations ---
def findFDRpValue(pValues, numberOfTests, fdrThreshold):
    for i in range(numberOfTests):
        if fdrThreshold * (i + 1) / numberOfTests < pValues[i]:
            print('Last position to pass FDR test', i - 1, 'pValue:', pValues[i - 1])
            return i, pValues[i - 1]

def plotFDRgraph(pValues, numberOfTests, numberOfSignificantHits):
    plotThreshold = int(numberOfSignificantHits * 1.5)
    x = list(range(plotThreshold))
    plt.plot(x, pValues[:plotThreshold])
    plt.plot(x, [FDRthreshold * (i + 1) / numberOfTests for i in x])
    plt.xlabel('Rank')
    plt.ylabel('p-value')
    plt.title('FDR Plot')
    plt.savefig(fdrPlotPath)
    plt.close()

# --- Scan and extract p-values from MR result directory ---
def extractPvaluesFromDir(MRresultDir, methodList):
    pValues = []
    totalTests = 0

    for fileName in os.listdir(MRresultDir):
        if not fileName.endswith('.csv'):
            continue

        filePath = os.path.join(MRresultDir, fileName)
        columnNames = getColumnNumber(filePath)
        df = pd.read_csv(filePath, header=None, names=columnNames)
        df = fixDfDelimiter(df, columnNames)
        df = removeCohorts(df)
        df = filterMedicallyRelevantResults(df)
        df = removeUnwantedMethods(df, methodList)

        totalTests += len(df)
        pValues.extend(df['p'].astype(float))

    pValues = np.sort(pValues)
    print('Number of tests:', totalTests)
    return pValues, totalTests

# --- Main Execution ---
pValues, numberOfTests = extractPvaluesFromDir(resultDir, methodList)
numberOfSignificantHits, thresholdPValue = findFDRpValue(pValues, numberOfTests, FDRthreshold)
plotFDRgraph(pValues, numberOfTests, numberOfSignificantHits)

# --- Save significant MR results ---
os.makedirs(outputDir, exist_ok=True)

for fileName in os.listdir(resultDir):
    if not fileName.endswith('.csv'):
        continue

    filePath = os.path.join(resultDir, fileName)
    columnNames = getColumnNumber(filePath)
    df = pd.read_csv(filePath, header=None, names=columnNames)
    df = fixDfDelimiter(df, columnNames)
    df = removeCohorts(df)
    df = removeUnwantedMethods(df, methodList)
    df = extractSignificant(df, thresholdPValue)

    if not df.empty:
        outputFileName = fileName.split('_')[-1]
        outputPath = os.path.join(outputDir, outputFileName)
        df.to_csv(outputPath, index=False)
