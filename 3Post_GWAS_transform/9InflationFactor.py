import os
import pandas as pd
import numpy as np
from scipy.stats import chi2

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
gwasDir = os.path.join(projectRoot, 'GWAS/Full/p05_comb_chr/')
outputDir = os.path.join(projectRoot, 'Visualization/InflationFactor/TempFiles/LambdaAll/')
csvExportFile = os.path.join(outputDir, 'GWASinflation.csv')

os.makedirs(outputDir, exist_ok=True)

# --- Function to calculate genomic inflation factor (lambda) ---
def extractInflation(dataFrame, fileName):
    medianPvalue = np.nanmedian(dataFrame['p'])
    chiValue = chi2.ppf(medianPvalue, df=1)
    inflationFactor = round(chiValue / chi2.ppf(0.5, df=1), 3)

    resultDf = pd.DataFrame([[fileName, inflationFactor]], columns=['File', 'Inflation factor'])

    # Append or write header if file does not exist
    resultDf.to_csv(csvExportFile, mode='a', header=not os.path.isfile(csvExportFile), index=False)

# --- Process all files ---
fileList = sorted(os.listdir(gwasDir))
for count, fileName in enumerate(fileList, start=1):
    print(f'Reading {fileName} ({count}/{len(fileList)})')
    filePath = os.path.join(gwasDir, fileName)
    df = pd.read_csv(filePath, sep='\t', compression='gzip', low_memory=False)
    extractInflation(df, fileName)
