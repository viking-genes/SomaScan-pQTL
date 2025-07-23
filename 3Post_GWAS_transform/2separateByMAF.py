import pandas as pd
import os

# Define project root and input/output paths
projectRoot = './prj_190_viking_somalogic/'
inputDirectory = os.path.join(projectRoot, 'GWAS/GWAShits/')
mafThresholds = [0.05] 

# Function to filter GWAS results by MAF threshold
def filterByMinorAlleleFrequency(mafThreshold, inputDirectory, outputDirectory):
    os.makedirs(outputDirectory, exist_ok=True)

    for fileName in os.listdir(inputDirectory):
        if not fileName.endswith('.csv'):
            continue

        inputFilePath = os.path.join(inputDirectory, fileName)
        df = pd.read_csv(inputFilePath)

        # Filter rows by MAF (freq1 must be between threshold and 1-threshold)
        filteredDf = df[(df['freq1'] > mafThreshold) & (df['freq1'] < 1 - mafThreshold)]

        # Drop unwanted columns
        filteredDf.drop(
            columns=[col for col in filteredDf.columns if 'unnamed' in col.lower() or 'protein' in col.lower()],
            inplace=True,
            errors='ignore'
        )

        # Write filtered data to output directory
        outputFilePath = os.path.join(outputDirectory, fileName)
        filteredDf.to_csv(outputFilePath, index=False)

# Loop over each threshold and apply filtering
for maf in mafThresholds:
    outputDir = os.path.join(projectRoot, f'GWAS/GWAShitsMAF{maf}/')
    filterByMinorAlleleFrequency(maf, inputDirectory, outputDir)
