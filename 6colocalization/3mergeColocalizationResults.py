import pandas as pd
import os

# Set threshold for colocalization posterior probability (H4)
hypothesis4Threshold = 0.8

# Define project root and relative paths
projectRoot = './prj_190_viking_somalogic/'
resultsDirectory = os.path.join(projectRoot, 'Scripts/colocalization/results/')
inputDataPath = os.path.join(projectRoot, 'Scripts/colocalization/1mergedTwoSampleMRResults.tsv')
outputExcelPath = os.path.join(projectRoot, 'Scripts/colocalization/3mergeColocalizationResults.xlsx')

# Load merged TwoSampleMR results
mergedResultsDataFrame = pd.read_csv(inputDataPath, sep='\t')

# Prepare a working copy and clean 'outcome' column to keep only study name
mergedResultsDataFrame['colocalizationH4'] = ''
resultsDataFrameCopy = mergedResultsDataFrame.copy()
resultsDataFrameCopy['outcome'] = resultsDataFrameCopy['outcome'].apply(lambda x: x.split(' || id:')[-1])

# Loop through all colocalization result files
for resultFileName in os.listdir(resultsDirectory):

    resultFilePath = os.path.join(resultsDirectory, resultFileName)

    # Read 7th line from result file and extract H4 value
    with open(resultFilePath, 'r') as resultFile:
        lines = resultFile.readlines()
        h4Value = float(lines[6].split('\t')[1].strip())

    # Extract protein and study name from filename
    fileNameParts = resultFileName.split('_')
    proteinName = fileNameParts[0]
    studyName = fileNameParts[1].replace('.txt', '')

    # Filter data for matching exposure (protein name)
    proteinMatches = mergedResultsDataFrame[mergedResultsDataFrame['exposure'] == proteinName]

    # Find index where study matches
    matchingIndexes = resultsDataFrameCopy.index[resultsDataFrameCopy['outcome'] == studyName]
    validIndexes = [idx for idx in matchingIndexes if idx in proteinMatches.index]

    # Handle mismatches and duplicates
    if len(validIndexes) == 0:
        print('No valid study match:', resultFileName, '| Study:', studyName)
        continue
    if len(validIndexes) > 1:
        print('Error: multiple valid matches for', studyName)
        print(proteinMatches)
        break

    # Write H4 value into main dataframe
    mergedResultsDataFrame.at[validIndexes[0], 'colocalizationH4'] = h4Value

# Export final DataFrame to Excel
mergedResultsDataFrame.to_excel(outputExcelPath, index=False)
