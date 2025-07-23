import pandas as pd
import re
import os

# Set project root directory
projectRoot = './prj_190_viking_somalogic/'

# Define data paths
phenotypeDataPath = os.path.join(projectRoot, 'Data/combinedCovariateDf.tsv')
proteinLevelsDataPath = os.path.join(projectRoot, 'Data_transformation/rank_transformed_covariate_unadjusted.csv')
outputPhenotypePath = os.path.join(projectRoot, 'GWAS_preparation/st02_01_qcd_phenotypes.txt')
gwasBaseDirectory = os.path.join(projectRoot, 'GWAS/Full/')

# Load input data
phenotypeDataFrame = pd.read_csv(phenotypeDataPath, delimiter='\t')
proteinLevelsDataFrame = pd.read_csv(proteinLevelsDataPath)

# Remove overlapping columns
def dropIdenticalColumns(baseDataFrame, columnsToCheck):
    for columnName in columnsToCheck:
        if columnName in baseDataFrame:
            baseDataFrame.drop(columns=columnName, inplace=True)
    return baseDataFrame

# Normalize column names for compatibility
def manageColumnNames(dataFrame, prefix, separator):
    updatedColumns = [prefix + re.sub(r'[-\s]', separator, colName) for colName in dataFrame.columns]
    dataFrame.columns = updatedColumns
    return dataFrame

# Apply cleaning functions
proteinLevelsDataFrame = dropIdenticalColumns(proteinLevelsDataFrame, phenotypeDataFrame.columns)
proteinLevelsDataFrame = manageColumnNames(proteinLevelsDataFrame, 'str', '__')
phenotypeDataFrame = manageColumnNames(phenotypeDataFrame, '', '_')

# Combine phenotype and protein data into a single DataFrame
combinedDataFrame = pd.concat([phenotypeDataFrame, proteinLevelsDataFrame], axis=1, join="inner")

# Save final phenotype table
os.makedirs(os.path.dirname(outputPhenotypePath), exist_ok=True)
combinedDataFrame.to_csv(outputPhenotypePath, index=False, sep='\t')

# Write GWAS analysis plan in chunks
def createAnalysisChunks(dataFrame, outputDirectory, chunkStartIndex, chunkSize):
    # Make control file directory
    controlFileDir = os.path.join(outputDirectory, 'control_files')
    os.makedirs(controlFileDir, exist_ok=True)

    # Define path for the YAML analysis plan
    analysisPlanPath = os.path.join(controlFileDir, 'analysis_plan.yml')
    if os.path.exists(analysisPlanPath):
        os.remove(analysisPlanPath)

    # Write YAML plan file
    with open(analysisPlanPath, 'w') as planFile:
        # Write standard covariates block
        planFile.write(
            'standard_covariates:\n'
            '  formula: sex + age + Plate_number + Plate_row + Plate_column + Days_since_venepuncture + pc1 + pc2 + pc3\n'
            '  analysis_type: NA\n'
            '  residual_transform: none\n'
            '  fit_mixed_model: NA\n'
        )
        # Write chunk of phenotype formulas
        chunkColumns = dataFrame.columns[chunkStartIndex:chunkStartIndex + chunkSize]
        for phenotype in chunkColumns:
            planFile.write(
                f'\n{phenotype}:\n'
                f'  formula: {phenotype} ~ standard_covariates\n'
                f'  analysis_type: quantitative_trait\n'
                f'  residual_transform: none\n'
                f'  fit_mixed_model: TRUE\n'
            )

# Split the protein dataframe into 18 chunks for GWAS execution to fit all protein measurements in SomaScan v4.1
chunkSize = 422
numberOfChunks = 18

for chunkIndex in range(numberOfChunks):
    startIndex = chunkIndex * chunkSize
    outputDir = os.path.join(gwasBaseDirectory, f'Partial{chunkIndex}')
    createAnalysisChunks(proteinLevelsDataFrame, outputDir, startIndex, chunkSize)
