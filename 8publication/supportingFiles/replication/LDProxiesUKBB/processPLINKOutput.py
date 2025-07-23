import pandas as pd
import os
import re

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Input and Output Paths ---
ldProxyInputDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/LDProxiesUKBB/PLINKOutput')
ldProxyOutputDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/LDProxiesUKBB/ProcessedOutput')

# --- Ensure Output Directory Exists ---
os.makedirs(ldProxyOutputDirectory, exist_ok=True)

# --- Convert .ld Files to .tsv ---
for filename in os.listdir(ldProxyInputDirectory):
    if filename.endswith('.ld'):
        inputFilePath = os.path.join(ldProxyInputDirectory, filename)
        outputFilePath = os.path.join(ldProxyOutputDirectory, f"{filename.split('.')[0]}.tsv")

        with open(inputFilePath, 'r') as file:
            rawText = file.read()

            # Replace multiple spaces with a single tab
            cleanedText = re.sub(r' +', '\t', rawText)
            # Remove trailing tabs
            cleanedText = re.sub(r'\t\n', '\n', cleanedText)
            # Remove leading tab from each line
            cleanedText = re.sub(r'\n\t', '\n', cleanedText)
            # Remove leading tab from the first line if present
            cleanedText = re.sub(r'^\t', '', cleanedText)

        with open(outputFilePath, 'w') as outputFile:
            outputFile.write(cleanedText)
