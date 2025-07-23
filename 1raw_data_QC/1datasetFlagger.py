# --- Setup ---
import os
import numpy as np
import pandas as pd

# Define project root and file paths
projectRoot = './prj_190_viking_somalogic/'
rawDataPath = os.path.join(
    projectRoot,
    'raw_data/viking/phenotypes/omics/2021_somalogic_proteomics/UED_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat'
)

# Parameters
allowedBlankOverlap = 0.05
metadataRowCount = 61
metadataColCount = 33

# --- Step 1: Load and clean raw data ---
# Load raw data (7630 columns, no headers)
rawData = pd.read_csv(rawDataPath, delimiter='\t', names=range(7630), low_memory=False)

# Drop metadata rows
rawData.drop(index=range(metadataRowCount), inplace=True)

# Extract sample IDs and protein identifiers
sampleIds = rawData[6][22:].reset_index(drop=True)
proteinIds = rawData.loc[metadataRowCount].reset_index(drop=True)

# Drop metadata columns
rawData.drop(columns=range(metadataColCount), inplace=True)
rawData.reset_index(drop=True, inplace=True)
rawData.columns = range(rawData.shape[1])

# Replace sample ID column values with clean ones
for pos in range(len(sampleIds)):
    rawData.at[22 + pos, 0] = sampleIds[pos]

# --- Step 2: Extract blank samples ---
blankRows = []
for index, sampleId in rawData[0].items():
    if index <= 22:
        continue
    if not str(sampleId).startswith('VIKI'):
        if str(sampleId) == '210230':
            blankRows.append(rawData.loc[index])
        rawData.drop(index, inplace=True)

blankData = pd.DataFrame(blankRows).astype(float)

# Add a new row for flags
rawData.at[22, 0] = 'Wilson flag'

# --- Step 3: Define blank overlap check function ---
def isBlankOverlap(blankValues, sampleValues, threshold):
    maxBlank = max(blankValues)
    overlapCount = sum(float(value) <= maxBlank for value in sampleValues)
    return (overlapCount / len(sampleValues)) >= threshold

# --- Step 4: Initialize flag collections ---
noProteinFlags = []
nonHumanFlags = []
somalogicQCFlags = []
blankOverlapFlags = []

# --- Step 5: Flag aptamers ---
for columnIndex in rawData.columns[1:]:
    flags = []

    # Flag: no protein target
    if rawData[columnIndex][3] == 'No protein':
        noProteinFlags.append(rawData[columnIndex])
        flags.append('No protein')

    # Flag: nonhuman target
    if rawData[columnIndex][8] != 'Human':
        nonHumanFlags.append(rawData[columnIndex])
        flags.append('Nonhuman')

    # Flag: SomaLogic internal QC flag
    if rawData[columnIndex][15] == 'FLAG':
        somalogicQCFlags.append(rawData[columnIndex])
        flags.append('Somalogic QC')

    # Flag: blank sample overlap
    sampleValues = rawData.loc[23:, columnIndex]
    if isBlankOverlap(blankData[columnIndex], sampleValues, allowedBlankOverlap):
        blankOverlapFlags.append(rawData[columnIndex])
        flags.append('Blank overlap')

    # Append flags to the flag row if any
    if flags:
        currentFlags = rawData[columnIndex][22]
        rawData.at[22, columnIndex] = ' '.join([str(currentFlags)] + flags) if pd.notnull(currentFlags) else ' '.join(flags)

# --- Step 6: Export of flagged data ---
flaggedDataPath = os.path.join(
    projectRoot,
    'processing/viking/phenotypes/omics/2021_somalogic_proteomics/viking1_proteomics_somalogic7k_2021_flagged.csv'
)
rawData.to_csv(flaggedDataPath, index=False, header=False)

print(f"Number of aptamers with blank overlap: {len(blankOverlapFlags)}")
