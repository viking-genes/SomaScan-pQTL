import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Input and Output Paths ---
pathToLDProxy = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/LDproxy')
outputFile = os.path.join(pathToLDProxy, 'LDproxyRSID.txt')

# --- Initialize Result List ---
result = []

# --- Iterate Through .xlsx Files and Collect RSIDs ---
for fileName in os.listdir(pathToLDProxy):
    if not fileName.endswith('.xlsx'):
        continue

    filePath = os.path.join(pathToLDProxy, fileName)
    df = pd.read_excel(filePath)

    rsIDs = df['RS_Number'].tolist()
    result.extend(rsIDs)

# --- Write Combined RSIDs to Output File ---
with open(outputFile, 'w') as f:
    for rsid in result:
        f.write(f"{rsid}\n")
