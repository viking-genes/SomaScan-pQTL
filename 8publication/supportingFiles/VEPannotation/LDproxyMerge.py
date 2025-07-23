import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File Paths ---
LDproxyDf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/LDproxy/LDproxy.xlsx')
VEPDf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/VEP/VEPoutput230705.txt')
outputFile = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/LDproxyVEP.xlsx')

# --- Load Data ---
LDproxyDf = pd.read_excel(LDproxyDf)
VEPDf = pd.read_csv(VEPDf, sep='\t')

# --- Annotate LDproxy DataFrame with VEP Consequence ---
LDproxyDf['VEP'] = ''

for i in range(len(LDproxyDf)):
    rsid = LDproxyDf.loc[i, 'RS_Number']
    indexes = VEPDf.index[VEPDf['#Uploaded_variation'] == rsid].tolist()

    if indexes:
        LDproxyDf.loc[i, 'VEP'] = VEPDf.loc[indexes[0], 'Consequence']
    else:
        LDproxyDf.loc[i, 'VEP'] = 'NA'

# --- Save Merged DataFrame ---
LDproxyDf.to_excel(outputFile, index=False)
