import requests
import pandas as pd
import io
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File Paths ---
inputFile = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
outputFile = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/LDproxy/LDproxy.xlsx')

# --- Load Supplementary SNPs ---
df = pd.read_excel(inputFile)

# --- Define Function to Query LDproxy ---
def LDproxyGet(rsID):
    params = {
        'var': str(rsID),
        'pop': 'CEU+TSI+GBR+IBS',
        'r2_d': 'r2',
        'window': '100000',
        'genome_build': 'grch37',
        'token': 'yourTokenHere'  
    }

    response = requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy', params=params, verify=True)
    data = response.content
    return pd.read_csv(io.StringIO(data.decode('utf-8')), sep='\t')

# --- Filter LD Results Based on R2 > 0.8 ---
def LDproxyPrune(df):
    result = []
    df['targetSNP'] = df.at[0, 'RS_Number']
    
    for _, row in df.iterrows():
        if row['R2'] > 0.8:
            result.append(row)
        else:
            break  # Assumes results are sorted by R2 descending

    return pd.DataFrame(result)

# --- Initialize ---
result = []
count = 0
excludeSNP = ['rs28895246', '9:140008750_G_C', 'rs77553517', 'rs2961544', 'rs78967641', 'rs9501393']

# --- Process Each Unique SNP ---
uniqueSnps = df['SNP'].unique()

for snp in uniqueSnps:
    count += 1
    if snp in excludeSNP:
        continue

    print(f'Processing {count}/{len(uniqueSnps)}: {snp}')
    
    try:
        ldProxy = LDproxyGet(snp)
        ldProxy = LDproxyPrune(ldProxy)

        result.append(ldProxy)
    except Exception as e:
        print(f'Failed to process {snp}: {e}')

# --- Concatenate and Export ---
if result:
    finalResult = pd.concat(result, ignore_index=True)
    os.makedirs(os.path.dirname(outputFile), exist_ok=True)
    finalResult.to_excel(outputFile, index=False)
