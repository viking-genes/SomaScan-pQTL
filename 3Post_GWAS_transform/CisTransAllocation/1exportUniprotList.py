import numpy as np
import pandas as pd
import os

def exportUniprotList():
    # --- Setup paths ---
    projectRoot = './prj_190_viking_somalogic/'
    inputPath = os.path.join(projectRoot, 'GWAS/viking1_proteomics_somalogic7k_2021_flagged.csv')
    outputPath = os.path.join(projectRoot, 'Post_GWAS_transform/CisTransAllocation/exportUniprotList.txt')

    # --- Load data ---
    df = pd.read_csv(inputPath, low_memory=False)

    # --- Extract UniProt IDs from row 5 (6th row, where protein mappings are stored) ---
    uniprotList = df.iloc[5, :].tolist()

    # --- Remove unwanted values ---
    unwanted = {5, 'UniProt', np.nan}
    filtered = [x for x in uniprotList if x not in unwanted and isinstance(x, str)]

    # --- Get unique values and export ---
    uniqueUniprots = sorted(set(filtered))
    pd.DataFrame(uniqueUniprots).to_csv(outputPath, header=False, index=False)

    print(f'✅ Exported {len(uniqueUniprots)} unique UniProt IDs to: {outputPath}')

# --- Run ---
exportUniprotList()
