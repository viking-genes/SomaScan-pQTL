import pandas as pd
import os

# --- File paths ---
projectRoot = './prj_190_viking_somalogic/'
pqtlFilePath = os.path.join(projectRoot, 'Post_GWAS_transform/novelpQTLs.xlsx')
uniprotToHugoPath = os.path.join(projectRoot, 'Post_GWAS_transform/uniprotHugoMap.xlsx')
snpnexusOutputPath = os.path.join(projectRoot, 'Post_GWAS_transform/SNPNexusImport.txt')

# --- Load data ---
df = pd.read_excel(pqtlFilePath)
dfHugo = pd.read_excel(uniprotToHugoPath)

# --- Add 'HUGO' column if missing ---
if 'HUGO' not in df.columns:
    columns = df.columns.tolist()
    columns.insert(14, 'HUGO')
    df = df.reindex(columns=columns)

# --- Map UniProt to HUGO ---
for idx, uniprot in df['Uniprot'].items():
    if pd.isna(uniprot):
        continue

    hugoSymbols = []
    for token in uniprot.split('|'):
        matches = dfHugo[dfHugo['From'] == token]['To'].tolist()
        if matches:
            hugoSymbols.append(matches[0])
        else:
            print(f'Warning! UniProt ID {token} not found in mapping table.')

    if hugoSymbols:
        df.at[idx, 'HUGO'] = '|'.join(hugoSymbols)

# --- Save annotated novel pQTLs back to Excel ---
df.to_excel(pqtlFilePath, index=False)
print(f"✅ Updated pQTLs with HUGO symbols saved to: {pqtlFilePath}")

# --- Prepare SNPnexus input ---
snpnexusInput = pd.DataFrame({
    0: ['dbsnp'] * len(df),
    1: df['SNP']
})
snpnexusInput.to_csv(snpnexusOutputPath, sep='\t', index=False, header=False)
print(f"✅ SNPnexus input saved to: {snpnexusOutputPath}")
