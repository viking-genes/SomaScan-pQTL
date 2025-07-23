import pandas as pd
import os

# --- File paths ---
projectRoot = './prj_190_viking_somalogic/'
pqtlPath = os.path.join(projectRoot, 'Post_GWAS_transform/novelpQTLs.xlsx')
nexusPath = os.path.join(projectRoot, 'Post_GWAS_transform/SNPNexusAnnotation.txt')

# --- Load data ---
pqtlDf = pd.read_excel(pqtlPath)
nexusDf = pd.read_csv(nexusPath, sep='\t')

# --- Annotate pQTLs with SNPnexus gene information ---
def annotateWithSnpnexus(pqtlDf, nexusDf):
    pqtlDf['Overlapped Gene'] = ''
    pqtlDf['Nearest Upstream Gene'] = 'NA'
    pqtlDf['Nearest Downstream Gene'] = 'NA'

    for idx, snpId in pqtlDf['SNP'].items():
        matchedRows = nexusDf[nexusDf['Variation ID'] == snpId]

        if not matchedRows.empty:
            overlappedGenes = matchedRows['Overlapped Gene'].dropna().tolist()
            overlappedGeneString = '|'.join(overlappedGenes) if overlappedGenes else 'None'
            pqtlDf.at[idx, 'Overlapped Gene'] = overlappedGeneString

            if overlappedGeneString == 'None':
                upstream = matchedRows['Nearest Upstream Gene'].dropna().tolist()
                downstream = matchedRows['Nearest Downstream Gene'].dropna().tolist()

                if upstream:
                    pqtlDf.at[idx, 'Nearest Upstream Gene'] = upstream[0]
                if downstream:
                    pqtlDf.at[idx, 'Nearest Downstream Gene'] = downstream[0]

    return pqtlDf

# --- Run annotation ---
pqtlDf = annotateWithSnpnexus(pqtlDf, nexusDf)

# --- Save updated file ---
pqtlDf.to_excel(pqtlPath, index=False)
print(f"✅ Annotated pQTLs saved to: {pqtlPath}")
