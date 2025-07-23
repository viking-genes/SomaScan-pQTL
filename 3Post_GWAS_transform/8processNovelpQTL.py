import pandas as pd
import numpy as np
import math
import requests
import io

# === File Paths ===
novelPQTLPath = './prj_190_viking_somalogic/Post_GWAS_transform/novelpQTLs.xlsx'
outputLDproxyPath = './prj_190_viking_somalogic/Post_GWAS_transform/LDproxy.xlsx'
outputFinalPQTLPath = './prj_190_viking_somalogic/Post_GWAS_transform/novelpQTLs_qc_annotated.xlsx'
qcFlagSourcePath = './prj_190_viking_somalogic/GWAS/viking1_proteomics_somalogic7k_2021_flagged.csv'

# === Parameters ===
ldlinkToken = 'yourTokenHere'
ldlinkPop = 'CEU+TSI+GBR+IBS'
ldlinkR2Threshold = 0.8
ldlinkWindow = 1_000_000


# === Functions ===

def extractQCAnnotationTable(rawQCDataFrame: pd.DataFrame) -> pd.DataFrame:
    """Extract SeqId → Flag QC mapping table from SomaLogic header."""
    headerRows = rawQCDataFrame.iloc[[0, 22]]
    transposed = headerRows.transpose()
    transposed.columns = transposed.iloc[0]
    qcTable = transposed[1:].reset_index(drop=True)
    return qcTable


def addQCFlags(pqtlDf: pd.DataFrame, qcReferenceDf: pd.DataFrame) -> pd.DataFrame:
    """Add SomaLogic QC flag annotations to pQTL dataframe."""
    pqtlDf['FlagQC'] = ''
    for index, somamerID in enumerate(pqtlDf['somamerID']):
        match = qcReferenceDf[qcReferenceDf['SeqId'] == somamerID]
        if len(match) != 1:
            print(f"Warning: somamerID {somamerID} not uniquely found in QC table.")
            continue

        try:
            flag = match.iloc[0, 1]
            pqtlDf.at[index, 'FlagQC'] = 'ok' if pd.isna(flag) else flag
        except Exception:
            pqtlDf.at[index, 'FlagQC'] = 'ok'

    return pqtlDf


def queryLDproxy(rsID: str) -> pd.DataFrame:
    """Query LDproxy API for a given rsID."""
    params = {
        'var': rsID,
        'pop': ldlinkPop,
        'r2_d': 'r2',
        'window': str(ldlinkWindow),
        'genome_build': 'grch37',
        'token': ldlinkToken
    }

    response = requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy', params=params)
    response.raise_for_status()
    data = response.content
    return pd.read_csv(io.StringIO(data.decode('utf-8')), sep='\t')


def pruneLDproxyResults(ldproxyDf: pd.DataFrame, r2Threshold: float) -> pd.DataFrame:
    """Filter LDproxy results by R² threshold and keep proxy SNPs only."""
    ldproxyDf['targetSNP'] = ldproxyDf.at[0, 'RS_Number']
    return ldproxyDf[ldproxyDf['R2'] > r2Threshold].copy()


def batchLDproxyAnalysis(pqtlDf: pd.DataFrame, outputPath: str) -> None:
    """Run LDproxy + pruning across all pQTLs."""
    allResults = []

    for count, rsID in enumerate(pqtlDf['SNP'], start=1):
        print(f"Querying LDproxy for SNP {count}/{len(pqtlDf)}: {rsID}")
        try:
            proxyDf = queryLDproxy(rsID)
            proxyDf = pruneLDproxyResults(proxyDf, ldlinkR2Threshold)
            allResults.append(proxyDf)
        except Exception as e:
            print(f"LDproxy failed for {rsID}: {e}")

    if allResults:
        combined = pd.concat(allResults, ignore_index=True)
        combined.to_excel(outputPath, index=False)
    else:
        print("No valid LDproxy results returned.")


# === Main Execution ===

# Load data
novelDf = pd.read_excel(novelPQTLPath)
qcRawDf = pd.read_csv(qcFlagSourcePath, low_memory=False)

# Export SNP list for SNPNexus
exportSNPsForSNPNexus(novelDf, outputSNPNexusPath)

# Add QC flags from SomaLogic
qcTable = extractQCAnnotationTable(qcRawDf)
novelDf = addQCFlags(novelDf, qcTable)

# Save QC-annotated novel pQTLs
novelDf.to_excel(outputFinalPQTLPath, index=False)

# Run LDproxy & prune by R²
batchLDproxyAnalysis(novelDf, outputLDproxyPath)
