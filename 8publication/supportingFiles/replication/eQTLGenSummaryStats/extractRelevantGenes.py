import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Output Paths ---
eqtlGenSummaryPath = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/eQTLGenSummaryStats.txt.gz')
eqtlGenFreqPath = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/eQTLGenAlleleFreq.txt.gz')
eqtlGenOutputDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract')

# --- Ensure Output Directory Exists ---
os.makedirs(eqtlGenOutputDirectory, exist_ok=True)

# --- Load Large Summary Files (Requires ~128 GB RAM) ---
eqtlGenSummaryDf = pd.read_csv(eqtlGenSummaryPath, sep='\t', compression='gzip', low_memory=False)
eqtlGenFreqDf = pd.read_csv(eqtlGenFreqPath, sep='\t', compression='gzip', low_memory=False)

# --- Define Gene Dictionary (CYP2C19 Not Present in This Dataset) ---
geneDictionary = {
    'ENSG00000156795': 'NTAQ1',
    'ENSG00000145781': 'COMMD10',
    'ENSG00000110987': 'BCL7A',
    'ENSG00000196290': 'NIF3L1',
    'ENSG00000173621': 'LRFN4',
    'ENSG00000121064': 'SCPEP1',
    'ENSG00000257594': 'GALNT4',
    'ENSG00000109956': 'B3GAT1',
    'ENSG00000087884': 'AAMDC',
    'ENSG00000133574': 'GIMAP4',
    'ENSG00000101000': 'PROCR',
    'ENSG00000062524': 'LTK',
    'ENSG00000115232': 'ITGA4|ITGB1'
}

# --- Extract and Annotate Gene-Specific eQTL Data ---
for ensemblId, geneSymbol in geneDictionary.items():
    geneDf = eqtlGenSummaryDf[eqtlGenSummaryDf['Gene'] == ensemblId].copy()
    geneDf = geneDf.sort_values(by='SNPPos')

    # Subset allele frequency data for relevant SNPs
    alleleFreqSubset = eqtlGenFreqDf[eqtlGenFreqDf['SNP'].isin(geneDf['SNP'])]
    
    # Merge and annotate frequency column
    geneDf = geneDf.merge(alleleFreqSubset[['SNP', 'AlleleB_all']], on='SNP', how='left')
    geneDf = geneDf.rename(columns={'AlleleB_all': 'AssessedAllele_freq'})

    # Save output
    outputPath = os.path.join(eqtlGenOutputDirectory, f"{geneSymbol}.tsv")
    geneDf.to_csv(outputPath, sep='\t', index=False)
