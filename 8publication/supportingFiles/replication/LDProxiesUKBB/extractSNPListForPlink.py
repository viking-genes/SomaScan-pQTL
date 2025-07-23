import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File Paths ---
supplementaryInputPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
snpListOutputDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/LDProxiesUKBB/SNPLists')
ukbbBasePath = './genotypes/10k_unrelated_white_british_reference'

# --- Create Output Directory If It Doesn't Exist ---
os.makedirs(snpListOutputDirectory, exist_ok=True)

# --- Load Supplementary SNP Data ---
supplementaryDf = pd.read_excel(supplementaryInputPath)

# --- Extract SNP Lists by Chromosome ---
for chromosome in range(1, 23):
    snpList = []
    print(f'Processing chromosome {chromosome}...')

    # Load UKBB reference RSID list for the chromosome
    bimFilePath = os.path.join(ukbbBasePath, f'ukbb_chr{chromosome}_10000_random_unrelated_white_british.bim')
    ukbbDf = pd.read_csv(bimFilePath, sep='\t', header=None)
    ukbbDf.columns = ['Chromosome', 'rsid', 'GeneticDistance', 'Position', 'EffectAllele', 'OtherAllele']
    ukbbRsids = set(ukbbDf['rsid'])

    # Filter supplementary data by chromosome
    chromDf = supplementaryDf[supplementaryDf['Chromosome'] == chromosome].reset_index(drop=True)

    # Collect SNPs that are present in the UKBB reference set
    for _, row in chromDf.iterrows():
        if row['SNP'] in ukbbRsids:
            snpList.append(row['SNP'])

    # Remove duplicates
    uniqueSnps = sorted(set(snpList))

    # Save to file
    outputFilePath = os.path.join(snpListOutputDirectory, f'snpListByChromosome{chromosome}.txt')
    with open(outputFilePath, 'w') as outFile:
        for snp in uniqueSnps:
            outFile.write(snp + '\n')
