import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
clumpedSnpPath = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/clumpResults5e-08/clumpedSNP.csv')
tssPath = os.path.join(projectRoot, 'Scripts/Post_GWAS_transform/CisTransAllocation/Database exports/biomartExportingGrCh37/TSSmartExport.csv')
outputPath = clumpedSnpPath.replace('.csv', '_cisTransGrCh37.csv')

# --- Load Data ---
snpDf = pd.read_csv(clumpedSnpPath)
tssDf = pd.read_csv(tssPath)

# --- Step 1: Add TSS and Chromosome info to SNPs ---
def addTssAnnotations(snpDf, tssDf):
    snpDf['Somamer target TSS'] = ''
    snpDf['Somamer target TSS Chr'] = ''

    def resolveTssAndChr(uniprotId):
        proteins = uniprotId.split('|')
        tssList, chrList = [], []
        for protein in proteins:
            try:
                subDf = tssDf[tssDf['Uniprot'] == protein]
                tssList.extend(subDf['Transcription start site (TSS)'].tolist())
                chrList.extend(subDf['Chromosome/scaffold name'].tolist())
            except ValueError:
                return 'NA', 'NA'
        return '|'.join(map(str, tssList)) if tssList else 'NA', '|'.join(map(str, chrList)) if chrList else 'NA'

    for idx, uniprotId in enumerate(snpDf['Uniprot']):
        if pd.isna(uniprotId):
            continue
        tss, chrom = resolveTssAndChr(str(uniprotId))
        snpDf.at[idx, 'Somamer target TSS'] = tss
        snpDf.at[idx, 'Somamer target TSS Chr'] = chrom

    return snpDf

# --- Step 2: Allocate cis/trans based on position and distance ---
def allocateCisTrans(snpDf):
    snpDf['cisTrans'] = ''
    snpDf['Genetic distance'] = ''

    for idx, row in snpDf.iterrows():
        snpChr = row['CHR']
        snpPos = row['BP']
        tssChr = row['Somamer target TSS Chr']
        tssList = str(row['Somamer target TSS']).split('|')

        # Skip if no TSS
        if tssChr == 'NA':
            snpDf.at[idx, 'cisTrans'] = 'NA'
            snpDf.at[idx, 'Genetic distance'] = 'Non-human'
            continue

        # Sex chromosomes
        if tssChr in ['X', 'Y']:
            snpDf.at[idx, 'cisTrans'] = 'Trans'
            snpDf.at[idx, 'Genetic distance'] = 'Different chr'
            continue

        # Handle multiple TSS chromosomes
        if '|' in tssChr:
            tssChrList = list(map(str, tssChr.split('|')))
            tssPosList = list(map(int, tssList))

            if str(snpChr) in tssChrList:
                # Match TSS positions on the same chromosome
                matchingIndexes = [i for i, x in enumerate(tssChrList) if x == str(snpChr)]
                matchingTssPos = [tssPosList[i] for i in matchingIndexes]

                distances = [abs(int(tssPos) - snpPos) for tssPos in matchingTssPos]
                closestDistance = min(distances)
                snpDf.at[idx, 'Genetic distance'] = closestDistance
                snpDf.at[idx, 'cisTrans'] = 'Cis' if closestDistance < 1_000_000 else 'Trans'
            else:
                snpDf.at[idx, 'cisTrans'] = 'Trans'
                snpDf.at[idx, 'Genetic distance'] = 'Different chr'
        else:
            # Single TSS case
            try:
                if int(tssChr) == snpChr:
                    distance = abs(int(tssList[0]) - snpPos)
                    snpDf.at[idx, 'Genetic distance'] = distance
                    snpDf.at[idx, 'cisTrans'] = 'Cis' if distance < 1_000_000 else 'Trans'
                else:
                    snpDf.at[idx, 'cisTrans'] = 'Trans'
                    snpDf.at[idx, 'Genetic distance'] = 'Different chr'
            except:
                snpDf.at[idx, 'cisTrans'] = 'NA'
                snpDf.at[idx, 'Genetic distance'] = 'Invalid'

    return snpDf

# --- Run pipeline ---
snpDf = addTssAnnotations(snpDf, tssDf)
snpDf = allocateCisTrans(snpDf)

# --- Save annotated DataFrame ---
snpDf.to_csv(outputPath, index=False)
print(f'Saved annotated results to: {outputPath}')
