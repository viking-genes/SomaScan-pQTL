import pandas as pd
import os
import numpy as np

# Suppress the SettingWithCopyWarning
pd.options.mode.chained_assignment = None

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

# Define input path
replicationPreprocessingPath = os.path.join(
    projectRoot, 'Scripts/Publication/supportingFiles/replication/supplementaryReplicationPreprocessing.xlsx'
)

# Load data
df = pd.read_excel(replicationPreprocessingPath)

naValues = ['NULL', 'NaN', 'nan', '-nan', 'n/a']
df = pd.read_excel(df, na_values=naValues, keep_default_na=False)

def matchAlleleFlip(df, i, totalCount, dropCount, assay):

    #Returns the in phase alleles in the form of {'EffectAlleleVIKING': 'X', 'InphaseAllele': 'Y', 'OtherAlleleVIKING': 'A', 'OutOfPhaseAllele': 'B'}
    def extractInPhaseAlleles(filepath):
        result = {'EffectAlleleVIKING': '', 'InphaseAllele': '', 'OtherAlleleVIKING': '', 'OutOfPhaseAllele': ''}
        # Open the file for reading
        with open(filepath, 'r') as file:

            haplotypeLineStart = len(file.readlines())

            #Extract the haplotypes from the LD file and return the in phase alleles
            file.seek(0)
            for num, line in enumerate(file):

                # Check if the line contains 'In phase alleles'
                if 'In phase alleles' in line:
                    #Return the line that contains the 'In phase alleles'
                    resultAlleles = line.strip().split(' ')[-1]
                    resultAlleles = resultAlleles.split('/')
                    
                    #Split the alleles if they are single nucleotide
                    if len(resultAlleles[0]) == 2 and len(resultAlleles[1]) == 2:
                        result['EffectAlleleVIKING'] = resultAlleles[0][0]
                        result['InphaseAllele'] = resultAlleles[0][1]
                        result['OtherAlleleVIKING'] = resultAlleles[1][0]
                        result['OutOfPhaseAllele'] = resultAlleles[1][1]
                        return result
                        
                    else:
                        print('Error: In phase alleles are not single nucleotide', resultAlleles, filepath)
                        return None

        print('Error: No in phase alleles found in file', filepath)
        exit()

    def flipAlleles(effectAllele, otherAllele, allelePhase, effectSize):
        # {'EffectAlleleVIKING': 'X', 'InphaseAllele': 'Y', 'OtherAlleleVIKING': 'A', 'OutOfPhaseAllele': 'B'}

        #If the alleles out of phase are the same as the reference alleles, flip the alleles and effect size
        if effectAllele == allelePhase['OutOfPhaseAllele'] and otherAllele == allelePhase['InphaseAllele']:
            return otherAllele, effectAllele, -1 * float(effectSize)

        #Assume that NA other allele is wildcard for any allelekeyerror: 'Somalogic Replicatio
        elif effectAllele == allelePhase['InphaseAllele'] and otherAllele == 'NA':
            return effectAllele, 'NA', float(effectSize)

        elif effectAllele == allelePhase['OutOfPhaseAllele'] and otherAllele == 'NA':
            return 'NA', otherAllele, -1 * float(effectSize)
        
        else:
            print('Error: Allele not found in allele phase. Likely strand issues, these are not accounted for, dropping.', effectAllele, otherAllele, allelePhase, refRsid, replicationRsid)
            return None, None, None

    #Not all effect alleles are reported as minor alleles, this function matches the phase to the one reported in the VIKING study
    def matchEffectToMinor(allelePhase, effectVIKING, otherVIKING):
        if allelePhase['EffectAlleleVIKING'] == effectVIKING:
            return allelePhase
        elif allelePhase['EffectAlleleVIKING'] == otherVIKING:
            allelePhase['EffectAlleleVIKING'], allelePhase['OtherAlleleVIKING'] = allelePhase['OtherAlleleVIKING'], allelePhase['EffectAlleleVIKING']
            allelePhase['InphaseAllele'], allelePhase['OutOfPhaseAllele'] = allelePhase['OutOfPhaseAllele'], allelePhase['InphaseAllele']
            return allelePhase


    # Define output directory for allele phasing results
    flipDir = os.path.join(
        projectRoot,
        'Scripts/Publication/supportingFiles/replication/allelePhasingOutput'
    )

    #Skip rows that do not have replication data
    if df.loc[i, assay + ' SNP'] == '':
        return df, totalCount, dropCount
    
    # Extract relevant data from dataframe for row 'i'
    refRsid, refEffectAllele, refAltAllele = df.loc[i, ['SNP', 'Effect allele', 'Other allele']]
    
    # Splitting string data in columns into lists
    replication_columns = [
        assay + ' Effect allele', 
        assay + ' Other allele', 
        assay + ' Effect Size', 
        assay + ' Study', 
        assay + ' SNP', 
        assay + ' SNP LD R2', 
        assay + ' Variance Explained'
    ]

    if assay == 'Somalogic':
        replication_columns.insert(3, assay + ' SomamerID')

        (
            replicationEffectAlleleList, 
            replicationOtherAlleleList, 
            replicationEffectSizeList, 
            replicationSomamerList, 
            replicationStudyList, 
            replicationSNPList, 
            replicationLDR2List, 
            replicationVarEplList
        ) = (df.loc[i, col].split('|') for col in replication_columns)

    if assay == 'Olink':
        (
            replicationEffectAlleleList, 
            replicationOtherAlleleList, 
            replicationEffectSizeList, 
            replicationStudyList, 
            replicationSNPList, 
            replicationLDR2List, 
            replicationVarEplList
        ) = (df.loc[i, col].split('|') for col in replication_columns)

    toDrop = []

    for num, replicationRsid in enumerate(df.loc[i, assay + ' SNP'].split('|')):

        #Counts the number of SNPs that have been processed
        totalCount += 1

        # Access elements from lists by index 'num'
        replicationEffectAllele, replicationOtherAllele, replicationEffectSize = (
            replicationEffectAlleleList[num],
            replicationOtherAlleleList[num],
            replicationEffectSizeList[num]
        )

        #Check if the alleles are flipped when rsid is the same
        if replicationRsid == refRsid:

            if refEffectAllele != replicationEffectAllele:
                #Flip the alleles
                if refEffectAllele == replicationOtherAllele and refAltAllele == replicationEffectAllele:
                    replicationEffectAllele = refEffectAllele
                    replicationOtherAllele = refAltAllele
                    replicationEffectSize = -1 * float(replicationEffectSize)

                #Try matching the alleles with NA as wildcard
                elif refAltAllele == replicationEffectAllele and replicationOtherAllele == 'NA':
                    replicationOtherAllele = replicationEffectAllele
                    replicationEffectAllele = 'NA'
                    replicationEffectSize = float(replicationEffectSize) * -1

                else:
                    print('Error: Same allele and alleles are flipped but do not match. Likely strand issues, these are not accounted for, dropping.', refRsid, refEffectAllele, refAltAllele, replicationEffectAllele, replicationOtherAllele)
                    toDrop.append(num)
        
        #Otherwise, find the in phase alleles and match them and the associated effect sizes
        else:
            #Extract the in phase alleles from a PLINK log file for UKBB 10k reference panel
            allelePhase = extractInPhaseAlleles(os.path.join(flipDir, refRsid + '_' + replicationRsid + '.log'))

            #If the alleles are non-SNV, remove them. PLINK's output does not support non-SNV alleles
            if allelePhase == None:
                toDrop.append(num)
                continue

            #Match the effect allele to the minor allele
            allelePhase = matchEffectToMinor(allelePhase, refEffectAllele, refAltAllele)

            #If the alleles are already in phase, skip
            if refEffectAllele == allelePhase['EffectAlleleVIKING'] and replicationEffectAllele == allelePhase['InphaseAllele']:
                continue

            replicationEffectAllele, replicationOtherAllele, replicationEffectSize = flipAlleles(replicationEffectAllele, replicationOtherAllele, allelePhase, replicationEffectSize)

            #If the alleles are not found in the allele phase, drop the row
            if replicationEffectAllele == None:
                toDrop.append(num)
                continue

        replicationEffectAlleleList[num] = replicationEffectAllele
        replicationOtherAlleleList[num] = replicationOtherAllele
        replicationEffectSizeList[num] = replicationEffectSize

    #Drop rows with errors
    dropCount += len(toDrop)

    if toDrop != []:
        toDrop.reverse()
        for k in toDrop:
            replicationEffectAlleleList.pop(k)
            replicationOtherAlleleList.pop(k)
            replicationEffectSizeList.pop(k)
            replicationStudyList.pop(k)
            replicationSNPList.pop(k)
            replicationLDR2List.pop(k)
            replicationVarEplList.pop(k)

            if assay == 'Somalogic':
                replicationSomamerList.pop(k)

    #Convert to str
    replicationEffectSizeList = [str(x) for x in replicationEffectSizeList]

    df.at[i, assay + ' Effect allele'] = '|'.join(replicationEffectAlleleList)
    df.at[i, assay + ' Other allele'] = '|'.join(replicationOtherAlleleList)
    df.at[i, assay + ' Effect Size'] = '|'.join(replicationEffectSizeList)
    df.at[i, assay + ' Study'] = '|'.join(replicationStudyList)
    df.at[i, assay + ' SNP'] = '|'.join(replicationSNPList)
    df.at[i, assay + ' SNP LD R2'] = '|'.join(replicationLDR2List)
    df.at[i, assay + ' Variance Explained'] = '|'.join(replicationVarEplList)

    if assay == 'Somalogic':
        df.at[i, assay + ' SomamerID'] = '|'.join(replicationSomamerList)
    
    return df, totalCount, dropCount

#Merge rows with the same protein
def mergeSomamer(df):
    
    df['VIKING SomamerID'] = ''
    df['VIKING Effect Size'] = ''
    
    result = []

    uniqueUniprotList = df['Uniprot'].unique()

    #Exclude rows with multiple target proteins
    uniqueUniprotList = [x for x in uniqueUniprotList if '|' not in x]

    for uniprotID in uniqueUniprotList:

        #Retrieve all rows with the current UniprotID
        dfRows = df[df['Uniprot'] == uniprotID]

        #If only one SNP per Protein, add the data
        if len(dfRows) == 1:
            dfRows['VIKING SomamerID'] = df.loc[dfRows.index[0], 'somamerID']
            dfRows['VIKING Effect Size'] = str(round(df.loc[dfRows.index[0], 'Effect Size (beta)'], 4))
            result.append(dfRows)
            continue

        #Drop rows without any replication if there are any rows with replication. If there are no rows with replication data, append the row with the lowest p-value instead

        #Check if there is any replication data
        replicationPresent = False
        for index in dfRows.index:
            if df.loc[index, 'Somalogic Study'] != '' or df.loc[index, 'Olink Study'] != '':
                replicationPresent = True
                break
        
        #If no replication data is present, append the row with the lowest p-value
        if replicationPresent == False:
            minPvalue = float('inf')
            minIndex = 0
            for index in dfRows.index:
                if df.loc[index, 'P-value'] < minPvalue:
                    minPvalue = df.loc[index, 'P-value']
                    minIndex = index

            dfRows.loc[minIndex, 'VIKING SomamerID'] = df.loc[minIndex, 'somamerID']
            dfRows.loc[minIndex, 'VIKING Effect Size'] = str(round(df.loc[minIndex, 'Effect Size (beta)'], 4))
            
            dfRows = dfRows.drop(dfRows.index[dfRows.index != minIndex])
            result.append(dfRows)
            continue

        #Drop rows without replication data
        if replicationPresent == True:
            for index in dfRows.index:
                if df.loc[index, 'Somalogic Study'] == '' and df.loc[index, 'Olink Study'] == '':
                    dfRows.drop(index, inplace=True)

        #If only one row remains, add the data
        if len(dfRows) == 1:
            result.append(dfRows)
            dfRows.loc[dfRows.index[0], 'VIKING SomamerID'] = df.loc[dfRows.index[0], 'somamerID']
            dfRows.loc[dfRows.index[0], 'VIKING Effect Size'] = str(round(df.loc[dfRows.index[0], 'Effect Size (beta)'], 4))
            continue

        # If there are multiple sentinel SNP reported, choose the one with Somalogic replication. Then, if there are multiple, the highest average LD in Somalogic studies
        if len(dfRows['SNP'].unique()) > 1:
            maxLD = 0
            maxIndex = 0
            for index in dfRows.index:
                LDvalues = dfRows.loc[index, 'Somalogic SNP LD R2'].split('|')

                #Skip rows without Somalogic LD data
                if LDvalues == ['']:
                    continue

                LDvalues = [float(x) for x in LDvalues]
                LDvalues = np.mean(LDvalues)
                if LDvalues > maxLD:
                    maxLD = LDvalues
                    maxIndex = index
            
            #Leave only rows with sentinel SNP that has the highest average LD
            for index in dfRows.index:
                if dfRows.loc[index, 'SNP'] != dfRows.loc[maxIndex, 'SNP']:
                    dfRows.drop(index, inplace=True)
            
        if len(dfRows) == 1:
            dfRows.loc[dfRows.index[0], 'VIKING SomamerID'] = df.loc[dfRows.index[0], 'somamerID']
            dfRows.loc[dfRows.index[0], 'VIKING Effect Size'] = str(round(df.loc[dfRows.index[0], 'Effect Size (beta)'], 4))
            result.append(dfRows)
            continue

        #If our sentinel SNP is the same in all VIKING protein measurements with different aptamers, merge rows and average the effect size
        if len(dfRows['SNP'].unique()) == 1:

            dfRows.at[dfRows.index[0], 'VIKING SomamerID'] = df.loc[dfRows.index[0], 'somamerID']
            dfRows.at[dfRows.index[0], 'VIKING Effect Size'] = str(df.loc[dfRows.index[0], 'Effect Size (beta)'])

            for index in dfRows.index[1:]:
                dfRows.at[dfRows.index[0], 'VIKING SomamerID'] += '|' + df.loc[index, 'somamerID']
                dfRows.at[dfRows.index[0], 'VIKING Effect Size'] += '|' + str(df.loc[index, 'Effect Size (beta)'])
            
            dfRows = dfRows.drop(dfRows.index[1:])
            
            #Average the effect size reported in VIKING
            effectSizes = dfRows.at[dfRows.index[0], 'VIKING Effect Size'].split('|')
            effectSizes = [float(x) for x in effectSizes]
            dfRows.at[dfRows.index[0], 'VIKING Effect Size'] = str(round(np.mean(effectSizes), 4))

            result.append(dfRows)
            continue
    
    result = pd.concat(result)
    result = result.reset_index(drop=True)
    
    #Add columns for replication data summary
    result['Somalogic Effect Size Average'] = ''
    result['Olink Effect Size Average'] = ''

    for index in result.index:

        if result.loc[index, 'Somalogic Study'] != '':
            effectSizes = result.loc[index, 'Somalogic Effect Size'].split('|')
            effectSizes = [float(x) for x in effectSizes]
            result.at[index, 'Somalogic Effect Size Average'] = str(round(np.mean(effectSizes), 4))

        if result.loc[index, 'Olink Study'] != '':
            effectSizes = result.loc[index, 'Olink Effect Size'].split('|')
            effectSizes = [float(x) for x in effectSizes]
            result.at[index, 'Olink Effect Size Average'] = str(round(np.mean(effectSizes), 4))

    return result


totalCount = 0
dropCount = 0

for i in range(len(df)):

    #Skip rows without replication
    if df.loc[i, 'Somalogic SomamerID'] != df.loc[i, 'Somalogic SomamerID'] and df.loc[i, 'Olink Study'] != df.loc[i, 'Olink Study']:
        continue
    
    #Correct for flipped alleles in replication studies
    for assay in ['Olink', 'Somalogic']:
        df, totalCount, dropCount = matchAlleleFlip(df, i, totalCount, dropCount, assay)

print('Total SNPs processed:', totalCount)
print('Total SNPs dropped:', dropCount)

#Drop columns that are not needed
columnsToDrop = ['Variance Explained', 'Replication Assay', 'Somalogic Variance Explained', 'Olink Variance Explained']
df.drop(columns=columnsToDrop, inplace=True)

#Merge on protein
df = mergeSomamer(df)


supplementaryReplicationPath = os.path.join(
    projectRoot,
    'Scripts/Publication/supplementaryReplication.xlsx'
)

# Export the dataframe
df.to_excel(supplementaryReplicationPath, index=False)
