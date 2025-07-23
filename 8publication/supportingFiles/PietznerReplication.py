import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Directory Paths ---
df = pd.read_excel(os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx'))
outputDf = os.path.join(projectRoot, 'Scripts/Publication/supplementary1Replication.xlsx')
pietznerDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/Pietzner2021SummaryStats')


df['Replication p-value'] = ''
df['Replication Effect Size'] = ''
df['Replication Allele Frequency'] = ''
df['Replication SomamerID'] = ''

# The authors (Pietzner et al.) report values as 0 having used bgenie. Set these values as 1x10-300 
def assignMinPValue(df):
    df['Replication p-value'] = df['Replication p-value'].replace(str(0.0), str(1e-300))

    for index in df.index:
        if '|' in df.loc[index, 'Replication p-value']:
            pValueList = df.loc[index, 'Replication p-value'].split('|')
            pValueList = [str(1e-300) if x == str(0.0) else x for x in pValueList]
            df.at[index, 'Replication p-value'] = '|'.join(pValueList)
            
    return df


allReplicationFilesList = os.listdir(pietznerDir)

#Find last processed index in the Replication Effect Size column
lastProcessedIndex = 0
if os.path.exists(outputDf):
    df = pd.read_excel(outputDf, na_values=['NA'], keep_default_na=False)
    #Find the first index not equal to ''
    if '' in df['Replication Effect Size'].values:
        lastProcessedIndex = df[df['Replication Effect Size'] == ''].index[0]
        print('Proceeding from last stop point at index', lastProcessedIndex)
    else:
        print('Database already processed')
        lastProcessedIndex = len(df)

for index in range(lastProcessedIndex, len(df)):

    if index % 10 == 0 and index != 0:
        print('Processing row', index, 'of', len(df))
        #Save progress replacing nans in the three columns processed with this script with NA
        df['Replication p-value'] = df['Replication p-value'].fillna('NA')
        df['Replication Effect Size'] = df['Replication Effect Size'].fillna('NA')
        df['Replication Allele Frequency'] = df['Replication Allele Frequency'].fillna('NA')
        df.to_excel(outputDf, index=False)

    uniprot = df.loc[index, 'Uniprot']
    rsid = df.loc[index, 'SNP']

    #Some rsid in VIKING are named differently and were annotated manually in dbsnp by chr/pos/alleles to allow crossreferencing in other studies
    rsidCorrectionDict = {'9:140008750_G_C': 'rs10747049', 'rs1057114': 'rs72620874'}
    if rsid in rsidCorrectionDict:
        rsid = rsidCorrectionDict[rsid]

    proteinsMissingInPietzner = ['O00391', 'P21549', 'P78410', 'O15354', 'P16455', 'Q9BRF8', 'P50225', 'Q16853', 'Q6UX27', 'P04080']
    if '|' in uniprot or df.loc[index, 'cisTrans'] == 'Trans' or df.loc[index, 'novel'] == True or uniprot in proteinsMissingInPietzner:
        df.at[index, 'Replication p-value'] = 'NA'
        df.at[index, 'Replication Effect Size'] = 'NA'
        df.at[index, 'Replication Allele Frequency'] = 'NA'
        continue
    
    #Some SNP were not genotyped and do not have any proxies R2>0.6
    SNPMissingInPietzner = ['rs4149358']
    if rsid in SNPMissingInPietzner:
        df.at[index, 'Replication p-value'] = 'Not genotyped'
        df.at[index, 'Replication Effect Size'] = 'Not genotyped'
        df.at[index, 'Replication Allele Frequency'] = 'Not genotyped'
        continue

    #Load Pietzner data for all proteins matching the uniprot
    replicationUniprot = [file for file in allReplicationFilesList if uniprot in file]

    if len(replicationUniprot) == 0:
        print('No replication data found for', uniprot)
        exit()
    
    replicationPvalueList = []
    replicationEffectSizeList = []
    replicationAlleleFrequencyList = []
    replicationSomamerIDList = []
    
    for fileName in replicationUniprot:
            
        replicationDf = pd.read_csv(os.path.join(pietznerDir, fileName), sep='\t')

        #Find index of the pQTL rsid found in VIKING
        replicationRsidIndex = replicationDf[replicationDf['rsid'] == rsid].index.tolist()
        
        if len(replicationRsidIndex) == 0:
            print('Warning, rsid', rsid, 'not found in', fileName)
            exit()

        replicationRsidIndex = replicationRsidIndex[0]

        #Annotate p-value, effect size and allele frequency
        #Format p-value to 2 decimal places
        pValue = "{:.2e}".format(replicationDf.loc[replicationRsidIndex, 'Pvalue'])
        pValue = float(pValue)
        
        effectSize = replicationDf.loc[replicationRsidIndex, 'Effect']
        #Round allele frequency to 4 decimal places
        frequency = round(replicationDf.loc[replicationRsidIndex, 'Freq1'], 4)

        #Compare and flip the effect allele in VIKING to the effect allele in Pietzner in Allele1 and Allele2 columns
        if df.loc[index, 'Effect allele'] != replicationDf.loc[replicationRsidIndex, 'Allele1'].upper():
            #Check if allele is flipped
            if df.loc[index, 'Effect allele'] == replicationDf.loc[replicationRsidIndex, 'Allele2'].upper() and df.loc[index, 'Other allele'] == replicationDf.loc[replicationRsidIndex, 'Allele1'].upper():
                effectSize = -1 * effectSize
                frequency = round(1 - frequency, 4)
            else:
                print('Warning, effect allele does not match for', rsid, 'in', fileName)
                print(replicationDf.loc[replicationRsidIndex, :])
                exit()

        somamerID = fileName.split('.')[0]
        somamerID = '-'.join(somamerID.split('_')[-2:])

        replicationPvalueList.append(pValue)
        replicationEffectSizeList.append(effectSize)
        replicationAlleleFrequencyList.append(frequency)
        replicationSomamerIDList.append(somamerID)

    #Annotate supplementary1 with the replication data
    df.loc[index, 'Replication p-value'] = '|'.join([str(x) for x in replicationPvalueList])
    df.loc[index, 'Replication Effect Size'] = '|'.join([str(x) for x in replicationEffectSizeList])
    df.loc[index, 'Replication Allele Frequency'] = '|'.join([str(x) for x in replicationAlleleFrequencyList])
    df.loc[index, 'Replication SomamerID'] = '|'.join([x for x in replicationSomamerIDList])

#Save replacing nans in the three columns processed with this script with NA
df['Replication p-value'] = df['Replication p-value'].fillna('NA')
df['Replication Effect Size'] = df['Replication Effect Size'].fillna('NA')
df['Replication Allele Frequency'] = df['Replication Allele Frequency'].fillna('NA')

#Convert 0 p-values to 1x10-300
df = assignMinPValue(df)

df.to_excel(outputDf, index=False)
