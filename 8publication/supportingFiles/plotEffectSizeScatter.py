import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Output Paths ---
dfPietzner = os.path.join(projectRoot, 'Scripts/Publication/supplementary1Replication.xlsx')
dfSurapaneni = os.path.join(projectRoot, 'Scripts/Publication/supplementary1ReplicationSurapaneni.xlsx')
outputDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication')


dfPietzner = pd.read_excel(dfPietzner)
dfSurapaneni = pd.read_excel(dfSurapaneni)

def plotScatter(df, x, y, title):

    fig, ax = plt.subplots()
    ax.scatter(x, y, alpha=0.5)
    ax.set_xlabel('VIKING Effect Size')
    ax.set_ylabel('Replication Effect Size')
    ax.set_title(title)

    #Add Pearson R2 to plot
    plt.text(min(x), max(y), 'R2 = ' + str(round(np.corrcoef(x, y)[0][1]**2, 2)))

    #Make transparent
    fig.patch.set_alpha(0.0)

    #Add y = x line
    ax.plot([min(x), max(x)], [min(x), max(x)], color='black', linestyle='dashed')

    plt.savefig(os.path.join(outputDir, title + '.png'))


#Count number of novel cis SNP in VIKING and how many of them replicated in Surapaneni
novelSNP = 0
replicatedSNP = 0
for index in dfSurapaneni.index:
    if dfSurapaneni.loc[index, 'novel'] and dfSurapaneni.loc[index, 'cisTrans'] == 'Cis':
        novelSNP += 1
        if dfSurapaneni.loc[index, 'Replication AASK Effect Size'] == dfSurapaneni.loc[index, 'Replication AASK Effect Size']:
            replicatedSNP += 1

print('Number of novel SNP in VIKING:', novelSNP)
print('Number of novel SNP replicated in Surapaneni2022:', replicatedSNP)

#Drop rows with missing values (no replication)
for index in dfSurapaneni.index:
    if dfSurapaneni.loc[index, 'Replication AASK Effect Size'] != dfSurapaneni.loc[index, 'Replication AASK Effect Size']:
        dfSurapaneni.drop(index, inplace=True)
        
dfSurapaneni.reset_index(drop=True, inplace=True)
print('Number of replicated VIKING pQTL in Surapaneni2022:', len(dfSurapaneni))

#If replication effect size is a list, copy the row and split the effect size and other data into separate rows
for index in dfSurapaneni.index:
    if '|' in dfSurapaneni.loc[index, 'Replication AASK Effect Size']:
        effectSizeList = dfSurapaneni.loc[index, 'Replication AASK Effect Size'].split('|')
        effectSizeList = [float(x) for x in effectSizeList]

        rsidList = dfSurapaneni.loc[index, 'Replication AASK rsid'].split('|')
        LDList = dfSurapaneni.loc[index, 'Replication AASK LD'].split('|')
        pValueList = dfSurapaneni.loc[index, 'Replication AASK p-value'].split('|')
        alleleFrequencyList = dfSurapaneni.loc[index, 'Replication AASK Allele Frequency'].split('|')
        somamerIDList = dfSurapaneni.loc[index, 'Replication AASK SomamerID'].split('|')

        new_rows = []
        for i in range(len(effectSizeList))[1:]:
            new_row = dfSurapaneni.loc[index].copy()
            new_row['Replication AASK Effect Size'] = effectSizeList[i]
            new_row['Replication AASK rsid'] = rsidList[i]
            new_row['Replication AASK LD'] = LDList[i]
            new_row['Replication AASK p-value'] = pValueList[i]
            new_row['Replication AASK Allele Frequency'] = alleleFrequencyList[i]
            new_row['Replication AASK SomamerID'] = somamerIDList[i]
            new_rows.append(new_row)

        dfSurapaneni = pd.concat([dfSurapaneni, pd.DataFrame(new_rows)], ignore_index=True)

        dfSurapaneni.loc[index, 'Replication AASK Effect Size'] = effectSizeList[0]
        dfSurapaneni.loc[index, 'Replication AASK rsid'] = rsidList[0]
        dfSurapaneni.loc[index, 'Replication AASK LD'] = LDList[0]
        dfSurapaneni.loc[index, 'Replication AASK p-value'] = pValueList[0]
        dfSurapaneni.loc[index, 'Replication AASK Allele Frequency'] = alleleFrequencyList[0]
        dfSurapaneni.loc[index, 'Replication AASK SomamerID'] = somamerIDList[0]


dfSurapaneni.reset_index(drop=True, inplace=True)
print('Number of total associations in Surapaneni2022:', len(dfSurapaneni))

plotScatter(dfSurapaneni, dfSurapaneni['Effect Size (beta)'], dfSurapaneni['Replication AASK Effect Size'].astype(float), 'Surapaneni2022 Replication by Effect Size')


###Pietzner replication 
#Drop rows with missing values (not measured in Pietzner study)
for index in dfPietzner.index:
    if dfPietzner.loc[index, 'Replication Effect Size'] != dfPietzner.loc[index, 'Replication Effect Size'] or dfPietzner.loc[index, 'Replication Effect Size'] == 'Not genotyped':
        dfPietzner.drop(index, inplace=True)

dfPietzner.reset_index(drop=True, inplace=True)
print('Total number of VIKING pQTL matched by protein in Pietzner2021:', len(dfPietzner))

#Filter out Pietzner results that do not share the somamerID with VIKING
columns = ['Replication p-value', 'Replication Effect Size', 'Replication Allele Frequency', 'Replication SomamerID']
for index in dfPietzner.index:

    matchedNum = ''
    for num, somamerIDPietzner in enumerate(dfPietzner.loc[index, 'Replication SomamerID'].split('|')):
        if somamerIDPietzner == dfPietzner.loc[index, 'somamerID']:
            matchedNum = num
            break

    if matchedNum == '':
        print(dfPietzner.loc[index, 'somamerID'], dfPietzner.loc[index, 'Uniprot'])
        dfPietzner.drop(index, inplace=True)
        continue

    #Reannotate using only the matched somamerID
    for column in columns:
        dfPietzner.loc[index, column] = dfPietzner.loc[index, column].split('|')[matchedNum]

dfPietzner.reset_index(drop=True, inplace=True)
print('Number of VIKING pQTL with matching somamerID in Pietzner2021:', len(dfPietzner))

#Create a dictionary of uniprot and their associated somamerID
uniprotDict = {}

for index in dfPietzner.index:
    if dfPietzner.loc[index, 'Uniprot'] not in uniprotDict:
        uniprotDict[dfPietzner.loc[index, 'Uniprot']] = {'SomamerIDList': [dfPietzner.loc[index, 'somamerID']]}
    else:
        uniprotDict[dfPietzner.loc[index, 'Uniprot']]['SomamerIDList'].append(dfPietzner.loc[index, 'somamerID'])

print('Number of unique proteins in the replication dataset:', len(uniprotDict))

print('Proteins with multiple aptamers targeting them:')
#Print the proteins with multiple somamerID
for key in uniprotDict:
    if len(uniprotDict[key]['SomamerIDList']) > 1:
        print(key, uniprotDict[key]['SomamerIDList'])

#Annotate each protein with whether it's effect size was replicated in the same direction
def calculateReplicationStandardDeviation(df, uniprotDict):

    for key in uniprotDict:
        uniprotDict[key]['EffectSizeDirectionConsistent'] = {}

    for index in df.index:
        somamerID = df.loc[index, 'somamerID']
        uniprotID = df.loc[index, 'Uniprot']

        effectSizeVIKING = float(df.loc[index, 'Effect Size (beta)'])
        effectSizeSEVIKING = float(df.loc[index, 'Standard error of beta'])
        effectSizeReplication = float(df.loc[index, 'Replication Effect Size'])

        #Check if directionally consistent effect size
        if effectSizeVIKING * effectSizeReplication > 0:
            uniprotDict[uniprotID]['EffectSizeDirectionConsistent'][somamerID] = True
        
    return uniprotDict

uniprotDict = calculateReplicationStandardDeviation(dfPietzner, uniprotDict)

print('Number of proteins with replication directionally consistent:', len([key for key in uniprotDict if True in uniprotDict[key]['EffectSizeDirectionConsistent'].values()]))

#Plot scatter plot of effect size
plotScatter(dfPietzner, dfPietzner['Effect Size (beta)'], dfPietzner['Replication Effect Size'].astype(float), 'cis pQTL replication in Pietzner2021 et al.')