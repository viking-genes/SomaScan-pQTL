# --- Setup ---
import os
import pandas as pd

# --- Paths ---
projectRoot = './prj_190_viking_somalogic/'

olinkPath = os.path.join(projectRoot, 'Scripts/SNPnovelty/previouslyPublished/uniprotOlink1500Map.xlsx')
somalogicPanelPath = os.path.join(projectRoot, 'Scripts/SNPnovelty/previouslyPublished/SomaLogicv4.0panelFerkingstad2021.xlsx')
combinedPanelPath = os.path.join(projectRoot, 'Scripts/SNPnovelty/somalogic4Olink1500Combined.xlsx')
annotationOutPath = os.path.join(projectRoot, 'Scripts/SNPnovelty/somalogic7kAnnotation.xlsx')
flaggedPath = os.path.join(projectRoot, 'GWAS/viking1_proteomics_somalogic7k_2021_flagged.csv')

# --- Olink processing ---
olinkDf = pd.read_excel(olinkPath, engine='openpyxl')

def extractReviewedProteins(df):
    test = []
    result = df[df['Reviewed'] == 'reviewed']
    for i in df['From']:
        if i not in result['From'].values:
            test.append(i)
    if len(test) > 0:
        print('WARNING: The following proteins are not reviewed: ' + str(test))
    result.reset_index(drop=True, inplace=True)
    result.drop(columns=['Reviewed', 'Gene Names (synonym)'], inplace=True)
    return result

def combineDuplicates(df):
    duplicateRows = df[df.duplicated(['From'])]
    duplicateIndexes = duplicateRows.index.values.tolist()
    result = df.groupby(['From']).agg({'Entry': lambda x: '|'.join(x), 'Gene Names': lambda x: '|'.join(x)}).reset_index()
    result.to_excel(os.path.join(projectRoot, 'Scripts/SNPnovelty/previouslyPublished/uniprotOlink1500MapReviewedCombined.xlsx'), index=False)
    result['Gene Names'] = result['Gene Names'].apply(lambda x: '|'.join(set(x.split('|'))))
    for i in range(len(result)):
        geneNames = result['Gene Names'][i].split('|')
        for j in geneNames:
            if j == result['From'][i]:
                geneNames.remove(j)
                result.at[i, 'Gene Names'] = '|'.join(geneNames)
                continue
    result = result[['From', 'Entry']]
    result.rename(columns={'From': 'Gene', 'Entry': 'UniProt'}, inplace=True)
    return result

def addMultiTarget(df, dataset):
    OlinkMultitargets = ['DEFA1_DEFA1B', 'EBI3_IL27', 'LGALS7_LGALS7B', 'DEFB4A_DEFB4B', 'FUT3_FUT5', 'IL12A_IL12B', 'MICB_MICA', 'CKMT1A_CKMT1B']
    result = df.copy()
    result['Multi-target ' + dataset] = 'No'
    multiTarget = {}
    if dataset == 'Olink':
        for i in OlinkMultitargets:
            targets = i.split('_')
            for j in targets:
                if j not in multiTarget.keys():
                    save = targets.copy()
                    save.remove(j)
                    multiTarget[j] = save
    for i in range(len(result)):
        if result['Gene'][i] in multiTarget.keys(): 
            result.at[i, 'Multi-target ' + dataset] = '|'.join(multiTarget[result['Gene'][i]])
    return result

olinkDf = extractReviewedProteins(olinkDf)
olinkDf = combineDuplicates(olinkDf)
olinkDf = addMultiTarget(olinkDf, 'Olink')

# --- SomaLogic panel processing ---
somalogicDf = pd.read_excel(somalogicPanelPath, skiprows=2)

def extractProteinsSomalogic(df):
    result = df[df['Organism'].str.contains('human')]
    result = result[result['UniProt'].notna()]
    result = result[['Gene', 'UniProt']]
    result.drop_duplicates(subset=['Gene', 'UniProt'], inplace=True)
    df = df.reset_index(drop=True)
    return result

def separateMultiTarget(df):
    multitargetIndexes = []
    df = df.reset_index(drop=True)
    for i in range(len(df)):
        if '|' in df['UniProt'][i]:
            multitargetIndexes.append(i)
    df['Multi-target Somalogic'] = 'No'
    for i in multitargetIndexes:
        gene = df['Gene'][i].split('|')
        uniprot = df['UniProt'][i].split('|')
        for j in range(len(gene)):
            df = df.append({'Gene': gene[j], 'UniProt': uniprot[j], 'Multi-target Somalogic': '|'.join(gene)}, ignore_index=True)
    for i in range(len(df)):
        if df['Gene'][i] in df['Multi-target Somalogic'][i]:
            df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'].replace(df['Gene'][i], '')
            if df['Multi-target Somalogic'][i][0] == '|':
                df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][1:]
            elif df['Multi-target Somalogic'][i][-1] == '|':
                df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][:-1]
    df.drop(multitargetIndexes, inplace=True)
    df.reset_index(drop=True, inplace=True)
    duplicateRows = df[df.duplicated(['Gene'], keep=False)].copy()
    duplicateRows.reset_index(drop=True, inplace=True)
    df.drop_duplicates(subset=['Gene'], keep=False, inplace=True)
    for i in range(len(duplicateRows)):
        if duplicateRows['Gene'][i] not in df['Gene'].tolist():
            df = df.append(duplicateRows.loc[i, :], ignore_index=True)
        else:
            index = df[df['Gene'] == duplicateRows['Gene'][i]].index[0]
            if duplicateRows['UniProt'][i] not in df['UniProt'][index]:
                df.at[index, 'UniProt'] = df.loc[index, 'UniProt'] + '|' + duplicateRows.loc[i, 'UniProt']
            df.at[index, 'Multi-target Somalogic'] = df.loc[index, 'Multi-target Somalogic'] + '|' + duplicateRows.loc[i, 'Multi-target Somalogic']
    for i in range(len(df)):
        multitargets = df['Multi-target Somalogic'][i].split('|')
        multitargets = list(set(multitargets))
        df.at[i, 'Multi-target Somalogic'] = '|'.join(multitargets)
        if df['Multi-target Somalogic'][i][0] == '|':
            df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][1:]
        elif df['Multi-target Somalogic'][i][-1] == '|':
            df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][:-1]
    for i in range(len(df)):
        if '|' in df['Multi-target Somalogic'][i]:
            multitargets = df['Multi-target Somalogic'][i].split('|')
            multitargets = list(set(multitargets))
            df.at[i, 'Multi-target Somalogic'] = '|'.join(multitargets)
            if df['Multi-target Somalogic'][i][0] == '|':
                df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][1:]
            elif df['Multi-target Somalogic'][i][-1] == '|':
                df.at[i, 'Multi-target Somalogic'] = df.loc[i, 'Multi-target Somalogic'][:-1]
            multitargets = df['Multi-target Somalogic'][i].split('|')
            if 'No' in multitargets:
                multitargets.pop(multitargets.index('No'))
                multitargets.sort()
                multitargets.insert(0, 'No')
                df.at[i, 'Multi-target Somalogic'] = '|'.join(multitargets)
    return df

somalogicDf = extractProteinsSomalogic(somalogicDf)
somalogicDf = separateMultiTarget(somalogicDf)

def combineDataframes(dfSomalogic, dfOlink):
    def extractDuplicates(df1, df2):
        dfDuplicates = df1[df1['Gene'].isin(df2['Gene'])].index.tolist()
        df1.loc[dfDuplicates, 'Panel'] = 'SomaLogic|Olink'
        return df1

    def removeDuplicates(df1):
        duplicateIndexes = df1[df1['Panel'] == 'SomaLogic|Olink'].index.tolist()
        for i in duplicateIndexes:
            gene = df1['Gene'][i]
            uniprot = df1['UniProt'][i]
            dfIndex = df1[df1['Gene'] == gene].index.tolist()
            for j in dfIndex:
                if j != i:
                    if uniprot in df1['UniProt'][j].split('|'):
                        df1.at[i, 'Multi-target Olink'] = df1.loc[j, 'Multi-target Olink']
                        df1.drop(j, inplace=True)
        return df1

    dfSomalogic['Panel'] = 'SomaLogic'
    dfOlink['Panel'] = 'Olink'
    duplicates = extractDuplicates(dfSomalogic, dfOlink)
    result = pd.concat([dfSomalogic, dfOlink])
    result.reset_index(drop=True, inplace=True)
    result = removeDuplicates(result)
    return result

combinedDf = combineDataframes(somalogicDf, olinkDf)
combinedDf.rename(columns={'Gene': 'Entrez Gene Symbol', 'UniProt': 'UniProt', 'Multi-target Somalogic': 'Multi-target SomaLogic', 'Multi-target Olink': 'Multi-target Olink'}, inplace=True)
combinedDf.to_excel(combinedPanelPath, index=False)

# --- SomaLogic 7k Annotation from flagged file ---
df = pd.read_csv(flaggedPath, low_memory=False)

columnsToExtract = ['HUGO', 'Uniprot', 'Full Name', 'Entrez Gene Symbol', 'Sequence ID', 'Somamer ID']
result = [[] for _ in range(len(columnsToExtract))]

for i in df.columns[2:]:
    result[0].append(df.loc[4, i])
    result[1].append(df.loc[5, i])
    result[2].append(df.loc[3, i])
    result[3].append(df.loc[7, i])
    result[4].append(df.loc[0, i])
    result[5].append(df.loc[2, i])

annotationDf = pd.DataFrame(result).transpose()
annotationDf.columns = columnsToExtract
annotationDf = annotationDf.rename(columns={'Full Name': 'Full Protein Name'})
annotationDf.to_excel(annotationOutPath, index=None)
