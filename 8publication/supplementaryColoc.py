import pandas as pd
import os
import ieugwaspy

supplementaryColocColumns = ['HUGO', 'Uniprot', 'somamerID', 'Outcome', 'Open GWAS Dataset ID', 'Open GWAS Sample size', 'Number of SNP', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Paths ---
resultsDir = os.path.join(projectRoot, 'Scripts/colocalization/results/')
MRresultsDir = os.path.join(projectRoot, 'MRResults/')
dfSupplementary1 = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')

dfSupplementary1 = pd.read_excel(dfSupplementary1)

result = {}

for i in supplementaryColocColumns:
    result[i] = []

def getOutcome(MRdf, studyID):

    MRdf['studyID'] = ''
    for i in range(len(MRdf)):
        MRdf.at[i, 'studyID'] = MRdf.loc[i, 'outcome'].split(' || id:')[-1]
        MRdf.at[i, 'outcome'] = MRdf.loc[i, 'outcome'].split(' || id:')[0]

    resultIndex = MRdf.index[MRdf['studyID'] == studyID].tolist()[0]

    result = MRdf.loc[resultIndex, 'outcome']
    return result

def getColocalizationResults(resultsDir, HUGO, datasetID):

    with open(resultsDir + HUGO + '_' + datasetID + '.tsv.gz') as f:
        for num, line in enumerate(f):
            
            if num == 1:
                NSNP = line.split('\t')[1][:-1]
            if num == 2:
                H0 = line.split('\t')[1][:-1]
            elif num == 3:
                H1 = line.split('\t')[1][:-1]
            elif num == 4:
                H2 = line.split('\t')[1][:-1]
            elif num == 5:
                H3 = line.split('\t')[1][:-1]
            elif num == 6:
                H4 = line.split('\t')[1][:-1]   
                
    return [NSNP, H0, H1, H2, H3, H4]

#Uses ieugwaspy to retrieve sample size. Some studies do not have a sample size, then the function calculates it by summing cases and controls
def getSampleSize(datasetID):

    try:
        result = ieugwaspy.gwasinfo([datasetID])[0]['sample_size']

    except KeyError:
        result = ieugwaspy.gwasinfo([datasetID])[0]['ncase'] + ieugwaspy.gwasinfo([datasetID])[0]['ncontrol']
    
    return result
    

for filename in os.listdir(resultsDir):
    if filename.endswith('.tsv.gz'):

        # print(filename)
        supplementary1Index = dfSupplementary1.index[dfSupplementary1['HUGO'] == filename.split('_')[0]].tolist()
        
        result['HUGO'].append(filename.split('_')[0])
        result['Uniprot'].append(dfSupplementary1.loc[supplementary1Index, 'Uniprot'].tolist()[0])
        result['somamerID'].append(dfSupplementary1.loc[supplementary1Index, 'somamerID'].tolist()[0])
        result['Open GWAS Dataset ID'].append(filename.split('_')[-1][:-7])

        MRdf = pd.read_excel(MRresultsDir + result['somamerID'][-1] + '.xlsx')
        result['Outcome'].append(getOutcome(MRdf, result['Open GWAS Dataset ID'][-1]))

        result['Open GWAS Sample size'].append(getSampleSize(result['Open GWAS Dataset ID'][-1]))

        colocalizationResults = getColocalizationResults(resultsDir, result['HUGO'][-1], result['Open GWAS Dataset ID'][-1])

        for num, i in enumerate(['Number of SNP', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']):
            result[i].append(float(colocalizationResults[num]))

def annotateOutcomeClass(df, annotationDf):
    df['Outcome group'] = ''
    for i in range(len(df)):
        if df.loc[i, 'Outcome'] in annotationDf['outcome'].to_list():
            annotationIndex = annotationDf.index[annotationDf['outcome'] == df.loc[i, 'Outcome']][0]
            df.at[i, 'Outcome group'] = annotationDf.loc[annotationIndex, 'Type']

        else:
            print('Warning, outcome', df.loc[i, 'Outcome'], 'not found in the annotation dataframe!')

    return df 

def reorderCols(df, columnName, newPosition):
    columns = df.columns.to_list()
    oldPosition = columns.index(columnName)
    columns.insert(newPosition,columns.pop(oldPosition))
    df = df[columns]
    return df
    
result = pd.DataFrame(result)

# result = annotateOutcomeClass(result, groupedOutcomesAnnotation)
result = result.sort_values(by = ['PP.H4.abf'], ascending = False)

newOrder = ['HUGO', 'Uniprot', 'somamerID', 'Outcome', 'Open GWAS Dataset ID', 'Open GWAS Sample size', 'Number of SNP', 'PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']

for num, i in enumerate(newOrder):
    result = reorderCols(result, i, num)

# --- Define Output Path ---
colocSupplementaryPath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryColoc.xlsx')

# --- Save File ---
result.to_excel(colocSupplementaryPath, index=False)