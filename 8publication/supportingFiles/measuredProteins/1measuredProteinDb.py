import os
import re
import pandas as pd

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
measuredProteinDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/measuredProteins')

outputIntermediateTsv = os.path.join(measuredProteinDir, 'measuredProteinsIntermediateFile.tsv')
outputUniprotListTxt = os.path.join(measuredProteinDir, 'uniquePublishedUniprotIDs.txt')

# --- Parse measured proteins from each assay/study ---
def extractMeasuredProteins():
    proteinDict = {}

    for assay in ['Olink', 'Somalogic']:
        assayDir = os.path.join(measuredProteinDir, assay)
        for fileName in os.listdir(assayDir):
            studyName = fileName.split('ProteinInfo.xlsx')[0]
            df = pd.read_excel(os.path.join(assayDir, fileName))

            for index in df.index:
                uniprot = df.loc[index, 'UniprotID']

                if pd.isna(uniprot):
                    continue
                if assay == 'Somalogic' and str(uniprot).startswith('HCE'):
                    continue
                if uniprot in ['Family', 'None', None]:
                    continue

                # Split multi-target Uniprot field using known separators
                uniprotSplit = re.split(',| |\\||;|\\.', str(uniprot))
                uniprotClean = [x for x in uniprotSplit if x]
                uniprotClean.sort()
                uniprotKey = '~'.join(uniprotClean)

                if uniprotKey not in proteinDict:
                    proteinDict[uniprotKey] = {'Assay': [], 'Study': []}
                if assay not in proteinDict[uniprotKey]['Assay']:
                    proteinDict[uniprotKey]['Assay'].append(assay)
                proteinDict[uniprotKey]['Study'].append(studyName)

    return proteinDict

# --- Format and export data ---
def formatAndExport(proteinDict):
    for uniprot in proteinDict:
        proteinDict[uniprot]['Assay'] = '|'.join(proteinDict[uniprot]['Assay'])
        proteinDict[uniprot]['Study'] = '|'.join(proteinDict[uniprot]['Study'])

    resultDf = pd.DataFrame.from_dict(proteinDict, orient='index')
    resultDf.index.name = 'UniprotID'
    resultDf.to_csv(outputIntermediateTsv, sep='\t')

    # Collect unique Uniprot IDs
    uniprotList = []
    for combinedId in resultDf.index:
        ids = combinedId.split('~')
        for uid in ids:
            if uid not in uniprotList:
                uniprotList.append(uid)

    with open(outputUniprotListTxt, 'w') as f:
        for uid in uniprotList:
            f.write(uid + '\n')

# --- Run ---
measuredProteinDict = extractMeasuredProteins()
formatAndExport(measuredProteinDict)
