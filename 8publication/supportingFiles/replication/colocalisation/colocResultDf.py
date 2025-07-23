import pandas as pd
import os
import numpy as np

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Directory Paths ---
colocResultDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocResults')
aarcDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/AARC_GRCh37')
ukbbDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/UKBB_GRCh37')
eqtlGenDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/eQTLGenSummaryStats/extract')
supplementaryTable1Path = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')

# --- Load Supplementary Table ---
supplementary1 = pd.read_excel(supplementaryTable1Path)

resultDf = pd.DataFrame(columns=['Gene name', 'somamerID', 'comparison assay', 'comparison region max -logp', 'nsnps', 'H0', 'H1', 'H2', 'H3', 'H4'])

for filename in os.listdir(resultDir):

    # Read the file content
    with open(os.path.join(resultDir, filename), 'r') as f:
        colocResult = f.read()

    #Create the result dictionary
    colocResultDict = {'nsnps': '', 'H0': '', 'H1': '', 'H2': '', 'H3': '', 'H4': ''}
    
    #Used in case the output file was split into separate lines
    resultList = []

    #Read lines in reverse order
    for line in colocResult.split('\n')[::-1]:
        
        #Skip empty lines
        if line == '':
            continue
        
        #Split the line by empty space
        lineSplit = line.split(' ')
        lineSplit = [x for x in lineSplit if x != '']
        
        #Check if all data is represented here or was the output file split into separate lines
        if len(lineSplit) == 6:
            colocResultDict['nsnps'] = lineSplit[0]
            colocResultDict['H0'] = lineSplit[1]
            colocResultDict['H1'] = lineSplit[2]
            colocResultDict['H2'] = lineSplit[3]
            colocResultDict['H3'] = lineSplit[4]
            colocResultDict['H4'] = lineSplit[5]
            break
        
        #If not all data is represented here, then add the data in reverse order to the list
        else:
            if len(resultList) != 6:
                #Check if the line is a header by numeric values
                if not lineSplit[0][0].isnumeric():
                    continue

                for i in lineSplit[::-1]:
                    resultList.append(i)

    #Populate the dictionary with the data from the reversed list
    if resultList:
        colocResultDict['nsnps'] = resultList.pop()
        colocResultDict['H0'] = resultList.pop()
        colocResultDict['H1'] = resultList.pop()
        colocResultDict['H2'] = resultList.pop()
        colocResultDict['H3'] = resultList.pop()
        colocResultDict['H4'] = resultList.pop()

    #Add the somamerID to the dictionary
    colocResultDict['somamerID'] = filename.split('.')[0].split('_')[0]

    #Add the comparison assay to the dictionary
    if 'UKBB' in filename:
        colocResultDict['comparison assay'] = 'olink'
    elif 'eQTLGen' in filename:
        colocResultDict['comparison assay'] = 'eQTLGen'
    else:
        colocResultDict['comparison assay'] = 'somalogic'

    #Add the max -logp value to the dictionary
    def referenceSummaryStats(somamerID, assay, olinkDir, eQTLdir, somalogicDir, supplementary1):
        
        if assay == 'olink':
            #Get the uniprotID from the supplementary1 dataframe
            uniprotID = supplementary1[supplementary1['somamerID'] == somamerID]['Uniprot'].values[0]
            #Get the reference dataframe
            referenceDf = pd.read_csv(os.path.join(olinkDir, f'{uniprotID}.tsv'), sep='\t')
            #Get the max -logp value
            return referenceDf['LOG10P'].max()

        elif assay == 'somalogic':
            referenceDf = pd.read_csv(os.path.join(somalogicDir, f'{somamerID}.tsv'), sep='\t')
            #Get the min p-value
            minp = referenceDf['p_value'].min()
            #Convert it to -log10
            minp = -np.log10(minp)
            return minp

        elif assay == 'eQTLGen':
            #Get the gene name from the supplementary1 dataframe
            geneName = supplementary1[supplementary1['somamerID'] == somamerID]['HUGO'].values[0]
            #Get the reference dataframe
            referenceDf = pd.read_csv(os.path.join(eQTLdir, f'{geneName}.tsv'), sep='\t')
            #Get the min p-value
            minp = referenceDf['Pvalue'].min()
            #Convert it to -log10
            minp = -np.log10(minp)
            return minp

    pValue = referenceSummaryStats(colocResultDict['somamerID'], colocResultDict['comparison assay'], UKBBDir, eQTLDir, AARCDir, supplementary1)
    colocResultDict['comparison region max -logp'] = round(pValue, 1)

    #Add the gene name to the dictionary
    colocResultDict['Gene name'] = supplementary1[supplementary1['somamerID'] == colocResultDict['somamerID']]['HUGO'].values[0]

    #Concat the dictionary to the result dataframe
    resultDf = pd.concat([resultDf, pd.DataFrame([colocResultDict])], ignore_index=True)


#Convert all columns to numeric
resultDf = resultDf.apply(pd.to_numeric, errors='ignore')

#Convert nsnps to integer
resultDf['nsnps'] = resultDf['nsnps'].astype(int)

#Sort by H4
resultDf = resultDf.sort_values('H4')

#Save the result dataframe
resultDf.to_excel(os.path.join(colocResultDirectory, 'colocResultDf.xlsx'), index=False)

