import pandas as pd
import os
import numpy as np

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Directory and File Paths ---
colocResultDirectory = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocResultsReplication')
supplementaryTable1Path = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
mrReplicationTablePath = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/5replicationDf.xlsx')

# --- Load Data ---
supplementary1 = pd.read_excel(supplementaryTable1Path)
mrReplicationDf = pd.read_excel(mrReplicationTablePath)


#Checks if the MR was significant and returns True if it was
def checkIfMRsignificant(filename, pValueThreshold = 0.000472):
    
    exposureName = '_'.join(filename.split('_')[:2])
    outcomeName = '_'.join(filename.split('_')[2:]).split('.')[0]
    
    #Extract the MR row in Stage 1
    MRrow = MRdf[(MRdf['Replication MR Exposure'] == exposureName) & (MRdf['Replication MR Stage 1 Study ID'] == outcomeName)]
    #Convert to series
    MRrow = MRrow.squeeze()

    #Extract the MR row in Stage 2 if not found in Stage 1
    if len(MRrow) == 0:
        MRStage2 = MRdf[(MRdf['Replication MR Exposure'] == exposureName)]
        #Find the index of the row with outcomeName in it
        for i, row in MRStage2.iterrows():
            if outcomeName in row['Discovery MR Study ID'].split('|'):
                MRrow = MRdf.loc[i]
                break

    #Check if the MR was significant in Stage 1
    #Check if not NA
    if MRrow['Replication MR Stage 1 p-value'] == MRrow['Replication MR Stage 1 p-value']:
        if MRrow['Replication MR Stage 1 p-value'] < pValueThreshold:
            return True

    #Check if the MR was significant in Stage 2
    pValueList = str(MRrow['Replication MR Stage 2 p-value']).split('|')
    outcomeIDList = MRrow['Discovery MR Study ID'].split('|')
    
    if outcomeName not in outcomeIDList:
        return False

    #Get the index of the outcome in the list
    outcomeIndex = outcomeIDList.index(outcomeName)

    if float(pValueList[outcomeIndex]) < pValueThreshold:
        return True

    return False



resultDf = pd.DataFrame(columns=['Gene name', 'exposure', 'outcome', 'nsnps', 'H0', 'H1', 'H2', 'H3', 'H4'])

for filename in os.listdir(resultDir):

    MRsignificant = 'No'
    #Check if the file was MR significant
    if checkIfMRsignificant(filename):
        MRsignificant = 'Yes'

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

    #Add the gene, exposure and outcome names to the dictionary
    colocResultDict['Gene name'] = filename.split('.')[0].split('_')[0]
    colocResultDict['exposure'] = filename.split('.')[0].split('_')[1]
    colocResultDict['outcome'] = ('_').join(filename.split('.')[0].split('_')[2:])

    #Add the MR significant column
    colocResultDict['MR significant'] = MRsignificant

    #Concat the dictionary to the result dataframe
    resultDf = pd.concat([resultDf, pd.DataFrame([colocResultDict])], ignore_index=True)


#Convert all columns to numeric
resultDf = resultDf.apply(pd.to_numeric, errors='ignore')

#Convert nsnps to integer
resultDf['nsnps'] = resultDf['nsnps'].astype(int)

#Sort by H4
resultDf = resultDf.sort_values('H4')

#Save the result dataframe
resultDf.to_excel(os.path.join(colocResultDirectory, 'colocReplicationResultDf.xlsx'), index=False)
