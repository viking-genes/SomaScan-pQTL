import pandas as pd
import json
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Load Main pQTL Table ---
df = os.path.join(projectRoot, 'Scripts/Post_GWAS_transform/novelpQTL3.xlsx')
df = pd.read_excel(df)

# --- Define Input Paths ---
alreadyPublishedDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/measuredProteins')

VEPDf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/LDproxyVEP.xlsx')
VEPDf = pd.read_excel(VEPDf)

VEPDfSeverity = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/VEPannotation/VEP_consequenceSeverityDescription_20230202.txt')
VEPDfSeverity = pd.read_csv(VEPDfSeverity, sep='\t')

gnomADAlleleFreq = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/gnomADAlleleFrequency/myVariant.info_gnomADeurNonFinAlleleFreq20230210.json')

noLDdirectory = os.path.join(projectRoot, 'GWAS/GWAShits')


def addExtraColumns(df, noLDdirectory):
    #Create new columns
    newColumns = ['a1', 'a0', 'n', 'freq1', 'beta1', 'se']
    for i in newColumns:
        df[i] = ''

    for i in range(len(df)):
        somamerID = df.loc[i, 'somamerID']
        tempDf = pd.read_csv(noLDdirectory + '/' + somamerID + '.csv')
        index = tempDf.index[tempDf['rsid'] == df.loc[i, 'SNP']].tolist()[0]

        for j in newColumns:
            df.at[i, j] = tempDf.loc[index, j]
    
    return df

#Reordering columns to follow A. Gilly's 2021 Cell paper
def reorderCols(df, columnName, newPosition):
    columns = df.columns.to_list()
    oldPosition = columns.index(columnName)
    columns.insert(newPosition,columns.pop(oldPosition))
    df = df[columns]
    return df

#Checks Uniprot names against published Somalogic and Olink assays in large scale studies
def addNovelty(df, alreadyPublishedDir):

    #Create a dictionary of already published proteins in the form of {study: [proteinList]}
    def createPublishedProteinDatabase(alreadyPublishedDir):
        import re

        result = {}

        for assay in ['Somalogic', 'Olink']:
            for study in os.listdir(alreadyPublishedDir + '/' + assay):

                studyName = study.split('.')[0][:-11]
                result[studyName] = []

                studyDf = pd.read_excel(os.path.join(alreadyPublishedDir, assay, study))
                for i in studyDf.index:
                    uniprotID = studyDf.loc[i, 'UniprotID']
                    if uniprotID == uniprotID:
                        #Check if it is does not start with HCE (hybridization control elution) used in Somalogic
                        if uniprotID[:3] != 'HCE':
                            result[studyName].append(uniprotID)

                #Split all protein names in the list by , | ; . or space in case there are multiple proteins in a single cell
                result[studyName] = [re.split(',| |\\||;|\\.', i) for i in result[studyName]]
                result[studyName] = [item for sublist in result[studyName] for item in sublist]
     
                #Leave only unique proteins
                result[studyName] = list(set(result[studyName]))

                #Remove empty strings and NaNs
                result[studyName] = [i for i in result[studyName] if i != '' and i != 'None' and i != 'Family']

                #Check if any uniprotID are not length 6
                for i in result[studyName]:
                    if len(i) != 6:
                        print('UniprotID length is not 6:', i, studyName)
                        exit()

        return result

    #Add studies that do not have a measured protein list but state which panels were used to the result dictionary
    def addOlinkStudies(result, sun2023Df):
        
        #Could not find panel data online so creating from local ORCADES data
        def extractPanelData(panelNameList):

            # --- Define Olink Data Directory ---
            OlinkDataDir = os.path.join(
                projectRoot,
                '/llod_data'
            )
            Olink800Df = pd.read_excel(os.path.join(OlinkDataDir, '20180621_Wilson_NPX.xlsx'))
            OlinkInfDf = pd.read_excel(os.path.join(OlinkDataDir, '2015132_Wilson_INF_NPX.xlsx'))
            OlinkCVD2Df = pd.read_excel(os.path.join(OlinkDataDir, '2015132_Wilson_CVD2_NPX.xlsx'))
            OlinkCVD3Df = pd.read_excel(os.path.join(OlinkDataDir, '2015132_Wilson_CVD3_NPX.xlsx'))

            result = {}
            for i in panelNameList:
                result[i] = {}

            for df in [Olink800Df, OlinkInfDf, OlinkCVD2Df, OlinkCVD3Df]: 
                #Truncate the protein measurement data leaving only the first 5 rows
                df = df.iloc[:5]
                #Transpose the dataframe
                df = df.T
                df.columns = df.iloc[0]
                df = df[1:]

                panelDict = {'MET': 'Olink METABOLISM(v.3402)', 'NEU': 'Olink NEUROLOGY(v.8011)', 'NEX': 'Olink NEURO EXPLORATORY(v.3901)', 'CVD3': 'Olink CARDIOVASCULAR III(v.6001)', 'CVD2': 'Olink CARDIOVASCULAR II(v.5001)', 'INF': 'Olink INFLAMMATION(v.3001)'}

                for i in panelNameList:
                    panelName = panelDict[i]
                    panelDf = df.loc[df.index[df['Panel'] == panelName].tolist()]
                    
                    if len(panelDf) > 0:
                        result[i] = panelDf['Uniprot ID'].to_list()

            return result

        #Retrieve all proteins measured in an Olink panel
        def retrievePanelData(panelDict, panelNameList):
            result = []
            for i in panelNameList:
                result += panelDict[i]

            #Remove nan values
            result = [i for i in result if i == i]

            #Split by , in case there are multiple proteins in a single cell
            result = [i.split(',') for i in result]
            result = [item for sublist in result for item in sublist]

            #Remove duplicates
            result = list(set(result))

            return result

        #Gudjonsson used the same panel as Emilsson
        result['Gudjonsson2022'] = result['Emilsson2022']

        #Add Olink studies using our internal ORCADES panel measurement as reference
        # Bretherick2020 - CVD2, CVD3, INF (249 Proteins)
        # Carland2023 - MET (92 Proteins)
        # Dunlop2021 - CVD2, CVD3 (184 Proteins)
        # Repetto2023 - NEU, NEX (184 Proteins)
        # Zhao2023 - INF (91 Proteins)
        olinkPanelDict = extractPanelData(['CVD2', 'CVD3', 'INF', 'MET', 'NEU', 'NEX'])

        result['Bretherick2020'] = retrievePanelData(olinkPanelDict, ['CVD2', 'CVD3', 'INF'])
        result['Carland2023'] = retrievePanelData(olinkPanelDict, ['MET'])
        result['Dunlop2021'] = retrievePanelData(olinkPanelDict, ['CVD2', 'CVD3'])
        result['Repetto2023'] = retrievePanelData(olinkPanelDict, ['NEU', 'NEX'])
        result['Zhao2023'] = retrievePanelData(olinkPanelDict, ['INF'])

        return result


    df['Protein measured before'] = ''
    df['Protein measured in study'] = ''

    publishedProteinDatabase = createPublishedProteinDatabase(alreadyPublishedDir)
    publishedProteinDatabase = addOlinkStudies(publishedProteinDatabase, os.path.join(alreadyPublishedDir, 'Olink', 'Sun2023ProteinInfo.xlsx'))

    assayStudyDict = {'Olink': ['Sun2023', 'Carland2023', 'Zhao2023', 'Bretherick2020', 'Dunlop2021', 'Repetto2023'], 'Somalogic': ['Emilsson2022', 'Pietzner2020', 'Ferkingstad2021', 'Folkersen2020', 'Gudjonsson2022']}

    for i in df.index:
        uniprotIDList = df.loc[i, 'Uniprot'].split('|')

        #If there are multiple targets, split by | and check each one. Assign novel if at least one of the targets has not been measured before
        for uniprotID in uniprotIDList:
            measured = False
            for study in publishedProteinDatabase.keys():
                if uniprotID in publishedProteinDatabase[study]:
                    measured = True
                    #Find key in assayStudyDict that contains the study name
                    for assay in assayStudyDict.keys():
                        if study in assayStudyDict[assay]:
                            
                            if assay not in df.loc[i, 'Protein measured before']:
                                df.at[i, 'Protein measured before'] += assay + '|'

                            df.at[i, 'Protein measured in study'] += study + '|'
                            break
            
            if measured == False:
                df.at[i, 'Protein measured before'] = ''
                df.at[i, 'Protein measured in study'] = ''
                break

        #Remove trailing |
        if df.loc[i, 'Protein measured before'] != '':
            df.at[i, 'Protein measured before'] = df.loc[i, 'Protein measured before'][:-1]
            df.at[i, 'Protein measured in study'] = df.loc[i, 'Protein measured in study'][:-1]
        
        #If there are no proteins measured before, put NA
        else:
            df.at[i, 'Protein measured before'] = 'No'
            df.at[i, 'Protein measured in study'] = 'NA'
            
        #Order the assays alphabetically and leave only unique values
        df.at[i, 'Protein measured before'] = '|'.join(sorted(df.loc[i, 'Protein measured before'].split('|')))
        df.at[i, 'Protein measured in study'] = '|'.join(sorted(list(set(df.loc[i, 'Protein measured in study'].split('|')))))
        
    return df

def addVEP(df, VEPDf, VEPDfSeverity):

    def getMostSevereConsequence(VEPDf, indexList, VEPDfSeverity):
        #Severity dictionary
        severity = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}

        consequences = {}

        for j in indexList:
            consequences[VEPDf.loc[j, 'VEP']] = VEPDfSeverity.loc[VEPDfSeverity.index[VEPDfSeverity['SO term'] == VEPDf.loc[j, 'VEP']], 'IMPACT'].to_list()

        #Remove empty entries to deal with variants without annotation
        for j in consequences.copy():
            if len(consequences[j]) == 0:
                consequences.pop(j)

        for j in consequences.keys():
            consequences[j] = severity[consequences[j][0]]

        maxSeverity = max(consequences.values())
        result = [key for key in consequences if consequences[key] == maxSeverity]

        return ' | '.join(result)

    df['Most severe consequence'] = ''

    for i in range(len(df)):
        indexes = VEPDf.index[VEPDf['targetSNP'] == df.loc[i, 'SNP']].tolist()

        #Checks for worst consequence if there are multiple matching entries in the VEP database for a particular SNP
        if len(indexes) > 1:
            df.at[i, 'Most severe consequence'] = getMostSevereConsequence(VEPDf, indexes, VEPDfSeverity)
            continue
        
        if len(indexes) != 0:
            df.at[i, 'Most severe consequence'] = VEPDf.loc[indexes[0], 'VEP']
            continue

        #If there is no match in the VEP database, put NA
        df.at[i, 'Most severe consequence'] = 'NA'

    return df

#Adding gnomAD non-finnish European allele frequencies for all rsids in the dataset. json file is given as input, downloaded from myVariant.info with the gnomad_genome.af.af_nfe field selected
def addgnomADAlleleFreq(df, gnomADjson):
    df['Effect allele freq in non-Finnish Eur'] = ''
    gnomADjson = json.load(open(gnomADjson))

    for i in gnomADjson:
        rsid = i['query']
        if 'gnomad_genome' in i.keys():
            chromosome = i['_id'].split(':')[0][3:]
            position = i['_id'].split(':')[1][2:-3]
            otherAllele = i['_id'].split(':')[1][-3]
            effectAllele = i['_id'].split(':')[1][-1]
            effectAlleleFrequency = i['gnomad_genome']['af']['af_nfe']
        
            matchedRows = df.loc[df.index[df['SNP'] == rsid].tolist()]
            for idx in matchedRows.index:
                
                if str(int(matchedRows.loc[idx, 'CHR'])) == str(int(chromosome)) and str(int(matchedRows.loc[idx, 'BP'])) == str(int(position)) and str(matchedRows.loc[idx, 'a0']) == str(otherAllele) and str(matchedRows.loc[idx, 'a1']) == str(effectAllele):
                    df.at[idx, 'Effect allele freq in non-Finnish Eur'] = effectAlleleFrequency

    #Manually assigning the following rsids as they were not exported to json but are present on gnomAD website as of 20230210.
    #rs200418001 https://gnomad.broadinstitute.org/variant/6-31847995-T-G?dataset=gnomad_r2_1
    #9:140008750_G_C https://gnomad.broadinstitute.org/variant/9-140008750-G-C?dataset=gnomad_r2_1
    #rs9269173 doesn't exist https://gnomad.broadinstitute.org/variant/rs9269173?dataset=gnomad_r2_1

    df.at[df.index[df['SNP'] == 'rs200418001'].tolist()[0], 'Effect allele freq in non-Finnish Eur'] = 0.005089
    df.at[df.index[df['SNP'] == '9:140008750_G_C'].tolist()[0], 'Effect allele freq in non-Finnish Eur'] = 0.6737
    df.at[df.index[df['SNP'] == 'rs9269173'].tolist()[0], 'Effect allele freq in non-Finnish Eur'] = 'NA'

    return df

def convertColumnToAbsoluteValues(df, column):
    for i in range(len(df)):
        try:
            df.at[i, column] = abs(df.loc[i, column])
        except TypeError:
            continue

    return df

def removeMAF(df):
    df = df.drop(df[df.freq1 < 0.05].index)
    df = df.drop(df[df.freq1 > 0.95].index)
    return df

#Many SNPs were found to be in long distance LD with each other after clumping in the 250kb region. This function filters out SNPs whose clumping distances overlap with each other and only keeps the SNP with the lowest P-value
def filterForIndependence(df):
    def checkChromosome(df):
        toDrop = []
        #Check if it is the same chromosome
        for i in range(len(df)):
            if df.iloc[i, 1] != df.iloc[0, 1]:
                toDrop.append(i)

        df = df.drop(df.index[toDrop])
        return df
    
    #Adds a column to the dataframe to keep track of which SNPs have been used in the filtering process
    df['Used'] = False

    result = []
    #List of somamers
    somamers = df['somamerID'].unique()

    #For each somamer, filter out SNPs that are in LD with each other
    for i in somamers:

        #Get all SNPs for the somamer
        somamerDf = df.loc[df.index[df['somamerID'] == i].tolist()]
        #Loops while there are still unused SNPs
        while len(somamerDf.loc[somamerDf.index[somamerDf['Used'] == False].tolist()]) > 0:
            #Get all SNPs for the somamer
            somamerDf = df.loc[df.index[df['somamerID'] == i].tolist()]

            #Leaves only unused SNPs
            somamerDf = somamerDf.loc[somamerDf.index[somamerDf['Used'] == False].tolist()]
            if len(somamerDf) == 0:
                break
            
            #Extracts SNPs that are on the same chromosome as the one in the first row
            somamerDf = checkChromosome(somamerDf)

            #If there is only one SNP left, add it to the result and continue to the next somamer
            if len(somamerDf) == 1:
                result.append(somamerDf)
                break
            
            #For each SNP in the dataframe, check if it overlaps in LD distance with the the next SNP
            somamerDf = somamerDf.sort_values(by = 'Position')
            somamerDfIndexes = somamerDf.index.tolist()

            for j in range(len(somamerDfIndexes) - 1):
                if somamerDf.loc[somamerDfIndexes[j + 1], 'Position'] - somamerDf.loc[somamerDfIndexes[j], 'Position'] >= 500000:
                    somamerDf = somamerDf.drop(somamerDfIndexes[j+1:])
                    break
            
            for j in somamerDf.index:
                df.at[j, 'Used'] = True

            #Get the index with the lowest P-value
            index = somamerDf.index[somamerDf['P-value'] == somamerDf['P-value'].min()].tolist()[0]
            result.append(df.loc[[index]])

    result = pd.concat(result)
    result = result.reset_index(drop = True)

    result.drop('Used', axis = 1, inplace = True)

    return result

def addVarianceExplained(df):
    
    df['Variance Explained'] = ''

    for i in range(len(df)):
        df.at[i, 'Variance Explained'] = 2 * df.loc[i, 'freq1'] * (1 - df.loc[i, 'freq1']) * df.loc[i, 'beta1']**2  

    return df

#Calculate F-statistic for forward MR instrument strength estimation
def addFStatistic(df):
    
    df['F-statistic'] = ''

    #Formula to calculate the F-statistic for a single SNP https://academic.oup.com/ije/article/40/3/755/745918?login=true
    #F=VarianceExplained(Nsamples-2)/(1-VarianceExplained)
    for index in df.index:

        pQTLVarianceExplained = df.loc[index, 'Variance Explained']
        pQTLn = df.loc[index, 'n']

        fStatistic = pQTLVarianceExplained * (pQTLn - 2) / (1 - pQTLVarianceExplained)

        df.at[index, 'F-statistic'] = round(fStatistic, 2)

    return df


#Checks both HUGO and Uniprot names against a published SomaLogic list and Olink 1500 panel
def addNovelty2(df):
    def findMatch(df1, df2, columnName):
        result = df2.index[df2[columnName] == df1.loc[i, columnName]].tolist()
        return result

    # --- Define Path to Previously Published SNPs ---
    alreadyPublishedDf = os.path.join(
        projectRoot,
        'Scripts/SNPnovelty2/1previouslyPublished.xlsx'
    )
    alreadyPublishedDf = pd.read_excel(alreadyPublishedDf)

    df['novel'] = ''

    for i in range(len(df)):

        if len(findMatch(df, alreadyPublishedDf, 'HUGO')) == 0 and len(findMatch(df, alreadyPublishedDf, 'Uniprot')) == 0:
            df.at[i, 'novel'] = 'TRUE'
            continue
        df.at[i, 'novel'] = 'FALSE'

    #Doublecheck for pipes; if a somamer is targeting multiple proteins and both match existing ones, it is not novel

    for i in range(len(df)):
        if df.loc[i, 'novel'] == 'TRUE':
            if '|' in str(df.loc[i, 'HUGO']):
                toMatch = df.loc[i, 'HUGO'].split('|')
                for j in alreadyPublishedDf['HUGO']:
                    for k in toMatch.copy():
                        if k in str(j):
                            toMatch.remove(k)
                if len(toMatch) == 0:
                    df.at[i, 'novel'] = 'FALSE'

    #Then adds an additional column to indicate if the protein is novel or not based on referencing the full somalogic 4k panel
    # --- Load Supplementary Tables ---
    publishedDf = pd.read_excel(os.path.join(
        projectRoot, 'Scripts/SNPnovelty/somalogic4Olink1500Combined.xlsx'))

    somaLogicDfSeqID = pd.read_excel(os.path.join(
        projectRoot, 'Scripts/SNPnovelty/previouslyPublished/SomaLogicv4.0panelFerkingstad2021.xlsx'), skiprows=2)
    somalogicDfSeqId = somaLogicDfSeqID['SeqId'].tolist()

    #In all 5 mismatch cases Uniprot was correct because it seems that SomaLogic changed their gene names

    df['novelNewUniprot'] = ''
    publishedUniprot = publishedDf['UniProt'].tolist()
    publishedUniprot = [i.split('|') for i in publishedUniprot]
    #Expand list of lists
    publishedUniprot = [item for sublist in publishedUniprot for item in sublist]
    publishedUniprot = list(set(publishedUniprot))

    for i in range(len(df)):
        #If multiple Uniprot IDs are present, split them and check if all of them match the published ones
        if '|' in df.loc[i, 'Uniprot']:
            uniprot = df.loc[i, 'Uniprot'].split('|')
            if set(uniprot).issubset(publishedUniprot):
                df.at[i, 'novelNewUniprot'] = False
            else:
                df.at[i, 'novelNewUniprot'] = True
        #If only one Uniprot ID is present, check if it matches the published ones
        else:
            if df.loc[i, 'Uniprot'] in publishedUniprot:
                df.at[i, 'novelNewUniprot'] = False
            else:
                df.at[i, 'novelNewUniprot'] = True
        #Check if the somamerID is present in the v4.0 assay since the annotation changed but the somamerID stayed the same
        seqID = '_'.join(df.loc[i, 'somamerID'].split('-'))
        if seqID in somalogicDfSeqId:
            df.at[i, 'novelNewUniprot'] = False
                
    #Join the two novelty columns together only if both are true
    df['novelNew'] = ''
    for i in range(len(df)):
        if df.loc[i, 'novel'] == 'TRUE' and df.loc[i, 'novelNewUniprot'] == True:
            df.at[i, 'novelNew'] = 'TRUE'
        else:
            df.at[i, 'novelNew'] = 'FALSE'

    #Remove the intermediate columns
    df.drop(columns = ['novelNewUniprot', 'novel'], inplace = True)
    df.rename(columns = {'novelNew': 'novel'}, inplace = True)
    
    return df

df = addExtraColumns(df, noLDdirectory)
df = addNovelty(df, alreadyPublishedDir)
df = addVEP(df, VEPDf, VEPDfSeverity)
df = addgnomADAlleleFreq(df, gnomADAlleleFreq)
df = convertColumnToAbsoluteValues(df, 'Genetic distance')
df = addVarianceExplained(df)
df = addFStatistic(df)
df = addNovelty2(df)

newOrder = ['HUGO', 'Uniprot', 'SNP', 'CHR', 'BP', 'n', 'a1', 'a0', 'freq1', 'Effect allele freq in non-Finnish Eur', 'beta1', 'se', 'F-statistic', 'P', 'Variance Explained', 'Somamer target TSS Chr', 'Somamer target TSS', 'Genetic distance', 'cisTrans', 'Overlapped Gene', 'Nearest Upstream Gene', 'Nearest Downstream Gene', 'Most severe consequence', 'novel', 'Protein measured before', 'Protein measured in study']

for num, i in enumerate(newOrder):
    df = reorderCols(df, i, num)

#Remove columns
removeCols = ['F', 'S05', 'S01', 'S001', 'S0001', 'SP2', 'Frequency', 'TOTAL', 'Effect Size', 'NSIG']

for i in removeCols:
    df = df.drop(i, axis = 1)

#Exclude MAF < 0.05
df = removeMAF(df)

#Rename columns
df = df.rename(columns = {'beta1': 'Effect Size (beta)', 'n': 'Individuals included in the study', 'a1': 'Effect allele', 'a0': 'Other allele', 'freq1': 'Effect allele frequency in VIKING', 'se': 'Standard error of beta', 'BP': 'Position', 'CHR': 'Chromosome', 'P': 'P-value', 'Somamer target TSS Chr': 'Somamer targeted protein chromosome', 'Somamer target TSS': 'Somamer targeted protein TSS', 'Genetic distance': 'Bp distance between targeted protein TSS and SNP'})
            
df = filterForIndependence(df)

df = df.sort_values(by = ['Chromosome', 'Position'])
df = df.reset_index(drop = True)

# --- Define Output Path ---
supplementary1Path = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')

# --- Export DataFrame ---
df.to_excel(supplementary1Path, index=False)

