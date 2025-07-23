#grch37 - Emilsson, Ferkingstad, Gudjonsson, Pietzner
#grch38 - Sun 2022

import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Directory Paths ---
ourpQTL = pd.read_excel(os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx'))
publishedDir = os.path.join(projectRoot, 'Scripts/SNPnovelty2/previouslyPublished')
outputDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/retrieveWithinLDForPublishedpQTL')

#Extract the pQTLs that are in our data that are non-novel
ourpQTL = ourpQTL[ourpQTL['novel'] == 'No']
ourpQTL = ourpQTL['SNP'].unique().tolist()

#Extract the rsids from the pQTLs
def extract_rsids(df):

    rsidColumnNames = ['SNP', 'rsid']
    chrColumnNames = ['Chr. pQTLs', 'chr\n(var.)', 'Chromosome', 'Chr', 'CHROM']

    #Create a dictionary to store the results
    result = {i:[] for i in range(1, 23)}

    #Attribute the correct chr column name to the dataset
    for i in chrColumnNames:
        if i in df.columns:
            chrColumnNames = i
            break

    #Fix the chr column if it is Gudjonsson's data
    if chrColumnNames == 'chr\n(var.)':
        df[chrColumnNames] = df[chrColumnNames].str.replace('chr', '')


    #Attribute the correct rsid column name to the dataset
    for i in rsidColumnNames:
        if i in df.columns:
            
            #Check if the column is rsXXXXX or chrX:XXXXX
            if df.loc[0, i].startswith('rs'):
                for j in range(len(df)):
                    try:
                        chromosome = int(df.loc[j, chrColumnNames])
                        if chromosome >= 1 and chromosome <=22:
                            #Some have two rsids separated by a comma
                            rsids = df.loc[j, i].split(',')
                            for k in rsids:
                                #Some don't have rsid at all
                                if k == '-':
                                    continue
                                #Transforms chrX:XXXXX to X:XXXXX
                                if k.startswith('chr'):
                                    k = k[3:]
                                result[chromosome].append(k)

                    #Exclude the non-autosomal chromosomes
                    except ValueError:
                        continue

    return result

def createPublishedDf(directory):

    def attributeAndFixUniprot(df):

        uniprotColumnNames = ['UniProt', 'Target UniProt']
        uniprotColumn = ''

        #Attribute the correct uniprot column name to the dataset
        for j in uniprotColumnNames:
            if j in df.columns:
                uniprotColumn = j
                break
        
        #Create empty column if there is no uniprot column
        if uniprotColumn == '':
            df['UniProt'] = ''
            uniprotColumn = 'UniProt'

        #Standardize column name
        df.rename(columns = {uniprotColumn: 'UniProt'}, inplace = True)

        return df, uniprotColumn

    def attributeAndFixSomamerID(df, filename):

        somamerColumnNames = ['SOMAmer', 'SeqId', 'Seq-Id', 'SomaScan.id']
        somamerColumn = ''

        #Attribute the correct somamer column name to the dataset
        for j in somamerColumnNames:
            if j in df.columns:
                somamerColumn = j
                break

        #Create empty column if there is no somamer column
        if somamerColumn == '':
            df['SOMAmer'] = ''
            somamerColumn = 'SOMAmer'

        #Ferkingstad's type is 14151_4
        if filename == 'Ferkingstad2021.xlsx':
            df[somamerColumn] = df[somamerColumn].str.replace('_', '-')

        #Gudjonsson's type is 14151_4_3
        if filename == 'Gudjonsson2022.xlsx':
            df[somamerColumn] = df[somamerColumn].str.replace('_', '-')
            #Split the values by '-' and join the first two values
            df[somamerColumn] = df[somamerColumn].apply(lambda x: '-'.join(x.split('-')[:2]))

        #Pietzner's type is SeqId_14151_4
        if filename == 'Pietzner2021.xlsx':
            df[somamerColumn] = df[somamerColumn].str.replace('_', '-')
            #Split the values by '-' and join the last two values
            df[somamerColumn] = df[somamerColumn].apply(lambda x: '-'.join(x.split('-')[1:]))
        
        #Emmisson's type is 7124-18_3
        if filename == 'Emilsson2020.xlsx':
            df[somamerColumn] = df[somamerColumn].apply(lambda x: x.split('_')[0])

        #Standardize column name
        df.rename(columns = {somamerColumn: 'somamerID'}, inplace = True)

        return df, somamerColumn

    def attributeAndFixRsid(df, filename):

        rsidColumnNames = ['rsid']
        rsidColumn = ''

        #Attribute the correct rsid column name to the dataset
        for j in rsidColumnNames:
            if j in df.columns:
                rsidColumn = j
                break

        #Remove chr from the rsid column
        df[rsidColumn] = df[rsidColumn].str.replace('chr', '')

        #Remove rows where there is no rsid
        df = df[df[rsidColumn] != '-']

        return df, rsidColumn


    #Create an empty dataframe to store the results with the columns: rsid, uniprot, somamer, study
    result = pd.DataFrame(columns = ['rsid', 'UniProt', 'somamerID', 'study'])

    for i in os.listdir(directory):

        if not i.startswith('Olink1536Panel'):
            df = pd.read_excel(os.path.join(directory, i))

            df, uniprotColumn = attributeAndFixUniprot(df)
            df, somamerColumn = attributeAndFixSomamerID(df, i)
            df, rsidColumn = attributeAndFixRsid(df, i)
            df['study'] = i.split('.')[0]

            #Add the df to the result
            result = pd.concat([result, df[['rsid', 'UniProt', 'somamerID', 'study']]])

    #Duplicate pQTL in studies are grouped together, adding their study name to the study column
    def combineDuplicates(df):

        #Fix index so that multiple studies don't overlap and have the same index
        df.reset_index(drop = True, inplace = True)

        for num, i in enumerate(df['rsid']):
            print(num)
            #Checks if the row is not already removed
            if num in df.index.values:
                
                #Checks if the same rsid is reported multiple times
                if len(df[df['rsid'] == i]) > 1:
                    #Then checks if it is reported for the same protein
                    currentUniprot = df.loc[num, 'UniProt']
                    currentSomamer = df.loc[num, 'somamerID']
                    
                    for j in range(num + 1, len(df)):

                        #Checks if the row is not already removed
                        if j in df.index.values:
                            
                        #Finds all matching rows by rsid and either UniProt or somamerID and adds their study name to the study column, then removes the duplicate rows
                            if df.loc[j, 'rsid'] == i:
                                
                                if df.loc[j, 'UniProt'] == currentUniprot or df.loc[j, 'somamerID'] == currentSomamer:
                                    
                                    df.at[num, 'study'] += '|' + df.loc[j, 'study']
                                    df.drop(j, inplace = True)

        df.reset_index(drop = True, inplace = True)
        return df

    result = combineDuplicates(result)
    #Export the result
    result.to_csv(os.path.join(outputDir, 'publishedpQTLs2.csv'), index = False)

    return result


publishedpQTL = {i:[] for i in range(1, 23)}
#Read in the pQTLs from the published papers
for i in os.listdir(publishedDir):

    print(i)
    #Skips the Olink1536Panel file that is just the panel without pQTLs
    if i.startswith('Olink1536Panel'):
        continue

    #Read in the file
    publishedDf = pd.read_excel(os.path.join(publishedDir, i))

    #Extract the rsids
    rsids = extract_rsids(publishedDf)

    #Add the rsids to the dictionary
    for j in rsids:
        publishedpQTL[j].extend(rsids[j])
            
#Leave only unique rsids in the dictionary
for i in publishedpQTL:
    publishedpQTL[i] = list(set(publishedpQTL[i]))

#Export the results by chromosome separated by '\n'
for i in publishedpQTL:
    with open(os.path.join(outputDir, '/publishedRsidByChromosome', 'chr' + str(i) + '.txt'), 'w') as f:
        f.write('\n'.join(publishedpQTL[i]))

createPublishedDf(publishedDir)
