import pandas as pd
import os
import ieugwaspy
import sys

# --- Define project root and paths ---
projectRoot = './prj_190_viking_somalogic/'

MRresultDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/2resultsMR')
preMRdir = os.path.join(MRresultDir, 'harmonisedPreMRData')
vikingMRpath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryMR.xlsx')

# --- Load VIKING MR results ---
VIKINGMRdf = pd.read_excel(vikingMRpath)


outcomeMatchDict = {'Diastolic_Blood_Pressure': ['ieu-b-39', 'ukb-b-7992'], 'Forced_Vital_Capacity': ['ieu-b-105', 'ukb-b-7953'], 'Systolic_Blood_Pressure': ['ieu-b-38'], 'eLife_Maternal_longevity': ['ebi-a-GCST006696'], 'GIANT_HEIGHT': ['ieu-a-89'], 'Educational_Attainment': ['ieu-a-1239'], 'ieu-a-1006': ['ebi-a-GCST004599'], 'ieu-b-4809': ['ieu-b-85'], 'ebi-a-GCST90093322': ['ebi-a-GCST006097'], 'GIANT_BMI': ['ukb-b-8909', 'ukb-b-19953', 'ukb-b-12854', 'ukb-b-20531'], 'type_2_diabetes': ['ebi-a-GCST006867', 'ukb-b-10753', 'ukb-b-14609'], 'smoking_cessation': ['ukb-b-19842'], 'ieu-a-302': ['ieu-b-111']}

#Query the OpenGWAS database for study metadata
def queryOpenGWASmetadata(studyID):

    #OpenGWAS API key is saved in .ieugwaspy.json {"jwt": "myKey"}

    #Query the OpenGWAS database for metadata if it's sourced from OpenGWAS
    if '-' in studyID[0]:
        metadata = ieugwaspy.query.gwasinfo(studyID)

        #Fix selected metadata
        if studyID[0] == 'ieu-b-105':
            metadata[0]['pmid'] = 34294062
        if studyID[0] == 'ieu-b-85':
            metadata[0]['sample_size'] = 140254
        if studyID[0] in ['ukb-b-7953', 'ieu-b-4809', 'ukb-b-8909', 'ukb-b-19953', 'ukb-b-12854', 'ukb-b-20531', 'ukb-b-19842', 'ukb-b-7992']:
            metadata[0]['pmid'] = 'NA'

    #If the study is outsourced, use the metadata from manual list
    else:
        metadataDict = {
                        'Diastolic_Blood_Pressure' : pd.DataFrame([{'author': 'Keaton JM', 'year':  2024, 'sample_size': 1028980, 'pmid': 38689001}]),
                        'Forced_Vital_Capacity' : pd.DataFrame([{'author': 'Shrine N', 'year': 2023, 'sample_size': 588452, 'pmid': 36914875}]),
                        'Systolic_Blood_Pressure' : pd.DataFrame([{'author': 'Keaton JM', 'year': 2024, 'sample_size': 1028980, 'pmid': 38689001}]),
                        'eLife_Maternal_longevity' : pd.DataFrame([{'author': 'Timmers P', 'year': 2019, 'sample_size': 502501, 'pmid': 30642433}]),
                        'GIANT_HEIGHT' : pd.DataFrame([{'author': 'Yengo L', 'year': 2022, 'sample_size': 2039467, 'pmid': 36224396}]),
                        'Educational_Attainment' : pd.DataFrame([{'author': 'Okbay A', 'year': 2022, 'sample_size': 3037499, 'pmid': 35361970}]),
                        'GIANT_BMI' : pd.DataFrame([{'author': 'Locke AE', 'year': 2015, 'sample_size': 339224, 'pmid': 25673413}]),
                        'type_2_diabetes' : pd.DataFrame([{'author': 'Suzuki K', 'year': 2024, 'sample_size': 2535601, 'pmid': 38374256}]),
                        'smoking_cessation' : pd.DataFrame([{'author': 'Saunders GRB', 'year': 2022, 'sample_size': 1400535, 'pmid': 36477530}])
        }

        metadata = [metadataDict[studyID[0]]]

    return metadata

#Retrieve all rows from the VIKING MR df that match the outcome
def retrieveVIKINGdata(geneName, outcome):
    result = []

    for outcomeID in outcomeMatchDict[outcome]:
        result.append(VIKINGMRdf[(VIKINGMRdf['Study ID'] == outcomeID) & (VIKINGMRdf['MR Exposure'] == geneName)])
    
    result = pd.concat(result)

    #Drop unnecessary columns
    result = result.drop(columns=['Reverse MR Number of SNPs', 'Reverse MR Effect size (beta)',
       'Reverse MR Standard Error of beta', 'Reverse MR p-value',
       'Reverse MR Method'])

    #Add the MR Study Type column as the first column with Discovery as the value
    result['MR Study Type'] = 'Discovery'
    result = result[['MR Study Type'] + [col for col in result.columns if col != 'MR Study Type']]

    result = result.reset_index(drop=True)

    #Add missing metadata to the VIKING data
    for i in range(len(result)):
        metadata = queryOpenGWASmetadata([result.loc[i, 'Study ID']])
        result.at[i, 'Outcome Study Author'] = metadata[0]['author']
        result.at[i, 'Outcome Study Publication Year'] = metadata[0]['year']
        result.at[i, 'Outcome Study Sample Size'] = metadata[0]['sample_size']
        result.at[i, 'Outcome Study pmid'] = str(metadata[0]['pmid'])

    return result

#Harmonise the MR data
def harmoniseReplicationData(df):

    #Calculate F-statistic for Wald Ratio MR instrument strength estimation
    #Formula to calculate the F-statistic for a single SNP https://academic.oup.com/ije/article/40/3/755/745918?login=true
    #F=VarianceExplained(Nsamples-2)/(1-VarianceExplained)
    def calculateFStatistic(df, singleSNP=False):
        
        #Open the file with the pQTL data for GWAS summary statistics
        exposure = df.loc[0, 'MR Exposure']
        outcomeID = df.loc[0, 'Study ID']
        filename = f'MRHarmonised_{exposure}_{outcomeID}.tsv'
        preMRdf = pd.read_csv(os.path.join(preMRdir, filename), sep='\t')
        
        #Process single SNP data
        if singleSNP:
            #Calculate the variance explained by the exposure
            exposureVarianceExplained = 2 * preMRdf.loc[0, 'eaf.exposure'] * (1 - preMRdf.loc[0, 'eaf.exposure']) * preMRdf.loc[0, 'beta.exposure']**2  
            exposureSampleSize = preMRdf.loc[0, 'samplesize.exposure']

            fStatistic = exposureVarianceExplained * (exposureSampleSize - 2) / (1 - exposureVarianceExplained)

            fStatistic = round(fStatistic, 2)
        
        #Process multiple SNP data
        if not singleSNP:

            FStatisticList = []

            for i in range(len(preMRdf)):
                #Calculate the variance explained by the exposure
                exposureVarianceExplained = 2 * preMRdf.loc[i, 'eaf.exposure'] * (1 - preMRdf.loc[i, 'eaf.exposure']) * preMRdf.loc[i, 'beta.exposure']**2  
                exposureSampleSize = preMRdf.loc[i, 'samplesize.exposure']

                fStatistic = exposureVarianceExplained * (exposureSampleSize - 2) / (1 - exposureVarianceExplained)

                FStatisticList.append(round(fStatistic, 2))
            
            #Return a value of average [min, max]
            fStatistic = f'{round(sum(FStatisticList) / len(FStatisticList), 2)} [{min(FStatisticList)}-{max(FStatisticList)}]'

        return fStatistic

    df.reset_index(drop=True, inplace=True)
    #Leave only necessary columns
    df = df[['id.outcome', 'id.exposure', 'outcome', 'method', 'nsnp.x', 'b', 'se', 'pval', 'year', 'author', 'pmid', 'sample_size']]

    #Rename columns
    df = df.rename(columns={'id.outcome': 'Study ID', 'id.exposure': 'MR Exposure', 'outcome': 'MR Outcome', 'nsnp.x': 'MR Number of SNPs', 'b': 'MR Effect size (beta)', 'se': 'MR Standard Error of beta', 'pval': 'MR p-value', 'year': 'Outcome Study Publication Year', 'author': 'Outcome Study Author', 'pmid': 'Outcome Study pmid', 'sample_size': 'Outcome Study Sample Size', 'method': 'MR Method'})

    #Fix the MR Outcome column
    df['MR Outcome'] = df['MR Outcome'].apply(lambda x: x.split(' || ')[0])

    #Retrieve the metadata for the replication study
    metadata = queryOpenGWASmetadata([df.loc[0, 'Study ID']])
    df['Outcome Study Author'] = metadata[0]['author']
    df['Outcome Study Publication Year'] = metadata[0]['year']
    df['Outcome Study Sample Size'] = metadata[0]['sample_size']
    df['Outcome Study pmid'] = metadata[0]['pmid']
    
    #Add the MR Study Type column as the first column with Replication as the value
    df['MR Study Type'] = 'Replication'
    df = df[['MR Study Type'] + [col for col in df.columns if col != 'MR Study Type']]

    #Add the F-statistic column
    if df.loc[0, 'MR Method'] == 'Wald ratio':
        df.at[0, 'MR Method'] == 'Wald Ratio'
        df.at[0, 'MR Instrumental variable F-statistic'] = calculateFStatistic(df, singleSNP=True)

    # F-statistic for IVW
    if df.loc[0, 'MR Method'] == 'Inverse variance weighted':
        df.at[0, 'MR Method'] == 'IVW'
        df.at[0, 'MR Instrumental variable F-statistic'] = calculateFStatistic(df, singleSNP=False)

    #Round the effect size and p-value 
    df['MR Effect size (beta)'] = df['MR Effect size (beta)'].apply(lambda x: round(x, 5))
    df['MR p-value'] = df['MR p-value'].apply(lambda x: float(f"{x:.5g}"))

    return df

#Places the data in the MR discovery rows into columns of the appropriate replication rows
def separateDiscoveryAndReplication(df):

    df = df.reset_index(drop=True)

    #Create new columns for the discovery data
    df['Discovery MR Outcome'] = ''
    df['Discovery MR Study ID'] = ''
    df['Discovery MR Effect Size (beta)'] = ''
    df['Discovery MR p-value'] = ''
    df['Discovery MR Outcome Study Publication Year'] = ''
    df['Discovery MR Outcome Study Sample Size'] = ''
    df['Discovery MR Outcome Study pmid'] = ''

    #Combine multiple discovery column data into one cell

    #Split the dataframe into chunks of a single outcome
    outcomeList = [[[], []]]
    for i in range(len(df)):
        if i == 0:
            outcomeList[-1][0].append(df.loc[i, 'MR Study Type'])
            outcomeList[-1][1].append(df.loc[i, :])
            continue
        
        if df.loc[i, 'MR Study Type'] == 'Discovery' and 'Replication' in outcomeList[-1][0]:
            outcomeList.append([[], []])
            outcomeList[-1][0].append(df.loc[i, 'MR Study Type'])
            outcomeList[-1][1].append(df.loc[i, :])
        else:
            outcomeList[-1][0].append(df.loc[i, 'MR Study Type'])
            outcomeList[-1][1].append(df.loc[i, :])

    #Make dataframes from the chunks
    for i in range(len(outcomeList)):
        outcomeList[i] = pd.DataFrame(outcomeList[i][1])

    #For each chunk, create a new row for combined Discovery data
    for i in range(len(outcomeList)):
        
        #Create placeholder discovery row
        discoveryRow = pd.DataFrame(columns=outcomeList[i].columns)

        #Find all rows with the discovery data
        discoveryRows = outcomeList[i][outcomeList[i]['MR Study Type'] == 'Discovery'].copy()
        
        #Convert the discovery data to proper formats
        for k in discoveryRows.index:
            discoveryRows.at[k, 'MR Effect size (beta)'] = round(discoveryRows.loc[k, 'MR Effect size (beta)'], 5)
            discoveryRows.at[k, 'Replication MR Stage 1 Standard Error of beta'] = round(discoveryRows.loc[k, 'MR Effect size (beta)'], 5)
            discoveryRows.at[k, 'MR p-value'] = float(f"{discoveryRows.loc[k, 'MR p-value']:.5g}")
            discoveryRows.at[k, 'Outcome Study Publication Year'] = int(discoveryRows.loc[k, 'Outcome Study Publication Year'])
            discoveryRows.at[k, 'Outcome Study Sample Size'] = int(discoveryRows.loc[k, 'Outcome Study Sample Size'])
        discoveryRows['Outcome Study Publication Year'] = discoveryRows['Outcome Study Publication Year'].astype(int)
        discoveryRows['Outcome Study Sample Size'] = discoveryRows['Outcome Study Sample Size'].astype(int)

        #Populate the discovery row
        discoveryRow.loc[0, 'Discovery MR Outcome'] = '|'.join(discoveryRows['MR Outcome'])
        discoveryRow.loc[0, 'Discovery MR Study ID'] = '|'.join(discoveryRows['Study ID'])
        discoveryRow.loc[0, 'Discovery MR Effect Size (beta)'] = '|'.join(discoveryRows['MR Effect size (beta)'].astype(str))
        discoveryRow.loc[0, 'Discovery MR p-value'] = '|'.join(discoveryRows['MR p-value'].astype(str))
        discoveryRow.loc[0, 'Discovery MR Outcome Study Publication Year'] = '|'.join(discoveryRows['Outcome Study Publication Year'].astype(str))
        discoveryRow.loc[0, 'Discovery MR Outcome Study Sample Size'] = '|'.join(discoveryRows['Outcome Study Sample Size'].astype(str))
        discoveryRow.loc[0, 'Discovery MR Outcome Study pmid'] = '|'.join(discoveryRows['Outcome Study pmid'].astype(str))

        #Remove the discovery rows from the chunk
        outcomeList[i] = outcomeList[i][outcomeList[i]['MR Study Type'] == 'Replication']

        #Add the data to the chunk
        outcomeList[i]['Discovery MR Outcome'] = discoveryRow.loc[0, 'Discovery MR Outcome']
        outcomeList[i]['Discovery MR Study ID'] = discoveryRow.loc[0, 'Discovery MR Study ID']
        outcomeList[i]['Discovery MR Effect Size (beta)'] = discoveryRow.loc[0, 'Discovery MR Effect Size (beta)']
        outcomeList[i]['Discovery MR p-value'] = discoveryRow.loc[0, 'Discovery MR p-value']
        outcomeList[i]['Discovery MR Outcome Study Publication Year'] = discoveryRow.loc[0, 'Discovery MR Outcome Study Publication Year']
        outcomeList[i]['Discovery MR Outcome Study Sample Size'] = discoveryRow.loc[0, 'Discovery MR Outcome Study Sample Size']
        outcomeList[i]['Discovery MR Outcome Study pmid'] = discoveryRow.loc[0, 'Discovery MR Outcome Study pmid']

    #Combine the chunks back into one dataframe
    df = pd.concat(outcomeList)
    df.reset_index(drop=True, inplace=True)

    #Rename existing columns
    df.rename(columns={'MR Outcome': 'Replication MR Stage 1 Outcome', 'Study ID': 'Replication MR Stage 1 Study ID', 'MR Exposure': 'Replication MR Exposure', 'MR Number of SNPs': 'Replication MR Number of SNPs', 'MR Effect size (beta)': 'Replication MR Stage 1 Effect size (beta)', 'MR Standard Error of beta': 'Replication MR Stage 1 Standard Error of beta', 'MR p-value': 'Replication MR Stage 1 p-value', 'MR Method': 'Replication MR Stage 1 Method', 'Outcome Study Author': 'Replication Stage 1 Outcome Study Author', 'Outcome Study Publication Year': 'Replication Stage 1 Outcome Study Publication Year', 'Outcome Study Sample Size': 'Replication Stage 1 Outcome Study Sample Size', 'Outcome Study pmid': 'Replication Stage 1 Outcome Study pmid', 'MR Instrumental variable F-statistic': 'Replication MR Stage 1 Instrumental variable F-statistic'}, inplace=True)

    #Remove unnecessary columns
    df.drop(columns=['MR Study Type'], inplace=True)

    return df


# replicationDict = {geneName: {outcome: {VIKINGData: dfRow, replicationData1: dfRow, replicationData2: dfRow}}}
replicationDict = {}

for filename in os.listdir(MRresultDir):
    if filename.endswith('.tsv'):

        df = pd.read_csv(os.path.join(MRresultDir, filename), sep='\t')

        geneName = filename.split('_')[1]
        assay = filename.split('_')[2]
        outcome = ('_').join(filename.split('_')[3:])[:-4]

        #Check if geneName already exists in replicationDict
        if not geneName in replicationDict:
            replicationDict[geneName] = {}
        #Check if outcome already exists in replicationDict
        if not outcome in replicationDict[geneName]:
            replicationDict[geneName][outcome] = {}

        #Retrieve all rows from the VIKING MR df that match the outcome
        replicationDict[geneName][outcome]['VIKINGData'] = retrieveVIKINGdata(geneName, outcome)

        #Add the replication data
        replicationDict[geneName][outcome][assay] = harmoniseReplicationData(df)


#Combine all the dataframes into one
result = pd.DataFrame()
for geneName in replicationDict:
    for outcome in replicationDict[geneName]:
        for MRstudy in replicationDict[geneName][outcome]:
            result = pd.concat([result, replicationDict[geneName][outcome][MRstudy]])

#Reformat the dataframe
result = separateDiscoveryAndReplication(result)

#Reorder the columns
result = result[['Replication MR Exposure', 'Replication MR Stage 1 Outcome', 'Replication MR Stage 1 Effect size (beta)', 'Replication MR Stage 1 Standard Error of beta', 'Replication MR Stage 1 p-value', 'Replication MR Number of SNPs', 'Discovery MR Outcome', 'Discovery MR Effect Size (beta)', 'Discovery MR p-value', 'Replication MR Stage 1 Study ID', 'Replication MR Stage 1 Method', 'Replication MR Stage 1 Instrumental variable F-statistic', 'Replication Stage 1 Outcome Study Author', 'Replication Stage 1 Outcome Study Publication Year', 'Replication Stage 1 Outcome Study Sample Size', 'Replication Stage 1 Outcome Study pmid', 'Discovery MR Study ID', 'Discovery MR Outcome Study Publication Year', 'Discovery MR Outcome Study Sample Size', 'Discovery MR Outcome Study pmid']]

#Save the result
result.to_excel(os.path.join(MRresultDir, '3replicationDf.xlsx'), index=False)
