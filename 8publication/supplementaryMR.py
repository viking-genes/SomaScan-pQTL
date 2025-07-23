import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Input and Output Paths ---
forwardMRdir = os.path.join(projectRoot, 'Scripts/twoSampleMr/3significantResultsMR')
reverseMRdir = os.path.join(projectRoot, 'MRResults')
significantReverseMRdf = os.path.join(projectRoot, 'Scripts/twoSampleMr/5combineReverseMRResults.xlsx')
supplementary1Df = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')


significantReverseMRdf = pd.read_excel(significantReverseMRdf)
supplementary1Df = pd.read_excel(supplementary1Df)

def splitOutcomeColumn(df):
    df['study ID'] = ''
    for i in range(len(df)):
        df.at[i, 'study ID'] = df.loc[i, 'outcome'].split(' || ')[1][3:]
        df.at[i, 'outcome'] = df.loc[i, 'outcome'].split(' || ')[0]
    return df


#Caclulate F-statistic for forward MR
def extractFStatistic(df, supplementary1Df):
    df['MR Instrumental variable F-statistic'] = ''

    for index in df.index:
        
        geneName = df.loc[index, 'MR Exposure']

        #Find index of gene in supplementary1Df
        df1Index = supplementary1Df.index[supplementary1Df['HUGO'] == geneName].tolist()
        
        if len(df1Index) > 1:
            print('Error: Multiple entries for gene in supplementary1Df. Only single SNP F-statistic reference is implemented')
            exit()

        df1Index = df1Index[0]

        df.at[index, 'MR Instrumental variable F-statistic'] = supplementary1Df.loc[df1Index, 'Instrumental variable F-statistic']

    return df

df = splitOutcomeColumn(significantReverseMRdf)
#Remove unnecessary columns
columnsToRemove = ['inverse_variance_weighted_nsnp', 'inverse_variance_weighted_beta', 'inverse_variance_weighted_se', 'inverse_variance_weighted_p']
df = df.drop(columnsToRemove, axis = 1)

#Rename columns
columnsToRename = {'exposure': 'MR Exposure', 'outcome': 'MR Outcome', 'method': 'MR Method', 'study ID': 'Study ID', 'nsnp': 'MR Number of SNPs', 'beta': 'MR Effect size (beta)', 'se': 'MR Standard Error of beta', 'p': 'MR p-value', 'inverse_variance_weighted__fixed_effects__nsnp': 'Reverse MR Number of SNPs', 'inverse_variance_weighted__fixed_effects__beta': 'Reverse MR Effect size (beta)', 'inverse_variance_weighted__fixed_effects__se': 'Reverse MR Standard Error of beta', 'inverse_variance_weighted__fixed_effects__p': 'Reverse MR p-value'}
df = df.rename(columns = columnsToRename)

#Add column for reverse MR study type
df['Reverse MR Method'] = 'Inverse variance weighted fixed effects'
#Assign Wald Ratio if only one SNP was tested
df.loc[df['wald_ratio_nsnp'] == 1, 'Reverse MR Method'] = 'Wald Ratio'
df.loc[df['wald_ratio_nsnp'] == 1, 'Reverse MR Number of SNPs'] = 1
df.loc[df['wald_ratio_nsnp'] == 1, 'Reverse MR Effect size (beta)'] = df.loc[df['wald_ratio_nsnp'] == 1, 'wald_ratio_beta']
df.loc[df['wald_ratio_nsnp'] == 1, 'Reverse MR Standard Error of beta'] = df.loc[df['wald_ratio_nsnp'] == 1, 'wald_ratio_se']
df.loc[df['wald_ratio_nsnp'] == 1, 'Reverse MR p-value'] = df.loc[df['wald_ratio_nsnp'] == 1, 'wald_ratio_p']

#Extract F-statistic for forward MR
df = extractFStatistic(df, supplementary1Df)

#Reorder columns
columns = ['MR Exposure', 'MR Outcome', 'Study ID', 'MR Method', 'MR Number of SNPs', 'MR Instrumental variable F-statistic', 'MR Effect size (beta)', 'MR Standard Error of beta', 'MR p-value', 'Reverse MR Number of SNPs', 'Reverse MR Effect size (beta)', 'Reverse MR Standard Error of beta', 'Reverse MR p-value', 'Reverse MR Method']
df = df[columns]

#Rename MR method 
df.loc[df['MR Method'] == 'wald_ratio', 'MR Method'] = 'Wald Ratio'

# --- Define Output Path ---
outputPath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryMR.xlsx')

# --- Export DataFrame ---
df.to_excel(outputPath, index=False)