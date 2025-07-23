import os
import pandas as pd

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

vikingIdPath = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/calculateCommonAlleles/VIKING200SampleFile.txt')
phenotypeTablePath = '/data/base_data/viking/phenotypes/viking_base_phenotypes.tsv'
cohortSummaryExportPath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryCohortSummary.xlsx')


# --- Functions ---
def retrieveVikingIds():
    vikingIds = pd.read_csv(vikingIdPath, sep=' ')
    vikingIds = vikingIds['ID_1'].tolist()[1:]
    return vikingIds

def retrieveCohortTable(idList):
    phenotypeDf = pd.read_csv(phenotypeTablePath, sep='\t')
    filteredDf = phenotypeDf[phenotypeDf['iid'].isin(idList)].reset_index(drop=True)
    return filteredDf

def calculateSummaries(statisticList, cohortDf):
    summaryDict = {key: {'male': None, 'female': None} for key in statisticList}

    maleDf = cohortDf[cohortDf['sex'] == 1]
    femaleDf = cohortDf[cohortDf['sex'] == 2]

    # Remove 'sex' from statistics to prevent conflict
    statisticList.remove('sex')
    summaryDict['sex']['male'] = round(len(maleDf) / len(cohortDf), 3)
    summaryDict['sex']['female'] = round(len(femaleDf) / len(cohortDf), 3)

    for stat in statisticList:
        summaryDict[stat]['male'] = round(maleDf[stat].mean(), 2)
        summaryDict[stat]['female'] = round(femaleDf[stat].mean(), 2)

    return summaryDict

def formatSummariesIntoDf(summaryDict):
    resultDf = pd.DataFrame(columns=summaryDict.keys(), index=['male', 'female'])

    for key in summaryDict:
        resultDf.at['male', key] = summaryDict[key]['male']
        resultDf.at['female', key] = summaryDict[key]['female']

    resultDf.rename(columns={
        'sex': 'Sex',
        'age': 'Age',
        'edu': 'Education Level, years',
        'diabetic': 'Diabetic',
        'bmi': 'Body Mass Index (BMI)',
        'bp_sys': 'Systolic Blood Pressure',
        'bp_dia': 'Diastolic Blood Pressure',
        'bp_med': 'Blood Pressure Medicated',
        'fev1': 'Forced Expiratory Volume (1 sec) (fev1)',
        'tot_chol': 'Total Cholesterol'
    }, inplace=True)

    return resultDf


# --- Run ---
vikingIdList = retrieveVikingIds()
cohortDf = retrieveCohortTable(vikingIdList)

statList = ['sex', 'age', 'edu', 'diabetic', 'bmi', 'bp_sys', 'bp_dia', 'bp_med', 'fev1', 'tot_chol']
cohortSummaryDict = calculateSummaries(statList, cohortDf)
cohortSummaryDf = formatSummariesIntoDf(cohortSummaryDict)

cohortSummaryDf.to_excel(cohortSummaryExportPath)
