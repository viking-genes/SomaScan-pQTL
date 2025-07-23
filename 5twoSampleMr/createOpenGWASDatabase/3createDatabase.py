import pandas as pd

# Read in the data
df = pd.read_csv('/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_190_viking_somalogic/Scripts/twoSampleMr/createOpenGWASDatabase/startingDataset.tsv', sep='\t')
chatGPTDf = pd.read_excel('/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_190_viking_somalogic/Scripts/twoSampleMr/createOpenGWASDatabase/ChatGPTCategorisation.xlsx')

#outcome_name column contains the lowercase name of the trait. This will be used to extract unique studies

#Removes studies that have shown to be unreliable due to strong associations with multiple traits when doing MR
def removeFaultyStudies(df):

    cohortsToRemove = ['ebi-a-GCST007800', 'ebi-a-GCST007799', 'ebi-a-GCST007236', 'ieu-b-5070']
    for cohort in cohortsToRemove:
        df = df[df['id.outcome'] != cohort]

    df.reset_index(drop=True, inplace=True)
    return df

#Removes studies that double as one side and the other, e.g. Leg predicted mass (left) and (right)
def removeSide(df, side):

    # Add brackets to side as it is in brackets in the outcome_name column
    side = '(' + side + ')'

    # All outcomes
    outcomes = df['outcome_name'].tolist()

    for i in range(len(df)):
        if side in df['outcome_name'][i]:

            outcomeName = ' '.join(df['outcome_name'][i].split(' ')[:-1])

            #Check if there is a matching opposite side
            if outcomeName + ' (left)' not in outcomes or outcomeName + ' (right)' not in outcomes:
                continue

            #Remove selected side from outcomes
            df.drop(i, inplace=True)

    df.reset_index(drop=True, inplace=True)
    return df

#Alters the outcome_name column to match the name of related studies. e.g. treatment/medication code: vitamin c product and vitamin c
def alterOutcomeName(df, nameList):

    replacementDict = {
    "high density lipoprotein cholesterol levels": "hdl cholesterol",
	"low density lipoprotein cholesterol levels": "ldl cholesterol",
	"total cholesterol levels": "total cholesterol",
	"cholesterol, total": "total cholesterol",
	"triglyceride levels": "triglycerides",
	"c-reactive protein levels": "c-reactive protein",
	"c-reactive protein level": "c-reactive protein",
	"body fat percentage": "body fat",
	"eosinophil counts": "eosinophil cell count",
	"monocyte count": "monocyte cell count",
	"neutrophil count": "neutrophil cell count",
	"neurociticism": "neuroticism",
	"peak expiratory flow (pef)": "peak expiratory flow",
	"urate levels": "urate",
    "serum 25-hydroxyvitamin d levels": "vitamin d levels",
    "25 hydroxyvitamin d level": "vitamin d levels",
    "sepsis (28 day death in critical care)": "sepsis (28 day death)",
    "parental longevity (father's age at death)": "father's age at death",
    "parental longevity (mother's age at death)": "mother's age at death",
    "parental longevity (combined parental age at death)": "parents' age at death",
    "medication for cholesterol, blood pressure, diabetes, or take exogenous hormones: insulin": "medication for cholesterol, blood pressure or diabetes: insulin",
    "treatment/medication code: insulin product": "medication for cholesterol, blood pressure or diabetes: insulin",
    "blood clot, dvt, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: asthma": "doctor diagnosed asthma",
    "diagnoses - main icd10: j45.9 asthma, unspecified": "doctor diagnosed asthma",
    "age asthma diagnosed by doctor": "age asthma diagnosed",
    "alcohol consumption": "alcohol consumed",
    "coffee consumed": "coffee intake",
    "filtered coffee intake": "coffee intake",
    "instant coffee intake": "coffee intake",
    "hemoglobin a1c levels": "glycated hemoglobin levels",
    "hba1c": "glycated hemoglobin levels",
    "subjective well-being": "subjective well being",
    "waist-hip ratio": "waist-to-hip ratio",
    }

    #Replace all treatment speciality of consultant (recoded): XXXX into main speciality of consultant (recoded): XXXX because they are all just duplicates and will be filtered out in the next step
    for i in range(len(df)):
        if df['outcome_name'][i].startswith('treatment '):
            df.at[i, 'outcome_name'] = df['outcome_name'][i].replace('treatment ', 'main ')

    #Make all outcome names lowercase
    df['outcome_name'] = df['outcome_name'].str.lower()

    #Replace outcome names with synonyms
    for key, value in replacementDict.items():
        df['outcome_name'] = df['outcome_name'].replace(key, value)

    #Replace outcome names matching those in nameList
    for i in range(len(df)):
        for name in nameList:
            if name in df['outcome_name'][i] and df['outcome_name'][i] not in replacementDict.values():
                df.at[i, 'outcome_name'] = name

    return df

#Removes ouotcomes with none of the above in the outcome_name column as they are not useful
def removeNoneOfTheAbove(df):

    for i in range(len(df)):
        if 'none of the above' in df['outcome_name'][i]:
            df.drop(i, inplace=True)
    
    df.reset_index(drop=True, inplace=True)

    return df

#Removes gut microbiome studies if they are duplicate with prevalence and abundance studies
def removeGutMicrobiomeDuplicates(df):

    for i in range(len(df)):
        if df['outcome_name'][i].endswith('prevalence'):
            if df['outcome_name'][i][:-10] + 'abundance' in df['outcome_name'].tolist():
                df.drop(i, inplace=True)
                
    df.reset_index(drop=True, inplace=True)

    return df

def renameOutcomeName(df, nameDict):

    for name in nameDict.keys():
        for i in range(len(df)):
            if name in df['outcome_name'][i]:
                df.at[i, 'outcome_name'] = nameDict[name]

    return df

#Leaves only one row per outcome_name by removing duplicate studies. European study with the largest case or sample size is selected
def removeDuplicates(df):

    #Returns an index list of rows in the dataframe that are of European origin or returns all of them if there are no european cohorts for that trait
    def filterOutNonEuropean(df, indexList):

        #Select non-european cohorts
        selectIndexes = []
        for index in indexList:
            if df['population'][index] != 'European':
                selectIndexes.append(index)

        #Check if any european cohorts remain. If not, check for mixed cohorts
        if len(selectIndexes) == len(indexList):
            selectIndexes = []
            for index in indexList:
                if df['population'][index] != 'Mixed':
                    selectIndexes.append(index)
            
            #If no european or mixed cohorts remain, return all cohorts
            if len(selectIndexes) == len(indexList):
                return indexList
        
        #Remove non-european cohorts
        indexList = [index for index in indexList if index not in selectIndexes]
        return indexList

    #Checks if there is a case size for any of the rows in the index list
    def checkCase(df, indexList):
        
        for i in indexList:
            if df['ncase'][i] == df['ncase'][i]:
                return True
        
        return False

    #Find all rows with duplicate outcome_name
    duplicates = {}
    for i in range(len(df)):
        if df['outcome_name'][i] not in duplicates:
            duplicates[df['outcome_name'][i]] = [i]
        else:
            duplicates[df['outcome_name'][i]].append(i)

    for i in duplicates.copy():
        if len(duplicates[i]) == 1:
            duplicates.pop(i)

    #Remove duplicates
    for key, value in duplicates.items():
        #Filter out non-european cohorts or return all if there are no european cohorts
        duplicates[key] = filterOutNonEuropean(df, value)
        #Remove all non-matching studies
        for i in value:
            if i not in duplicates[key]:
                df.drop(i, inplace=True)

        #Of the remaining cohorts, select the one with the largest case or sample size
        largestCaseSize = 0
        largestCaseSizeIndex = 0

        #Check if there are any cohorts with a case size
        if checkCase(df, duplicates[key]):
            
            #Select the cohort with the largest case size
            for i in duplicates[key]:
                if df['ncase'][i] > largestCaseSize:
                    largestCaseSize = df['ncase'][i]
                    largestCaseSizeIndex = i
        
        else:  
            #Select the cohort with the largest sample size
            for i in duplicates[key]:
                if df['sample_size'][i] > largestCaseSize:
                    largestCaseSize = df['sample_size'][i]
                    largestCaseSizeIndex = i
        
        #Remove all other cohorts from the dataframe
        for i in duplicates[key]:
            if i != largestCaseSizeIndex:
                df.drop(i, inplace=True)

    df.reset_index(drop=True, inplace=True)
    return df


#Assigns the clinically relevant column from the chatGPTDf to the df. TRUE = 1, AMBIGUOUS = 0.5, FALSE = 0. If score >= 2.5, clinically relevant is yes. If score < 2.5, clinically relevant is no
def assignChatGPTAnnotation(df, chatGPTDf, threshold = 2.5):
    def calculateScore(chatGPTDf):
        chatGPTDf['Score'] = ''

        for i in range(len(chatGPTDf)):
            score = 0

            for j in range(1, 6):
                if chatGPTDf['ChatGPT Pass' + str(j)][i] == 'TRUE':
                    score += 1
                elif chatGPTDf['ChatGPT Pass' + str(j)][i] == 'AMBIGUOUS':
                    score += 0.5

            chatGPTDf.at[i, 'Score'] = score

        return chatGPTDf

    df['Clinically relevant'] = ''
    chatGPTDf = calculateScore(chatGPTDf)

    for i in range(len(df)):
        for j in range(len(chatGPTDf)):
            if df['outcome_name'][i].lower() == chatGPTDf['Outcome'][j] or df['trait'][i].lower() == chatGPTDf['Outcome'][j]:
                if chatGPTDf['Score'][j] >= threshold:
                    df.at[i, 'Clinically relevant'] = 'Yes'
                else:
                    df.at[i, 'Clinically relevant'] = 'No'

    return df

def modifyChatGPTAnnotation(df, listToTrue = [], listToFalse = []):
    for i in listToTrue:
        for j in range(len(df)):
            if i == df['outcome_name'][j] or i == df['trait'][j]:
                df.at[j, 'Clinically relevant'] = 'Yes'
                break
    
    for i in listToFalse:
        for j in range(len(df)):
            if i == df['outcome_name'][j] or i == df['trait'][j]:
                df.at[j, 'Clinically relevant'] = 'No'
                break

    return df
    
# Extract unique studies
df = removeFaultyStudies(df)
#Remove a side from the outcomes if both sides are present
df = removeSide(df, 'left')
#Remove none of the above outcomes
df = removeNoneOfTheAbove(df)
#Remove gut microbiome prevalence outcomes if there is also a gut microbiome abundance outcome
df = removeGutMicrobiomeDuplicates(df)
#Alter outcome names to match related studies
df = alterOutcomeName(df, ['vitamin c', 'vitamin d', 'vitamin b12', 'multivitamin', 'vitamin e', 'cheese', 'cholesterol lowering medication', 'blood pressure medication'])
#Rename outcome names if they all match the same term and can be removed in a later step
df = renameOutcomeName(df, {})

#Remove duplicate studies by name selecting the largest case or sample size
df = removeDuplicates(df)
#Assign clinically relevant column from chatGPTDf
df = assignChatGPTAnnotation(df, chatGPTDf)
#Modify clinically relevant column to true for all outcomes in the provided list
df = modifyChatGPTAnnotation(df, ['accelerometer-based physical activity measurement (average acceleration)', 'treatment speciality of consultant (recoded): rheumatology', 'fi1 : numeric addition test', 'main speciality of consultant (recoded): ophthalmology', 'extraversion', 'risk taking', 'years of schooling', 'number of children', 'methods of admission to hospital (recoded): maternity admission: post-partum', 'treatment speciality of consultant (recoded): cardiothoracic surgery', 'main speciality of consultant (recoded): neurosurgery', 'main speciality of consultant (recoded): gynaecology', 'treatment speciality of consultant (recoded): cardiac surgery', 'treatment speciality of consultant (recoded): medical oncology', 'main speciality of consultant (recoded): endocrinology', 'treatment speciality of consultant (recoded): upper gastrointestinal surgery', 'main speciality of consultant (recoded): radiology', 'treatment speciality of consultant (recoded): neurology', 'schizophrenia vs adhd (ordinary least squares (ols))', 'main speciality of consultant (recoded): obstetrics', 'father\'s age', 'main speciality of consultant (recoded): adult mental illness', 'treatment speciality of consultant (recoded): clinical haematology', 'mother\'s age', 'main speciality of consultant (recoded): gastroenterology', 'treatment speciality of consultant (recoded): cardiology', 'treatment speciality of consultant (recoded): infectious diseases', 'main speciality of consultant (recoded): neurology', 'intelligence'], [])


df.to_excel('/gpfs/igmmfs01/eddie/wilson-lab/projects/prj_190_viking_somalogic/Scripts/twoSampleMr/createOpenGWASDatabase/openGWASChatGPTAnnotated.xlsx', index=False)