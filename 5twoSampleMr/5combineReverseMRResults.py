# --- Setup ---
import pandas as pd
import os

# --- Paths ---
projectRoot = './prj_190_viking_somalogic/'

MRdir = os.path.join(projectRoot, 'Scripts/twoSampleMr/3significantResultsMR/')
reverseMRdir = os.path.join(projectRoot, 'Scripts/twoSampleMr/4resultsReverseMR/')
outputDir = os.path.join(projectRoot, 'Scripts/twoSampleMr/MRResults/')

referenceDfPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
annotationPath = os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/openGWASChatGPTAnnotated.xlsx')
finalOutputPath = os.path.join(projectRoot, 'Scripts/twoSampleMr/5combineReverseMRResults.xlsx')

# --- Parameters ---
pThresholdForwardMR = 1.532645e-5
pThresholdReverseMR = 0.01
reverseMRcolumns = [
    "inverse_variance_weighted__fixed_effects__nsnp", "inverse_variance_weighted__fixed_effects__beta",
    "inverse_variance_weighted__fixed_effects__se", "inverse_variance_weighted__fixed_effects__p",
    "wald_ratio_nsnp", "wald_ratio_beta", "wald_ratio_se", "wald_ratio_p"
]

# --- Load reference data ---
medicallyRelevantAnnotation = pd.read_excel(annotationPath)
referenceDf = pd.read_excel(referenceDfPath)

# --- Functions ---
def combineDf(df1, df2):
    target = df2.loc[0, 'exposure'].replace(',', '<><>')
    matchIndex = df1.index[df1['outcome'] == target].tolist()
    if not matchIndex:
        print('No match found for', df2.loc[0, 'exposure'])
        return None
    idx = matchIndex[0]
    for _, row in df2.iterrows():
        for col in ['nsnp', 'beta', 'se', 'p']:
            df1.at[idx, f"{row['method']}_{col}"] = row[col]
    return df1

def fixDfDelimiter(filePath):
    with open(filePath) as f:
        lines = f.readlines()
    header = lines[0].strip()
    rows = []
    for line in lines[1:]:
        idxs = [i for i, c in enumerate(line) if c == ',']
        pipePos = line.find('|')
        if pipePos > -1:
            for i in reversed(idxs):
                if i < pipePos:
                    line = line[:i] + '<><>' + line[i + 1:]
        rows.append(line.strip())
    df = pd.DataFrame([row.split(',') for row in [header] + rows])
    df.columns = df.iloc[0]
    df = df.drop(index=0).reset_index(drop=True)
    return df

def checkPValue(df, idx):
    if df.loc[idx, 'p'] > pThresholdForwardMR:
        return False
    for col in ["inverse_variance_weighted__fixed_effects__p", "wald_ratio_p"]:
        val = df.loc[idx, col]
        if isinstance(val, float) and val <= pThresholdReverseMR:
            return True
    return False

def removeCohorts(df):
    exclude = [
        'Asthma (childhood onset) || id:ebi-a-GCST007800',
        'Asthma (adult onset) || id:ebi-a-GCST007799',
        'Breast cancer || id:ebi-a-GCST007236',
        'Schizophrenia || id:ieu-b-5070'
    ]
    return df[~df['outcome'].isin(exclude)].reset_index(drop=True)

def addHugo(df):
    for idx, row in df.iterrows():
        somamerId = row['exposure'].split('_')[-1][:-4]
        hugoMatch = referenceDf[referenceDf['somamerID'] == somamerId]
        if not hugoMatch.empty:
            df.at[idx, 'exposure'] = hugoMatch.iloc[0]['HUGO']
    return df

def filterMedicallyRelevantResults(df, idx):
    outcome = df.loc[idx, 'outcome'].split(' || id:')[1]
    match = medicallyRelevantAnnotation[medicallyRelevantAnnotation['id.outcome'] == outcome]
    if len(match) == 1:
        return match.iloc[0]['Clinically relevant'] == 'Yes'
    if len(match) == 0:
        print('Not in annotation:', outcome)
        return False
    print('Multiple matches in annotation for', outcome)
    return None

# --- Create merged MR database ---
def createDb():
    for fileName in os.listdir(MRdir):
        if not fileName.endswith('.csv'):
            continue
        MRdf = pd.read_csv(os.path.join(MRdir, fileName))
        for col in reverseMRcolumns:
            MRdf[col] = ''
        matchingFiles = [f for f in os.listdir(reverseMRdir) if fileName in f]
        for matchFile in matchingFiles:
            reverseDf = fixDfDelimiter(os.path.join(reverseMRdir, matchFile))
            MRdf = combineDf(MRdf, reverseDf) or MRdf
        outputFile = os.path.join(outputDir, fileName.replace('.csv', '.xlsx'))
        MRdf.to_excel(outputFile, index=False)

# --- Run full pipeline ---
os.makedirs(outputDir, exist_ok=True)
createDb()

# --- Extract significant reverse MR results ---
finalResults = []
for fileName in os.listdir(outputDir):
    if not fileName.endswith('.xlsx'):
        continue
    df = pd.read_excel(os.path.join(outputDir, fileName))
    for idx in range(len(df)):
        if checkPValue(df, idx) and filterMedicallyRelevantResults(df, idx):
            finalResults.append(df.loc[idx, :])

# --- Final formatting and export ---
finalResults = pd.DataFrame(finalResults).reset_index(drop=True)
finalResults = removeCohorts(finalResults)
finalResults = addHugo(finalResults)
finalResults.to_excel(finalOutputPath, index=False)
