import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
chatGptDir = os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/ChatGPTCategorisation')
outputPath = os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/ChatGPTCategorisation.xlsx')

# --- Function: Parse ChatGPT classification text files ---
def extractChatGptResults(chatGptDir):
    resultsByPass = []

    for fileName in os.listdir(chatGptDir):
        if fileName == 'Prompt.txt':
            continue

        resultForOnePass = []
        with open(os.path.join(chatGptDir, fileName), 'r') as file:
            for line in file:
                tokens = line.strip().split(' ')
                if not tokens:
                    continue
                # Case: separated by " -"
                if len(tokens) > 2 and tokens[-2] == '-':
                    resultForOnePass.append([' '.join(tokens[:-2]), tokens[-1]])
                # Case: separated by ":"
                elif len(tokens) > 1 and tokens[-2].endswith(':'):
                    resultForOnePass.append([' '.join(tokens[:-1])[:-1], tokens[-1]])
                else:
                    print('Error in', fileName, line.strip())
        resultsByPass.append(resultForOnePass)

    # Check for duplicates in each pass
    for passIndex, results in enumerate(resultsByPass):
        seen = set()
        for outcome, _ in results:
            if outcome in seen:
                print(f'Duplicate in pass {passIndex + 1}: {outcome}')
            else:
                seen.add(outcome)

    return resultsByPass

# --- Function: Combine all passes into a single DataFrame ---
def combinePassesToDataFrame(resultsByPass):
    firstPass = resultsByPass[0]
    combinedDf = pd.DataFrame(firstPass, columns=['Outcome', 'ChatGPT Pass1'])

    for i in range(1, len(resultsByPass)):
        columnName = f'ChatGPT Pass{i + 1}'
        combinedDf[columnName] = ''

    for passIndex in range(1, len(resultsByPass)):
        for outcome, value in resultsByPass[passIndex]:
            matchIndex = combinedDf.index[combinedDf['Outcome'] == outcome]
            if not matchIndex.empty:
                combinedDf.loc[matchIndex[0], f'ChatGPT Pass{passIndex + 1}'] = value

    return combinedDf

# --- Function: Append prompt to final DataFrame ---
def attachPromptColumn(combinedDf, chatGptDir):
    promptPath = os.path.join(chatGptDir, 'Prompt.txt')
    with open(promptPath, 'r') as promptFile:
        promptText = promptFile.read()

    combinedDf[''] = ''
    combinedDf['Prompt'] = ''
    combinedDf.at[0, 'Prompt'] = promptText

    return combinedDf

# --- Main Execution ---
rawChatGptResults = extractChatGptResults(chatGptDir)
combinedResultsDf = combinePassesToDataFrame(rawChatGptResults)
finalDf = attachPromptColumn(combinedResultsDf, chatGptDir)
finalDf.to_excel(outputPath, index=False)
