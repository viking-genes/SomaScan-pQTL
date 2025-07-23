import os
import pandas as pd
import shutil
import time

# --- Config ---
projectRoot = './prj_190_viking_somalogic/'
inputDirectory = os.path.join(projectRoot, 'GWAS/GWAShitsMAF0.05/onlySignficantHits/')
plinkReferencePath = './10k_unrelated_white_british_reference/'
processedDirectory = os.path.join(inputDirectory, 'processed')

# --- Separate hits by chromosome ---
def separateByChromosome(inputDirectory, fileName):
    inputPath = os.path.join(inputDirectory, fileName)
    df = pd.read_csv(inputPath)

    # Allocate rows by chromosome into a dictionary
    chrDict = {i: [] for i in range(1, 23)}
    for i, row in df.iterrows():
        currentChr = row['chr']
        chrDict[currentChr].append(i)

    # Create temp folders
    tempDir = os.path.join(inputDirectory, fileName[:-4])
    chrDir = os.path.join(tempDir, 'separatedByChromosome')
    os.makedirs(chrDir, exist_ok=True)

    # Export chromosome-specific CSVs
    for chrNum, indices in chrDict.items():
        if not indices:
            continue
        chrDf = df.loc[indices]
        chrDf.to_csv(os.path.join(chrDir, f'{chrNum}.csv'), sep=' ', index=False)

    return tempDir, chrDir

# --- Create Qsub scripts per chromosome ---
def createQsubFiles(tempDir, chrDir):
    qsubDir = os.path.join(tempDir, 'qsub')
    os.makedirs(qsubDir, exist_ok=True)
    logDir = os.path.join(tempDir, 'logs')
    os.makedirs(logDir, exist_ok=True)

    for fileName in os.listdir(chrDir):
        chrNum = fileName[:-4]  # Remove .csv
        qsubPath = os.path.join(qsubDir, f'chr{chrNum}.sh')

        with open(qsubPath, 'w') as f:
            clumpOut = os.path.join(tempDir, 'clumpResults', fileName)
            plinkCommand = (
                f'plink --bfile {plinkReferencePath}ukbb_chr{chrNum}_10000_random_unrelated_white_british '
                f'--clump {os.path.join(chrDir, fileName)} '
                f'--clump-snp-field rsid --clump-field p --clump-kb 250 --clump-r2 0.01 '
                f'--out {clumpOut} --clump-p2 0.0000025 --clump-p1 0.00000005 --memory 30000'
            )

            f.write(
                f'#!/bin/sh\n'
                f'#$ -l h_vmem=32G\n'
                f'#$ -o {logDir}/$JOB_NAME.o\n'
                f'#$ -e {logDir}/$JOB_NAME.e\n'
                f'#$ -cwd\n'
                f'. /etc/profile.d/modules.sh\n'
                f'module load roslin/plink/1.90p\n'
                f'{plinkCommand}\n'
            )

# --- Submit Qsub jobs ---
def runPlinkClumpJobs(tempDir):
    clumpResultsDir = os.path.join(tempDir, 'clumpResults')
    os.makedirs(clumpResultsDir, exist_ok=True)
    qsubDir = os.path.join(tempDir, 'qsub')

    for scriptName in os.listdir(qsubDir):
        os.system(f'qsub {os.path.join(qsubDir, scriptName)}')

# --- Main Loop ---
def main():
    os.makedirs(processedDirectory, exist_ok=True)
    counter = 0

    for fileName in os.listdir(inputDirectory):
        if not fileName.endswith('.csv'):
            continue

        print(f'Processing {fileName}')
        tempDir, chrDir = separateByChromosome(inputDirectory, fileName)
        createQsubFiles(tempDir, chrDir)
        runPlinkClumpJobs(tempDir)

        # Move processed file
        shutil.move(os.path.join(inputDirectory, fileName), os.path.join(processedDirectory, fileName))

        # Sleep after every 200 files to avoid overloading the queue
        counter += 1
        if counter >= 200:
            counter = 0
            print('Sleeping for 20 minutes to respect job submission limits...')
            time.sleep(1200)


main()
