import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, shapiro

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

supplementaryFilePath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
statisticsFilePath = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/PAVBoxPlot/statistics.txt')
outputDir = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/PAVBoxPlot/')
os.makedirs(outputDir, exist_ok=True)

# --- Load data ---
supplementaryDf = pd.read_excel(supplementaryFilePath)

# --- Function ---
def plotBetaExplainedByPAV(inputDf, cisTrans='All', mafThreshold=None):
    def calculateStatisticalDifference(list1, list2, list1Name, list2Name):
        list1stat, list1p = shapiro(list1)
        list2stat, list2p = shapiro(list2)

        stat, p = mannwhitneyu(list1, list2)

        with open(statisticsFilePath, 'a') as f:
            f.write('All' + list1Name + ' ' + list2Name + '\n')
            f.write('Statistics=%.3f, p=%.3e\n' % (stat, p))
            f.write('Average effect size ' + list1Name + ': %.2f N = %d\n' % (np.mean(list1), len(list1)))
            f.write('Average effect size ' + list2Name + ': %.2f N = %d\n\n' % (np.mean(list2), len(list2)))

        return p

    df = inputDf.copy()

    if mafThreshold is not None:
        df = df[df['Effect allele frequency in VIKING'] > mafThreshold]

    if cisTrans != 'All':
        df = df[df['cisTrans'] == cisTrans]
    df.reset_index(inplace=True, drop=True)

    nonPAV = []
    moderatePAV = []
    highPAV = []

    moderateImpactTerms = ['inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant']
    highImpactTerms = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained',
                       'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification']

    for _, row in df.iterrows():
        consequence = row.get('Most severe consequence')
        if isinstance(consequence, str):
            firstConsequence = consequence.split('|')[0]
            varianceExplained = abs(row.get('Variance Explained', np.nan))
            if firstConsequence in moderateImpactTerms:
                moderatePAV.append(varianceExplained)
            elif firstConsequence in highImpactTerms:
                highPAV.append(varianceExplained)
            else:
                nonPAV.append(varianceExplained)

    allGroups = [highPAV, moderatePAV, nonPAV]
    groupLabels = ['High impact PAV', 'Moderate impact PAV', 'Non-PAV']

    for label, values in zip(groupLabels, allGroups):
        print(f'Average variance explained {label}: {np.mean(values):.2f} N = {len(values)}')

    # --- Plot ---
    plt.figure()
    plt.boxplot(
        allGroups,
        labels=[f'{label}\nN = {len(vals)}' for label, vals in zip(groupLabels, allGroups)]
    )
    plt.scatter([1]*len(highPAV), highPAV, color='black', facecolors='none')
    plt.ylabel('GWAS Variance Explained')
    plt.xlabel(f'{cisTrans} Protein Altering Variants')

    for i in range(len(allGroups)):
        for j in range(i + 1, len(allGroups)):
            pValue = calculateStatisticalDifference(allGroups[i], allGroups[j], groupLabels[i], groupLabels[j])
            if pValue < 0.05:
                y = max(max(allGroups[i]), max(allGroups[j])) + 0.04
                plt.plot([i + 1, j + 1], [y, y], color='black')
                plt.plot([i + 1, i + 1], [y, y - 0.02], color='black')
                plt.plot([j + 1, j + 1], [y, y - 0.02], color='black')
                plt.text((i + j + 2) / 2, y, '*', ha='center', va='center', fontsize=20)
                plt.ylim(0, y + 0.06)

    # --- Save outputs ---
    fileBase = f'VarianceExplained{cisTrans}'
    plt.savefig(os.path.join(outputDir, fileBase + '.png'))
    plt.savefig(os.path.join(outputDir, fileBase + '.svg'))

# --- Run ---
plotBetaExplainedByPAV(supplementaryDf, cisTrans='All', mafThreshold=None)
