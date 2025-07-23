import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
gwasHitDir = os.path.join(projectRoot, 'GWAS/GWAShits/')
outputDir = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/')
supplementaryPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')

os.makedirs(outputDir, exist_ok=True)

# --- Load supplementary annotation ---
supplementaryDf = pd.read_excel(supplementaryPath)
uniqueSomamerIDs = supplementaryDf['somamerID'].unique()

# --- Aggregate GWAS results with annotations ---
combinedDf = pd.DataFrame()
for somamerID in uniqueSomamerIDs:
    filePath = os.path.join(gwasHitDir, somamerID + '.csv')
    if not os.path.exists(filePath):
        continue

    gwasDf = pd.read_csv(filePath)
    gwasDf['novel'] = supplementaryDf.loc[supplementaryDf['somamerID'] == somamerID, 'novel'].values[0]
    gwasDf['cisTrans'] = supplementaryDf.loc[supplementaryDf['somamerID'] == somamerID, 'cisTrans'].values[0]

    combinedDf = pd.concat([combinedDf, gwasDf])

# Filter on MAF
combinedDf = combinedDf[(combinedDf['freq1'] >= 0.05) & (combinedDf['freq1'] <= 0.95)]
combinedDf = combinedDf.sort_values(by='chr').reset_index(drop=True)

# --- Plotting function ---
def plotMiami(df):
    chromosomeLengths = {
        1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
        6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
        11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
        16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
        21: 46709983, 22: 50818468
    }

    df['Xposition'] = df.apply(lambda row: row['chr'] + row['pos'] / chromosomeLengths[row['chr']] - 0.5, axis=1)
    df['NegativeLogP'] = -np.log10(df['p'])

    maxY = df['NegativeLogP'].max() + 5

    dfCis = df[df['cisTrans'] == 'Cis']
    dfTrans = df[df['cisTrans'] == 'Trans']

    dfNovelCis = dfCis[dfCis['novel']]
    dfKnownCis = dfCis[~dfCis['novel']]
    dfNovelTrans = dfTrans[dfTrans['novel']]
    dfKnownTrans = dfTrans[~dfTrans['novel']]

    dfEvensCis = dfKnownCis[dfKnownCis['chr'] % 2 == 0]
    dfOddsCis = dfKnownCis[dfKnownCis['chr'] % 2 != 0]
    dfEvensTrans = dfKnownTrans[dfKnownTrans['chr'] % 2 == 0]
    dfOddsTrans = dfKnownTrans[dfKnownTrans['chr'] % 2 != 0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 20), dpi=200, gridspec_kw={'height_ratios': [1, 0.8]})

    ax1.scatter(dfNovelCis['Xposition'], dfNovelCis['NegativeLogP'], color=(240/256, 77/256, 66/256))
    ax1.scatter(dfEvensCis['Xposition'], dfEvensCis['NegativeLogP'], color='#1C1C8A')
    ax1.scatter(dfOddsCis['Xposition'], dfOddsCis['NegativeLogP'], color='#37379C')

    ax2.scatter(dfNovelTrans['Xposition'], dfNovelTrans['NegativeLogP'], color=(240/256, 77/256, 66/256))
    ax2.scatter(dfEvensTrans['Xposition'], dfEvensTrans['NegativeLogP'], color='#1C1C8A')
    ax2.scatter(dfOddsTrans['Xposition'], dfOddsTrans['NegativeLogP'], color='#37379C')

    for ax in [ax1, ax2]:
        ax.set_ylabel('-log10(p-value)', fontsize=20)
        ax.set_ylim(5, maxY)
        ax.set_xticks(np.arange(1, 23))
        ax.set_xticklabels(np.arange(1, 23), fontsize=20, y=-0.01)
        ax.tick_params(axis='x', bottom=False, top=False)
        ax.set_yticklabels(np.arange(0, 90, 10), fontsize=20)
        ax.margins(x=0.01)

    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    ax1.text(0.9, 0.9, 'Cis', fontsize=30, transform=ax1.transAxes, fontstyle='italic', ha='center', va='center')
    ax2.text(0.9, 0.1, 'Trans', fontsize=30, transform=ax2.transAxes, fontstyle='italic', ha='center', va='center')

    ax2.set_ylim(5, maxY * 0.8)
    ax2.invert_yaxis()
    ax2.tick_params(axis='x', labelbottom=False)

    ax1.axhline(y=-np.log10(5e-8), color='black', linestyle='--')
    ax2.axhline(y=-np.log10(6.58e-12), color='black', linestyle='--')

    fig.tight_layout()

    plt.savefig(os.path.join(outputDir, 'Miami.png'))
    plt.savefig(os.path.join(outputDir, 'Miami.svg'))

# --- Run Plotting ---
plotMiami(combinedDf)
