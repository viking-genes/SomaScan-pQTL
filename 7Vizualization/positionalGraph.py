import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
supplementaryPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
outputDir = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/')
os.makedirs(outputDir, exist_ok=True)

chromosomeLengths = {
    1:248956422, 2:242193529, 3:198295559, 4:190214555, 5:181538259, 6:170805979,
    7:159345973, 8:145138636, 9:138394717, 10:133797422, 11:135086622, 12:133275309,
    13:114364328, 14:107043718, 15:101991189, 16:90338345, 17:83257441, 18:80373285,
    19:58617616, 20:64444167, 21:46709983, 22:50818468
}

# --- Color Map Helper ---
def createColorMap(colorList):
    def interFrom256(x):
        return np.interp(x=x, xp=[0, 255], fp=[0, 1])

    colorDict = {'red': [], 'green': [], 'blue': []}
    for channel in colorDict:
        for idx, rgb in enumerate(colorList):
            value = interFrom256(rgb[channel])
            position = idx / (len(colorList) - 1)
            colorDict[channel].append((position, value, value))

    return colors.LinearSegmentedColormap('effectGradient', segmentdata=colorDict)

# --- Plotting ---
def plotPositionalGraph(inputDf):
    # Drop ambiguous or non-autosomal entries
    df = inputDf.copy()
    df['X'] = np.nan
    df['Y'] = np.nan
    df['NegativeLogP'] = np.nan

    for index, row in df.iterrows():
        try:
            if 'X' in str(row['Somamer targeted protein chromosome']) or '|' in str(row['Somamer targeted protein TSS']):
                continue

            chrom = int(row['Chromosome'])
            tssChrom = int(row['Somamer targeted protein chromosome'])

            df.at[index, 'X'] = chrom + row['Position'] / chromosomeLengths[chrom] - 0.5
            df.at[index, 'Y'] = tssChrom + row['Somamer targeted protein TSS'] / chromosomeLengths[tssChrom] - 0.5
            df.at[index, 'NegativeLogP'] = -math.log10(row['P-value'])

        except Exception:
            df = df.drop(index)

    df = df.dropna().sort_values('NegativeLogP', ascending=False)

    fig, ax = plt.subplots(figsize=(12.00, 4.80), dpi=300)

    points = ax.scatter(
        df['X'],
        df['Y'],
        c=df['NegativeLogP'],
        s=[float(x) * 2 * math.sqrt(2) for x in df['NegativeLogP']],
        cmap=createColorMap([
            {'red': 204, 'green': 204, 'blue': 255},
            {'red': 28,  'green': 28,  'blue': 138}
        ])
    )

    ax.set_xlim(0, 23)
    ax.set_ylim(0, 23)
    ax.set_xticks(np.arange(1, 23, 1.0))
    ax.set_yticks(np.arange(1, 23, 1.0))
    ax.set_xticks(np.arange(0.5, 23, 1.0), minor=True)
    ax.set_yticks(np.arange(0.5, 23, 1.0), minor=True)
    ax.grid(which='minor', color='#E6E6E6')
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_axisbelow(True)

    # Outline style
    for spine in ax.spines.values():
        spine.set_color('#E5E5E5')

    ax.set_xlabel('pQTL Position')
    ax.set_ylabel('Gene TSS Position')

    colorbar = plt.colorbar(points)
    colorbar.outline.set_visible(False)
    colorbar.ax.set_ylabel('-log(P-value)', rotation=270, labelpad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(outputDir, 'positionalGraph.svg'), transparent=True)
    plt.savefig(os.path.join(outputDir, 'positionalGraph.png'), transparent=True)

    return df

# --- Run ---
supplementaryDf = pd.read_excel(supplementaryPath)
cleanedDf = plotPositionalGraph(supplementaryDf)
