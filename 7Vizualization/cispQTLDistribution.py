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

# --- Load Data ---
supplementaryDf = pd.read_excel(supplementaryPath)

# --- Helper Functions ---
def extractCisHits(dataframe):
    return dataframe[dataframe['cisTrans'] == 'Cis'].reset_index(drop=True)

def addNegativeLogPValues(dataframe):
    dataframe['-log(p)'] = [-math.log(pval, 10) for pval in dataframe['P-value']]
    return dataframe

def calculateTSSDistance(dataframe):
    distances = []
    for i in range(len(dataframe)):
        tss = dataframe.loc[i, 'Somamer targeted protein TSS']
        position = dataframe.loc[i, 'Position']
        distances.append(position - tss if isinstance(tss, int) else np.nan)
    dataframe['TSS distance'] = distances
    return dataframe

def computeAbsoluteEffectSizes(dataframe):
    dataframe['Effect size'] = dataframe['Effect Size (beta)'].abs()
    return dataframe

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

# --- Plotting Function ---
def plotTSSDistanceGraph(dataframe):
    dataframe = dataframe.sort_values(by='Effect size', ascending=True)
    fig, ax = plt.subplots(figsize=(12.00, 4.80), dpi=1000)

    colormap = createColorMap([
        {'red': 204, 'green': 204, 'blue': 255},
        {'red': 28,  'green': 28,  'blue': 138}
    ])

    points = ax.scatter(
        dataframe['TSS distance'],
        dataframe['-log(p)'],
        c=dataframe['Effect size'],
        cmap=colormap
    )

    ax.set_xlim(-1_000_000, 1_000_000)
    ax.set_xlabel('Distance from TSS, Mb')
    ax.set_ylabel('-log(P-value)')
    ax.grid(color='#E6E6E6')

    colorbar = plt.colorbar(points)
    colorbar.outline.set_visible(False)
    colorbar.ax.set_ylabel('Absolute Effect Size', rotation=270, labelpad=20)

    plt.tight_layout()
    plt.savefig(os.path.join(outputDir, 'TSSdistance.svg'))

# --- Run Pipeline ---
cisHitsDf = extractCisHits(supplementaryDf)
cisHitsDf = addNegativeLogPValues(cisHitsDf)
cisHitsDf = calculateTSSDistance(cisHitsDf)
cisHitsDf = computeAbsoluteEffectSizes(cisHitsDf)
plotTSSDistanceGraph(cisHitsDf)
