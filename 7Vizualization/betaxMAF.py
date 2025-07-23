import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import os
from matplotlib import colors

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
supplementaryPath = os.path.join(projectRoot, 'Scripts/Publication/supplementary1.xlsx')
outputDir = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/')
os.makedirs(outputDir, exist_ok=True)

# --- Load data ---
supplementaryDf = pd.read_excel(supplementaryPath)

# --- Utility Functions ---
def filterMinorAlleleFrequency(df):
    for i in range(len(df)):
        maf = abs(df.at[i, 'Effect allele frequency in VIKING'])
        if maf > 0.5:
            df.at[i, 'Effect allele frequency in VIKING'] = 1 - maf
    return df

def createColorMap(colorList):
    def interFrom256(x):
        return np.interp(x=x, xp=[0, 255], fp=[0, 1])
    
    colorDict = {'red': [], 'green': [], 'blue': []}
    for channel in colorDict:
        for idx, color in enumerate(colorList):
            value = interFrom256(color['rgb'][channel])
            pos = idx / (len(colorList) - 1)
            colorDict[channel].append((pos, value, value))
    
    return colors.LinearSegmentedColormap('customColormap', segmentdata=colorDict)

# --- Plot Functions ---
def plotScatter(df):
    cisDf = df[df['cisTrans'] == 'Cis']
    transDf = df[df['cisTrans'] == 'Trans']

    plt.scatter(cisDf['Effect allele frequency in VIKING'], cisDf['Effect Size (beta)'].abs(), color='blue')
    plt.scatter(transDf['Effect allele frequency in VIKING'], transDf['Effect Size (beta)'].abs(), color='red')

    plt.legend(['Cis p < 5e-8', 'Trans p < 6.6e-12'])
    plt.xlabel('Minor Allele Frequency (MAF)', fontsize=18)
    plt.ylabel('Absolute effect size', fontsize=18)
    plt.savefig(os.path.join(outputDir, 'betaxMAF.png'))

def plotScatterWithColorbar(df):
    df['NegativeLogP'] = np.nan

    for index, row in df.iterrows():
        # Skip ambiguous or X-chromosome entries
        if 'X' in str(row['Somamer targeted protein chromosome']) or '|' in str(row['Somamer targeted protein TSS']):
            df.drop(index, inplace=True)
        else:
            df.at[index, 'NegativeLogP'] = -math.log10(row['P-value'])

    df = df.sort_values('NegativeLogP', ascending=False)

    cisDf = df[df['cisTrans'] == 'Cis']
    transDf = df[df['cisTrans'] == 'Trans']

    fig, ax = plt.subplots(figsize=(12.00, 4.80), dpi=1000)

    # Define colormaps
    redMap = createColorMap([
        {'rgb': {'red': 251, 'green': 210, 'blue': 208}},
        {'rgb': {'red': 255, 'green': 40,  'blue': 26}}
    ])
    blueMap = createColorMap([
        {'rgb': {'red': 204, 'green': 204, 'blue': 255}},
        {'rgb': {'red': 28,  'green': 28,  'blue': 138}}
    ])

    # Plot
    pointsCis = ax.scatter(
        cisDf['Effect allele frequency in VIKING'],
        cisDf['Effect Size (beta)'].abs(),
        c=cisDf['NegativeLogP'], cmap=redMap
    )
    pointsTrans = ax.scatter(
        transDf['Effect allele frequency in VIKING'],
        transDf['Effect Size (beta)'].abs(),
        c=transDf['NegativeLogP'], cmap=blueMap
    )

    ax.set_xlim((0, 0.52))
    ax.set_xlabel('Minor allele frequency (MAF)')
    ax.set_ylabel('Effect size (beta)')
    ax.xaxis.labelpad = 8
    ax.yaxis.labelpad = 12

    # Add colorbars
    cisBar = plt.colorbar(pointsCis, pad=0.02)
    transBar = plt.colorbar(pointsTrans, pad=0, anchor=(-0.5, 0.167), shrink=0.684)

    cisBar.ax.set_xlabel('Cis')
    transBar.ax.set_xlabel('Trans', labelpad=20)
    cisBar.outline.set_visible(False)
    transBar.outline.set_visible(False)
    fig.axes[1].set_ylabel('-log(P-value)', rotation=270, labelpad=20)

    # Plot styling
    for spine in ax.spines.values():
        spine.set_color('#E6E6E6')
        spine.set_linewidth(0.8)

    plt.savefig(os.path.join(outputDir, 'betaxMAF.svg'))

# --- Run ---
supplementaryDf = filterMinorAlleleFrequency(supplementaryDf)
plotScatter(supplementaryDf)
plotScatterWithColorbar(supplementaryDf)
