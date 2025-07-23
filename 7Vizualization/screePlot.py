import os
import pandas as pd
import matplotlib.pyplot as plt

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
eigenvaluePath = os.path.join(projectRoot, 'GWAS/GRM/pca20.eigenval')
outputDir = os.path.join(projectRoot, 'Scripts/Visualization/Graphs/')
os.makedirs(outputDir, exist_ok=True)

# --- Load eigenvalues and prepare for scree plot ---
eigenvalueDf = pd.read_csv(eigenvaluePath, sep='\t', header=None)
eigenvalueDf.columns = ['eigenvalue']

# Keep first 10 components (0-indexed to 9)
eigenvalueDf = eigenvalueDf.iloc[:10]
eigenvalueDf.index = eigenvalueDf.index + 1  # PC1 to PC10

# --- Plot scree plot ---
plt.figure(figsize=(6, 4), dpi=300)
plt.plot(
    eigenvalueDf.index,
    eigenvalueDf['eigenvalue'],
    linewidth=2,
    color='black',
    marker='o',
    markerfacecolor='black',
    markersize=8
)

plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Eigenvalue')
plt.xticks(eigenvalueDf.index)

# Annotate inflection point (PC3)
plt.annotate(
    'Inflection Point',
    xy=(3, eigenvalueDf.loc[3, 'eigenvalue']),
    xytext=(3.5, eigenvalueDf.loc[3, 'eigenvalue'] + 1.5),
    arrowprops=dict(facecolor='black', shrink=0.05)
)

# --- Save figure ---
outputPath = os.path.join(outputDir, 'screePlot.svg')
plt.tight_layout()
plt.savefig(outputPath)
