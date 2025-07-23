import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ncx2, chi2
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'
powerOutputDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/powerCalculation/')
os.makedirs(powerOutputDir, exist_ok=True)
outputPlotPath = os.path.join(powerOutputDir, 'power_calculation.png')

# --- Parameters ---
gwasSignificanceThreshold = 5e-8
sampleSizes = np.arange(10, 401, 1)  # Sample sizes from 10 to 400
varianceExplainedList = [0.1, 0.2, 0.5]  # Heritability H² values


# --- Power Calculation Function ---
def calculateExactPower(sampleSize, varianceExplained, significanceThreshold=gwasSignificanceThreshold):
    """
    Calculate exact statistical power using non-central chi-squared distribution.

    Args:
        sampleSize (int): Number of individuals in the study.
        varianceExplained (float): Proportion of variance explained (H²).
        significanceThreshold (float): Genome-wide significance threshold (e.g., 5e-8).

    Returns:
        float: Power (1 - β) of detecting association.
    """
    significanceThresholdChi = chi2.ppf(1 - significanceThreshold, df=1)
    nonCentralityParameter = sampleSize * varianceExplained
    return ncx2.sf(significanceThresholdChi, df=1, nc=nonCentralityParameter)


# --- Power Curve Generation ---
def generatePowerCurves():
    """
    Plot and export power curves for different levels of variance explained (H²).
    """
    plt.figure(figsize=(10, 6))

    for h2 in varianceExplainedList:
        powerValues = [
            calculateExactPower(sampleSize=n, varianceExplained=h2)
            for n in sampleSizes
        ]
        plt.plot(sampleSizes, powerValues, label=f'$H^2 = {h2}$', linewidth=2)

    # Plot annotations and styling
    plt.title('Statistical Power vs. Sample Size\nfor Hypothetical SNP Variance Explained (H²)', fontsize=14, pad=20)
    plt.xlabel('Sample Size ($N$)', fontsize=12)
    plt.ylabel('Study Power (1 - $\\beta$)', fontsize=12)
    plt.axhline(y=0.8, color='gray', linestyle='--', label='80% Power (Reference)')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1.05)
    plt.legend(fontsize=12, framealpha=1)
    plt.tight_layout()

    # Save output
    plt.savefig(outputPlotPath, dpi=300, bbox_inches='tight')
    plt.close()


# --- Run ---
generatePowerCurves()
