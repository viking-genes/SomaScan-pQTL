import pandas as pd
import os

# --- Setup ---
projectRoot = './prj_190_viking_somalogic/'

# Define input paths
publishedProteinsPath = os.path.join(projectRoot, 'Scripts/SNPnovelty2/1previouslyPublished.xlsx')
somalogicAnnotationPath = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/SomalogicTargets/extractSomalogicUniprotHugoAnnotation.xlsx')

# Define output path
outputPath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryNovelty.xlsx')

# --- Load data ---
df = pd.read_excel(publishedProteinsPath)
df7kSomalogic = pd.read_excel(somalogicAnnotationPath)

# Remove empty rows in HUGO column
df = df[df['HUGO'].notna()]

# --- Write to Excel with two sheets ---
with pd.ExcelWriter(outputPath) as writer:
    df.to_excel(writer, sheet_name='previouslyPublishedProteins', index=False)
    df7kSomalogic.to_excel(writer, sheet_name='somaLogicv4.1AssayTargets', index=False)
