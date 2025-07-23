import openpyxl
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Output Folder ---
outputFolder = os.path.join(projectRoot, 'Scripts/Publication/outputFiles')

# --- Define Supplementary Table Paths ---
supplementaryTables = {
    os.path.join(projectRoot, 'Scripts/Publication/supplementary1Replication.xlsx'): "supplementary1.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supplementary1ReplicationSurapaneni.xlsx'): "supplementary2AASKCisReplication.xlsx",
    os.path.join(projectRoot, 'Scripts/SNPnovelty3/somalogic7kAnnotation.xlsx'): "supplementary3Somalogic7kTargets.xlsx",
    os.path.join(projectRoot, 'Scripts/SNPnovelty3/somalogic4Olink1500Combined.xlsx'): "supplementary4PreviouslyPublished.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supplementaryMRReplication.xlsx'): "supplementary5MR.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supplementaryColoc.xlsx'): "supplementary6Coloc.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/5replicationDf.xlsx'): "supplementary7MRreplication.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocResultDf.xlsx'): "supplementary8ColocalisationWithReplicationDatasets.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supplementaryGOEnrichment.xlsx'): "supplementary9GOEnrichment.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supplementaryCohortSummary.xlsx'): "supplementary10CohortSummary.xlsx",
    os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/ChatGPTCategorisation.xlsx'): "supplementary11ChatGPTCategorisation.xlsx",
    os.path.join(projectRoot, 'Scripts/twoSampleMr/createOpenGWASDatabase/openGWASChatGPTAnnotated.xlsx'): "supplementary12openGWASChatGPTAnnotated.xlsx",
    os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocReplicationResultDf.xlsx'): "supplementary13colocReplication.xlsx",
}

# --- Create Output Folder If It Doesn't Exist ---
os.makedirs(outputFolder, exist_ok=True)

# --- Freeze First Row and Save ---
for file, output in supplementaryTables.items():
    wb = openpyxl.load_workbook(file)
    ws = wb.active
    ws.freeze_panes = ws['A2']
    outputPath = os.path.join(outputFolder, output)
    wb.save(outputPath)
    wb.close()
    print(f"Saved {output}")
