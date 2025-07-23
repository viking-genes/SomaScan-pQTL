import synapseclient
import synapseutils
import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Directory Paths ---
pietznerMeasuredProteins = pd.read_excel(os.path.join(
    projectRoot, 'Scripts/Publication/supportingFiles/replication/measuredProteins/Somalogic/Pietzner2020ProteinInfo.xlsx'))

supplementary1 = pd.read_excel(os.path.join(
    projectRoot, 'Scripts/Publication/supplementary1.xlsx'))

downloadLocation = os.path.join(
    projectRoot, 'Scripts/Publication/supportingFiles/replication/Pietzner2021SummaryStats')

syn = synapseclient.Synapse() 
syn.login(authToken="yourTokenHere") 

# synapse.org ID of the Pietzner2021 study
#https://www.synapse.org/Synapse:syn51824537
projectId = 'syn51824537'


# Function to find files by partial name within a project
def findFilesByPartialName(projectId, partialName):
    fileList = []
    # Use syn.getChildren to iterate over items in the project/folder
    for item in syn.getChildren(projectId):
        # Check if the item is a file and if the name contains the partial name
        if item['type'] == 'org.sagebionetworks.repo.model.FileEntity' and partialName in item['name']:
            # Fetch more details about the file
            fileEntity = syn.get(item['id'], downloadFile=False)
            fileInfo = {
                'id': fileEntity.id,
                'name': fileEntity.name,
                'content type': fileEntity.contentType,
                'modified on': fileEntity.modifiedOn
            }
            fileList.append(fileInfo)

    return fileList

# Function to download the non-MD5 file
def downloadNonMd5File(files, UniprotID):
    # Identify the file that is not an MD5 checksum file
    for file in files:
        if not file['name'].endswith('.md5'):
            # Download the identified file
            print(f"Downloading file: {file['name']} with ID: {file['id']} to {downloadLocation}")
            downloadedFile = syn.get(file['id'], downloadLocation=downloadLocation, downloadFile=True)

            newFileName = f"{UniprotID}_{file['name']}"
            newFilePath = os.path.join(downloadLocation, newFileName)

            # Rename the file
            os.rename(downloadedFile.path, newFilePath)

            print(UniprotID, "download complete.")
    return

for index in supplementary1.index:

    #Skip aptamers with multiple protein targets
    if '|' in supplementary1.loc[index, 'Uniprot']:
        continue

    #Check if the protein was measured in Pietzner2020
    if supplementary1.loc[index, 'Uniprot'] in pietznerMeasuredProteins['UniprotID'].values:

        #Find index of all the proteins in the Pietzner2020 measured proteins
        pietznerIdxList = pietznerMeasuredProteins[pietznerMeasuredProteins['UniprotID'] == supplementary1.loc[index, 'Uniprot']].index.tolist()

        for pietznerIdx in pietznerIdxList:
            #Retrieve all aptamer supplementary data for the protein
            partialName = pietznerMeasuredProteins.loc[pietznerIdx, 'SeqId']

            #Split by - join by _ to match the file name format
            partialName = partialName.split('-')
            partialName = '_'.join(partialName)
            print(partialName)
            #Check if the aptamer file is already downloaded by checking the truncated file name
            downloadedFiles = os.listdir(downloadLocation)
            if any(partialName in file for file in downloadedFiles):
                continue

            # Get the list of files that match the partial name
            matchingFiles = findFilesByPartialName(projectId, partialName)

            # Call the function with the list of files
            downloadNonMd5File(matchingFiles, supplementary1.loc[index, 'Uniprot'])
