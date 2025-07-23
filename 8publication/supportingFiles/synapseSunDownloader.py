import synapseclient
import synapseutils
import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Download Location ---
downloadLocation = os.path.join(
    projectRoot,
    'Scripts/Publication/supportingFiles/replication/Sun2023SummaryStats'
)

syn = synapseclient.Synapse() 
syn.login(authToken="yourTokenHere") 

# synapse.org ID of the Sun2023 study
#https://www.synapse.org/Synapse:syn51364943
projectId = 'syn51365303'


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


# proteinsToDownload = ['Q9H7C9', 'Q4VC05', 'Q9HB40', 'Q9UNN8']
proteinsToDownload = ['P02766'] # This is for BridgeBio at ASHG2024

for protein in proteinsToDownload:

    # Get the list of files that match the partial name
    matchingFiles = findFilesByPartialName(projectId, protein)

    # Call the function with the list of files
    downloadNonMd5File(matchingFiles, protein)
