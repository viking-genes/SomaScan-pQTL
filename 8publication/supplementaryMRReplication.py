import pandas as pd
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Supplementary Files ---
MRdf = os.path.join(projectRoot, 'Scripts/Publication/supplementaryMR.xlsx')
colocDf = os.path.join(projectRoot, 'Scripts/Publication/supplementaryColoc.xlsx')
MRReplicationdf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/MR/5replicationDf.xlsx')
colocReplicationdf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocReplicationResultDf.xlsx')
VIKINGcolocdf = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/colocalisation/colocResultDf.tsv')


MRdf = pd.read_excel(MRdf)
MRdf = MRdf.fillna("NA")
colocDf = pd.read_excel(colocDf)
MRReplicationdf = pd.read_excel(MRReplicationdf)
colocReplicationdf = pd.read_excel(colocReplicationdf)
VIKINGcolocdf = pd.read_csv(VIKINGcolocdf, sep='\t')

#Add colocalisation results to MRdf with a Yes/No column
def addColocalisation(df, colocDf):

    df['Exposure and Outcome colocalises'] = ""

    for i in range(len(df)):

        exposure = df.loc[i, 'MR Exposure']
        outcome = df.loc[i, 'Study ID']

        #Find the row in the colocalisation df that matches the exposure and outcome
        colocalisationRow = colocDf[(colocDf['HUGO'] == exposure) & (colocDf['Open GWAS Dataset ID'] == outcome)]

        #Check if it colocalises
        if colocalisationRow['PP.H4.abf'].values[0] > 0.8:
            df.loc[i, 'Exposure and Outcome colocalises'] = "Yes"
        else:
            df.loc[i, 'Exposure and Outcome colocalises'] = "No"

    return df

# Add MR replication results to MRdf with a Yes/No column and additional data if Yes
def addMRReplication(df, replicationDf):
    
    df['MR replicates'] = ""
    df['MR replication significant'] = ""
    df['MR replication matching directionality'] = ""
    df['MR replication type'] = ""
    df['MR replication instrument'] = ""
    df['MR replication comment'] = ""

    for i in range(len(df)):
    
        #Skip if it did not colocalise in the discovery analysis
        if df.loc[i, 'Exposure and Outcome colocalises'] == "No":
            df.loc[i, 'MR replicates'] = "NA"
            df.loc[i, 'MR replication significant'] = "NA"
            df.loc[i, 'MR replication matching directionality'] = "NA"
            df.loc[i, 'MR replication type'] = "NA"
            df.loc[i, 'MR replication instrument'] = "NA"
            continue

        exposure = df.loc[i, 'MR Exposure']
        outcome = df.loc[i, 'Study ID']

        #No valid instrument for LTK replication
        if exposure == 'LTK':
            df.loc[i, 'MR replicates'] = "NA"
            df.loc[i, 'MR replication significant'] = "NA"
            df.loc[i, 'MR replication matching directionality'] = "NA"
            df.loc[i, 'MR replication type'] = "NA"
            df.loc[i, 'MR replication instrument'] = "NA"
            df.loc[i, 'MR replication comment'] = "No valid instrument for LTK replication"
            continue
        
        #No valid instrument for CYP2C19 replication
        if exposure == 'CYP2C19':
            df.loc[i, 'MR replicates'] = "NA"
            df.loc[i, 'MR replication significant'] = "NA"
            df.loc[i, 'MR replication matching directionality'] = "NA"
            df.loc[i, 'MR replication type'] = "NA"
            df.loc[i, 'MR replication instrument'] = "NA"
            df.loc[i, 'MR replication comment'] = "No valid instrument for CYP2C19 replication"
            continue

        #Find the row in the replication df that matches the exposure and outcome
        replicationRows = replicationDf[(replicationDf['Replication MR Exposure'].str.contains(exposure)) & (replicationDf['Discovery MR Study ID'].str.contains(outcome))]

        MRreplicates = []
        MRreplicationSignificant = []
        MRreplicationType = []
        MRreplicationInstrument = []
        MRreplicationffectSizeDirectionality = []

        #Iterate over MR instruments
        for j in replicationRows.index:
            MRreplicationSignificant.append(replicationRows.loc[j, 'MR Significant result'])

            #If it does replicate, identify the replication type and check if the effect size directionality is the same between discovery and replication
            if replicationRows.loc[j, 'MR Significant result'] == "Yes":
                
                #Take the average effect size in case of multiple instruments
                discoveryEffectSize = sum([float(x) for x in str(replicationRows.loc[j, 'Discovery MR Effect Size (beta)']).split("|")]) / len(str(replicationRows.loc[j, 'Discovery MR Effect Size (beta)']).split("|"))
                
                #Check if Stage 2 is NA, meaning it was significant in Stage 1
                if replicationRows.loc[j, 'Replication MR Stage 2 p-value'] != replicationRows.loc[j, 'Replication MR Stage 2 p-value']:
                    MRreplicationType.append("Independent")

                    replicationEffectSize = sum([float(x) for x in str(replicationRows.loc[j, 'Replication MR Stage 1 Effect size (beta)']).split("|")]) / len(str(replicationRows.loc[j, 'Replication MR Stage 1 Effect size (beta)']).split("|"))

                    #Check if the effect size directionality is the same
                    if replicationEffectSize * discoveryEffectSize > 0:
                        MRreplicationffectSizeDirectionality.append("Yes")
                    else:
                        MRreplicationffectSizeDirectionality.append("No")

                else:
                    MRreplicationType.append("Discovery MR Outcome")
                    
                    replicationEffectSize = sum([float(x) for x in str(replicationRows.loc[j, 'Replication MR Stage 2 effect size (beta)']).split("|")]) / len(str(replicationRows.loc[j, 'Replication MR Stage 1 Effect size (beta)']).split("|"))

                    if replicationEffectSize * discoveryEffectSize > 0:
                        MRreplicationffectSizeDirectionality.append("Yes")
                    else:
                        MRreplicationffectSizeDirectionality.append("No")

            else:
                MRreplicationType.append("NA")
                MRreplicationffectSizeDirectionality.append("NA")


            #Add the replication instrument
            MRinstrument = replicationRows.loc[j, 'Replication MR Exposure'].split("_")[1]
            MRreplicationInstrument.append(MRinstrument)

        #Create the MR replicates column
        for j in range(len(MRreplicationSignificant)):
            print(MRreplicationSignificant, MRreplicationffectSizeDirectionality)
            if MRreplicationSignificant[j] == "Yes" and MRreplicationffectSizeDirectionality[j] == "Yes":
                MRreplicates.append("Yes")
            else:
                MRreplicates.append("No")

        #Add the data to the MRdf
        df.loc[i, 'MR replicates'] = '|'.join(MRreplicates)
        df.loc[i, 'MR replication significant'] = '|'.join(MRreplicationSignificant)
        df.loc[i, 'MR replication matching directionality'] = '|'.join(MRreplicationffectSizeDirectionality)
        df.loc[i, 'MR replication type'] = '|'.join(MRreplicationType)
        df.loc[i, 'MR replication instrument'] = '|'.join(MRreplicationInstrument)

    return df

#Prioritises the display of MR replication results
def filterMRreplication(df):
    
    for i in df.index:

        if df.loc[i, 'MR replication significant'] == "NA":
            continue

        MRReplicatesList = df.loc[i, 'MR replication significant'].split('|')
        MRReplicationTypeList = df.loc[i, 'MR replication type'].split('|')
        MRReplicationInstrumentList = df.loc[i, 'MR replication instrument'].split('|')

        #Remove non-replicating MR instruments
        for j in range(len(MRReplicatesList)):
            if MRReplicatesList[j] == "No":
                MRReplicatesList[j] = ""
                MRReplicationTypeList[j] = ""
                #Only remove the instrument from the list if there is another replicating instrument
                if "Yes" in MRReplicatesList:
                    MRReplicationInstrumentList[j] = ""
        
        MRReplicatesList = list(filter(None, MRReplicatesList))
        MRReplicationTypeList = list(filter(None, MRReplicationTypeList))
        MRReplicationInstrumentList = list(filter(None, MRReplicationInstrumentList))

        #If it has been all Non-replicating, set to No
        if len(MRReplicatesList) == 0:
            MRReplicatesList = ['No']
            MRReplicationTypeList = ['NA']
            # MRReplicationInstrumentList = ['NA']

        if 'Yes' in MRReplicatesList:
            MRReplicatesList = ['Yes']

        if 'Independent' in MRReplicationTypeList:
            MRReplicationTypeList = ['Independent']

            #Find index of the independent replication
            index = MRReplicationTypeList.index('Independent')
            MRReplicationInstrumentList = [MRReplicationInstrumentList[index]]
        
        elif 'Discovery MR Outcome' in MRReplicationTypeList:
            MRReplicationTypeList = ['Discovery MR Outcome']

        #Save the filtered lists back to the df
        df.at[i, 'MR replication significant'] = '|'.join(MRReplicatesList)
        df.at[i, 'MR replication type'] = '|'.join(MRReplicationTypeList)
        df.at[i, 'MR replication instrument'] = '|'.join(MRReplicationInstrumentList)

    return df

#Adds a column to the MRdf that indicates if the MR instrument colocalises with VIKING
def addInstrumentColocVIKING(df, VIKINGcolocdf):

    df["Instrument colocalises with VIKING"] = ""

    for i in range(len(df)):

        #Skip if it there was no MR replication
        if df.loc[i, 'MR replication instrument'] != df.loc[i, 'MR replication instrument'] or df.loc[i, 'MR replication instrument'] == "NA":
            df.at[i, 'Instrument colocalises with VIKING'] = "NA"
            continue

        exposure = df.loc[i, 'MR Exposure']

        instrumentList = df.loc[i, 'MR replication instrument'].split('|')
        result = []
        #Check if the instrument colocalises with VIKING
        for instrument in instrumentList:

            #Rename the assay to match across the two dfs
            if instrument == "Somalogic":
                instrument = "somalogic"
            elif instrument == "Olink":
                instrument = "olink"

            colocVIKINGrow = VIKINGcolocdf[(VIKINGcolocdf['Gene name'] == exposure) & (VIKINGcolocdf['comparison assay'] == instrument)]
            
            if colocVIKINGrow['H4'].values[0] > 0.8:
                result.append("Yes")
            else:
                result.append("No")

        df.at[i, 'Instrument colocalises with VIKING'] = '|'.join(result)

    #Place the column in the second last position
    col = df.pop("Instrument colocalises with VIKING")
    df.insert(len(df.columns) - 1, "Instrument colocalises with VIKING", col)

    return df



MRdf = addColocalisation(MRdf, colocDf)
MRdf = addMRReplication(MRdf, MRReplicationdf)
MRdf = filterMRreplication(MRdf)
MRdf = addInstrumentColocVIKING(MRdf, VIKINGcolocdf)

#Sort by Exposure and Outcome colocalises with Yes on top
MRdf = MRdf.sort_values(by=['Exposure and Outcome colocalises'], ascending=False)

# --- Export Results ---
outputPath = os.path.join(projectRoot, 'Scripts/Publication/supplementaryMRReplication.xlsx')
MRdf.to_excel(outputPath, index=False)

