import pandas as pd
import os
import subprocess

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define File and Directory Paths ---
df = pd.read_excel(os.path.join(projectRoot, 'Scripts/Publication/supplementary1Replication.xlsx'))
outputDf = os.path.join(projectRoot, 'Scripts/Publication/supplementary1ReplicationSurapaneni.xlsx')

SurapaneniDf = pd.read_excel(os.path.join(
    projectRoot, 'Scripts/Publication/supportingFiles/replication/Surapaneni2022Data/Surapaneni2022.xlsx'))

LDdir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/LDProxiesUKBB/ProcessedOutput')


VIKINGdf = df.copy()
#Cohort name is AASK in the Surapaneni data (African American Study of Kidney Disease and Hypertension)
df['Replication AASK rsid'] = ''
df['Replication AASK LD'] = ''
df['Replication AASK p-value'] = ''
df['Replication AASK Effect Size'] = ''
df['Replication AASK Allele Frequency'] = ''
df['Replication AASK SomamerID'] = ''

def checkIfProxy(rsidVIKING, rsidSurapaneni, chromosome, R2 = 0.6):

    #Check if the rsid is a proxy for the rsid in UKBB 10k dataset with R2 >= 0.4
    LDfile = os.path.join(LDdir, f'ProxyList_chr{chromosome}.tsv')

    LDfileDf = pd.read_csv(LDfile, sep = '\t')

    #Extract all proxy dataframe indexes for the rsid in VIKING
    proxiesVIKING = LDfileDf[LDfileDf['SNP_A'] == rsidVIKING].index.tolist()

    if len(proxiesVIKING) == 0:
        print('Error: No proxies found for rsid', rsidVIKING, 'in UKBB 10k dataset')
        exit()

    #Leave only the proxies with the specified R2. Currently 0.6 is the minimum R2
    proxiesVIKING = LDfileDf.loc[proxiesVIKING]
    proxiesVIKING = proxiesVIKING[proxiesVIKING['R2'] >= R2]

    if rsidSurapaneni in proxiesVIKING['SNP_B'].tolist():
        #Find index of the proxy in the dataframe
        proxyIndex = proxiesVIKING[proxiesVIKING['SNP_B'] == rsidSurapaneni].index.tolist()[0]
        return proxiesVIKING.loc[proxyIndex, 'R2']

    else:
        return 'NA'


#Checks if the alleles are flipped between the two datasets and flips the beta and allele frequency if necessary
def checkAlleleFlip(rsidVIKING, effectAlleleVIKING, otherAlleleVIKING, rsidSurapaneni, codedAlleleSurapaneni, betaSurapaneni, alleleFrequencySurapaneni, chromosome):

    #PLINK command to extract the phasing of two SNPs
    def calculatePhasing(snp1, snp2, chr):
        # Define the path to the PLINK executable
        plinkPath = "./1.90p/plink"
        
       # --- Define the base file path ---
        bfilePath = os.path.join(
            baseDataRoot,
            '/genotypes/10k_unrelated_white_british_reference',
            f'ukbb_chr{chr}_10000_random_unrelated_white_british'
        )

        # --- Define the output file path ---
        outputPath = os.path.join(
            projectRoot,
            'Scripts/Publication/supportingFiles/replication/allelePhasingOutput'
        )
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
        outputPath += f"/{snp1}_{snp2}"
        
        #Check if the output file already exists
        if os.path.exists(outputPath + '.log'):
            return

        # Construct the PLINK command
        command = [
            plinkPath,
            "--bfile", bfilePath,
            "--ld", snp1, snp2,
            "--out", outputPath
        ]
        
        # Run the command using subprocess
        try:
            with open('/dev/null', 'w') as devnull:
                subprocess.run(command, stdout=devnull, stderr=devnull, check=True)
            print(f"{snp1, snp2} phasing completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running PLINK: {e}")
            exit()

    #Returns the in phase alleles in the form of {'EffectAlleleVIKING': 'X', 'InphaseAllele': 'Y', 'OtherAlleleVIKING': 'A', 'OutOfPhaseAllele': 'B'}
    def extractInPhaseAlleles(snp1, snp2):

        # --- Define the Path to the PLINK Output File ---
        filepath = os.path.join(
            projectRoot,
            'Scripts/Publication/supportingFiles/replication/allelePhasingOutput',
            f'{snp1}_{snp2}.log'
        )
        result = {'EffectAlleleVIKING': '', 'InphaseAllele': '', 'OtherAlleleVIKING': '', 'OutOfPhaseAllele': ''}

        # Open the file for reading
        with open(filepath, 'r') as file:

            haplotypeLineStart = len(file.readlines())

            #Extract the haplotypes from the LD file and return the in phase alleles
            file.seek(0)
            for num, line in enumerate(file):

                # Check if the line contains 'In phase alleles'
                if 'In phase alleles' in line:
                    #Return the line that contains the 'In phase alleles'
                    resultAlleles = line.strip().split(' ')[-1]
                    resultAlleles = resultAlleles.split('/')
                    
                    #Split the alleles if they are single nucleotide
                    if len(resultAlleles[0]) == 2 and len(resultAlleles[1]) == 2:
                        result['EffectAlleleVIKING'] = resultAlleles[0][0]
                        result['InphaseAllele'] = resultAlleles[0][1]
                        result['OtherAlleleVIKING'] = resultAlleles[1][0]
                        result['OutOfPhaseAllele'] = resultAlleles[1][1]
                        return result
                        
                    else:
                        print('In phase alleles are not single nucleotide. Removing', resultAlleles, snp1, snp2)
                        return None

        print('Error: No in phase alleles found in file', filepath)
        exit()

    #Check if the alleles in VIKING are flipped compared to UKBB
    def matchVIKINGPhasing(allelePhaseDict, effectVIKING, otherVIKING, rsid):

        if effectVIKING == allelePhaseDict['EffectAlleleVIKING'] and otherVIKING == allelePhaseDict['OtherAlleleVIKING']:
            return allelePhaseDict

        if effectVIKING == allelePhaseDict['OtherAlleleVIKING'] and otherVIKING == allelePhaseDict['EffectAlleleVIKING']:
            return {'EffectAlleleVIKING': allelePhaseDict['OtherAlleleVIKING'], 'InphaseAllele': allelePhaseDict['OutOfPhaseAllele'], 'OtherAlleleVIKING': allelePhaseDict['EffectAlleleVIKING'], 'OutOfPhaseAllele': allelePhaseDict['InphaseAllele']}
            return allelePhaseDict

        print('Error: Alleles do not match between VIKING and UKBB for', rsid, effectVIKING, otherVIKING, allelePhaseDict)
        exit()

    #Coded allele in Surapaneni seems to match the effect allele in VIKING
    #Compares the alleles between the two datasets and flips the effect size & frequency if necessary
    def flipAlleles(allelePhaseDict, codedAlleleSurapaneni, effectSize, alleleFrequency):
        # {'EffectAlleleVIKING': 'X', 'InphaseAllele': 'Y', 'OtherAlleleVIKING': 'A', 'OutOfPhaseAllele': 'B'}

        if codedAlleleSurapaneni == allelePhaseDict['InphaseAllele']:
            return effectSize, alleleFrequency

        if codedAlleleSurapaneni == allelePhaseDict['OutOfPhaseAllele']:
            return -1 * float(effectSize), 1 - float(alleleFrequency)

        else:
            print('Error: Allele not found in allele phase. Likely strand issues, these are not accounted for, dropping.', allelePhaseDict, codedAlleleSurapaneni)
            exit()

    #If the reported rsid in the two studies are the same, only check if it is flipped between the studies
    if rsidVIKING == rsidSurapaneni:
        
        if effectAlleleVIKING == codedAlleleSurapaneni:
            return betaSurapaneni, alleleFrequencySurapaneni

        if otherAlleleVIKING == codedAlleleSurapaneni:
            return -1 * float(betaSurapaneni), 1 - float(alleleFrequencySurapaneni)

        else:
            print('Error! Reported allele does not match neither effect nor other alleles in VIKING', rsidVIKING, effectAlleleVIKING, otherAlleleVIKING, codedAlleleSurapaneni)
            exit()

    #Check in-phase alleles between the two datasets
    calculatePhasing(rsidVIKING, rsidSurapaneni, chromosome)

    #Parse the PLINK phasing log file
    phaseAlleles = extractInPhaseAlleles(rsidVIKING, rsidSurapaneni)

    #Some alleles cannot be processed this way, so they are dropped. When they report CA as coded allele but there is no match in the UKBB data for the rsid rs3214535
    if phaseAlleles == None:
        return 'NA', 'NA'

    #Check if allele annotation from UKBB matches VIKING effect/other allele annotation. If not, flip the alleles to match
    phaseAlleles = matchVIKINGPhasing(phaseAlleles, effectAlleleVIKING, otherAlleleVIKING, rsidVIKING)

    #Flip the Surapaneni effect size and allele frequency if necessary
    betaSurapaneni, alleleFrequencySurapaneni = flipAlleles(phaseAlleles, codedAlleleSurapaneni, betaSurapaneni, alleleFrequencySurapaneni)

    return betaSurapaneni, alleleFrequencySurapaneni

def checkWhatProxy(rsidVIKING, rsidSurapaneni, chromosome):
    #PLINK command to extract the phasing of two SNPs
    def runPLINKphasing(snp1, snp2, chr, outputPath):
        # Define the path to the PLINK executable
        plinkPath = "/1.90p/plink"
        
        # --- Define the Base File Path ---
        bfilePath = os.path.join(
            projectRoot,
            '/genotypes/10k_unrelated_white_british_reference',
            f'ukbb_chr{chr}_10000_random_unrelated_white_british'
        )
        
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
        outputPath += f"/{snp1}_{snp2}"
        
        #Check if the output file already exists
        if os.path.exists(outputPath + '.log'):
            return

        # Construct the PLINK command
        command = [
            plinkPath,
            "--bfile", bfilePath,
            "--ld", snp1, snp2,
            "--out", outputPath
        ]
        
        # Run the command using subprocess
        try:
            with open('/dev/null', 'w') as devnull:
                subprocess.run(command, stdout=devnull, stderr=devnull, check=True)
            print(f"{snp1, snp2} phasing completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running PLINK: {e}")
            exit()

    # --- Define the Output File Path ---
    outputPath = os.path.join(
        projectRoot,
        'Scripts/Publication/supportingFiles/replication/allelePhasingOutput'
    )

    #Check if the file exists
    filePath = os.path.join(outputPath, f'{rsidVIKING}_{rsidSurapaneni}.log')
    if not os.path.exists(filePath):
        runPLINKphasing(rsidVIKING, rsidSurapaneni, chromosome, outputPath)

    #Read the file and extract the R2 value
    with open(filePath, 'r') as file:
        for line in file:
            if 'R-sq' in line:
                LD = float(line.strip().split(' ')[2])
                break

    return LD

#Create a dummy file with the R2 value of 1
def createDummyLDFile(rsid, effectAllele, otherAllele):
    # --- Define the Path to the Output File ---
    outputPath = os.path.join(
        projectRoot,
        'Scripts/Publication/supportingFiles/replication/allelePhasingOutput'
    )
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    outputPath += f"/{rsid}_{rsid}.log"

    #Create a dummy file with the R2 value of 1 and allele annotations
    with open(outputPath, 'w') as file:
        file.write('R-sq = 1')
        file.write('\n')
        file.write('In phase alleles are ' + effectAllele + effectAllele + '/' + otherAllele + otherAllele)

    return 1

#Extract reported gene names in the study as they do not have UniprotIDs annotated
replicationGeneList = SurapaneniDf['EntrezGeneSymbol'].tolist() + SurapaneniDf['Target'].tolist()
replicationGeneList = list(set(replicationGeneList))

for index in df.index:

    geneNameVIKING = df.loc[index, 'HUGO']
    rsidVIKING = df.loc[index, 'SNP']
    effectAlleleVIKING = df.loc[index, 'Effect allele']
    otherAlleleVIKING = df.loc[index, 'Other allele']

    #Find all indexes in Surapaneni with the gene name
    geneIndexList = SurapaneniDf[SurapaneniDf['EntrezGeneSymbol'] == geneNameVIKING].index.tolist() + SurapaneniDf[SurapaneniDf['Target'] == geneNameVIKING].index.tolist()
    geneIndexList = list(set(geneIndexList))

    #Remove indexes that are for trans associations.
    #Manually annotated cis associations are kept where they were erroneously annotated as trans in Surapaneni
    cisDict = {'C8A|C8B|C8G': 'rs12067507', 'ITGA4|ITGB1': 'rs1143676', 'HLA-C': 'rs3094680', 'MICA': 'rs1063631', 'MICB': 'rs3134900', 'C4A|C4B': 'rs7774197', 'HLA-DQA2': 'rs35062987', 'TAPBP': 'rs34132052', 'SCARF1': 'rs2272011', 'CCL18': 'rs2015086', 'IL12A|EBI3': 'rs4740', 'CGA|LHB': 'rs78537284', 'CGA|CGB3|CGB7': 'rs78537284', 'LILRA6': 'rs34810796', 'LILRB5': 'rs12975366', 'LILRB2': 'rs383369', 'LILRA3': 'rs103294', 'LILRB1': 'rs10427127', 'KIR2DL5A': 'rs625698', 'GSTT1': 'rs140196'}
    geneIndexList = [index for index in geneIndexList if SurapaneniDf.loc[index, 'Cis/trans'] == 'cis' or cisDict.get(SurapaneniDf.loc[index, 'EntrezGeneSymbol'], None) == SurapaneniDf.loc[index, 'rsid']]

    #Some require special treatment if they have multiple independent pQTL in Surapaneni
    multipleCisDict = {'VSTM1': ['rs2433724', 'rs612529']}
    if geneNameVIKING in multipleCisDict:
        for rsid in multipleCisDict[geneNameVIKING]:
            geneIndexList += SurapaneniDf[(SurapaneniDf['rsid'] == rsid) & (SurapaneniDf['EntrezGeneSymbol'] == geneNameVIKING)].index.tolist()

    #Check if this protein was reported by annotated gene name
    #Also do not check trans associations reported in VIKING
    if geneNameVIKING not in replicationGeneList or df.loc[index, 'cisTrans'] == 'Trans' or geneIndexList == []:

        if df.loc[index, 'cisTrans'] == 'Cis':
            print('Not replicated', geneNameVIKING, df.loc[index, 'Chromosome'], df.loc[index, 'Position'])

        df.at[index, 'Replication AASK rsid'] = 'NA'
        df.at[index, 'Replication AASK LD'] = 'NA'
        df.at[index, 'Replication AASK p-value'] = 'NA'
        df.at[index, 'Replication AASK Effect Size'] = 'NA'
        df.at[index, 'Replication AASK Allele Frequency'] = 'NA'
        df.at[index, 'Replication AASK SomamerID'] = 'NA'
        continue
    
    resultRsid = []
    resultLD = []
    resultPValue = []
    resultEffectSize = []
    resultAlleleFrequency = []
    resultSomamerID = []

    for surapaneniIndex in geneIndexList:
        rsidSurapaneni = SurapaneniDf.loc[surapaneniIndex, 'rsid']
        chromosome = df.loc[index, 'Chromosome']

        #Skip if the variant is monoallelic in EUR
        monoallelic = ['rs115671497', 'rs12144939', 'rs61034679', 'rs541505191', 'rs74363294', 'rs634013', 'rs116416318', 'rs143283905']
        if rsidSurapaneni in monoallelic:
            continue

        #Skip if the variant is not in UKBB
        notInPanel = ['rs396991', 'rs80045772', 'rs35904956', 'rs3775291', 'rs397723447', 'rs4759815', 'rs625698', 'rs35078252', 'rs11302931']
        if rsidSurapaneni in notInPanel:
            continue
        
        #Some variants are not assigned an rsid in Surapaneni. Full summary stats were downloaded for allele data and rsid mapped with dbsnp where possible
        noRsid = {'FOLR3': '', 'JAML': '', 'OAF': '', 'GP6': ''}
        if geneNameVIKING in noRsid:
            rsidSurapaneni = noRsid[geneNameVIKING]

            #Skip if the rsid was not possible to be determined
            if rsidSurapaneni == '':
                continue

        if rsidSurapaneni != rsidVIKING:
            #Check if the reported rsid is a proxy for the rsid in VIKING
            LD = checkIfProxy(rsidVIKING, rsidSurapaneni, chromosome)
            LD = checkWhatProxy(rsidVIKING, rsidSurapaneni, chromosome)
        else:
            #Create a dummy file as PLINK does not work with the same rsid
            LD = createDummyLDFile(rsidVIKING, effectAlleleVIKING, otherAlleleVIKING)

        #Check if the alleles are flipped between the two datasets and flips the reported beta and allele frequency if necessary
        betaSurapaneni, alleleFrequencySurapaneni = checkAlleleFlip(rsidVIKING, effectAlleleVIKING, otherAlleleVIKING, rsidSurapaneni, SurapaneniDf.loc[surapaneniIndex, 'Coded allele'], SurapaneniDf.loc[surapaneniIndex, 'Beta'], SurapaneniDf.loc[surapaneniIndex, 'CAF'], chromosome)

        #Some rsid even though with proxies could not be processed this way due to mismatch in reported alleles
        if betaSurapaneni == 'NA':
            continue

        pValueSurapaneni = SurapaneniDf.loc[surapaneniIndex, 'P-value']
        somamerID = '-'.join(SurapaneniDf.loc[surapaneniIndex, 'seqid_in_sample'].split('_')[1:])
    
        resultRsid.append(rsidSurapaneni)
        resultLD.append(str(LD))
        resultPValue.append(str(pValueSurapaneni))
        resultEffectSize.append(str(betaSurapaneni))
        resultAlleleFrequency.append(str(alleleFrequencySurapaneni))
        resultSomamerID.append(somamerID)

    if len(resultRsid) == 0:
        df.at[index, 'Replication AASK rsid'] = 'NA'
        df.at[index, 'Replication AASK LD'] = 'NA'
        df.at[index, 'Replication AASK p-value'] = 'NA'
        df.at[index, 'Replication AASK Effect Size'] = 'NA'
        df.at[index, 'Replication AASK Allele Frequency'] = 'NA'
        df.at[index, 'Replication AASK SomamerID'] = 'NA'
        continue

    df.at[index, 'Replication AASK rsid'] = '|'.join(resultRsid)
    df.at[index, 'Replication AASK LD'] = '|'.join(resultLD)
    df.at[index, 'Replication AASK p-value'] = '|'.join(resultPValue)
    df.at[index, 'Replication AASK Effect Size'] = '|'.join(resultEffectSize)
    df.at[index, 'Replication AASK Allele Frequency'] = '|'.join(resultAlleleFrequency)
    df.at[index, 'Replication AASK SomamerID'] = '|'.join(resultSomamerID)


#Filter out trans associations
df = df[df['cisTrans'] == 'Cis']

#Leave only the necessary columns
df = df[['HUGO', 'SNP', 'Effect allele', 'Other allele', 'Effect allele frequency in VIKING', 'Effect Size (beta)', 'P-value', 'somamerID', 'Replication AASK rsid', 'Replication AASK LD', 'Replication AASK p-value', 'Replication AASK Effect Size', 'Replication AASK Allele Frequency', 'Replication AASK SomamerID']]

#Rename columns
df.columns = ['HUGO', 'VIKING rsid', 'Phased Effect allele', 'Phased Other allele', 'VIKING Effect Allele Frequency', 'VIKING Effect Size', 'VIKING p-value', 'VIKING SomamerID', 'Replication AASK rsid', 'Replication AASK LD', 'Replication AASK p-value', 'Replication AASK Effect Size', 'Replication AASK Effect Allele Frequency', 'Replication AASK SomamerID']

df.to_excel(outputDf, index = False)

def calculateReplicationStats():

    def calculateSurapaneniCis():
        #Calculate the number of proteins with cis associations
        cisAssociationsDf = SurapaneniDf[SurapaneniDf['Cis/trans'] == 'cis']
        result = cisAssociationsDf['EntrezGeneSymbol'].to_list()

        #Add the wrongly annotated trans associations
        result += cisDict.keys()

        print('Number of proteins with genome-wide significant cis pQTL in Surapaneni:', len(list(set(result))))
        return result

    def calculateCisOverlap(cisReplication):
        #Retrieve the proteins with cis associations in VIKING
        vikingCis = VIKINGdf[VIKINGdf['cisTrans'] == 'Cis']['HUGO'].to_list()

        #Find the overlap between the two datasets
        result = list(set(cisReplication).intersection(vikingCis))

        print('Number of proteins with cis pQTL in Surapaneni and VIKING:', len(result))

        return result

    def calculateLDpass(LDThresholdR2 = 0.6):
        
        result = 0

        #List for not double counting proteins
        proteinsAnalysed = []

        for index in df.index:
            
            if df.loc[index, 'HUGO'] in proteinsAnalysed:
                continue

            proteinsAnalysed.append(df.loc[index, 'HUGO'])
            
            if df.loc[index, 'Replication AASK LD'] == 'NA':
                continue

            LDvalues = df.loc[index, 'Replication AASK LD'].split('|')

            if any(float(LD) >= LDThresholdR2 for LD in LDvalues):
                result += 1

        print('Number of proteins with cis pQTL in Surapaneni and VIKING with LD R2 >=', LDThresholdR2, ':', result)
        return result

    def calculateTotalReplicatedcispQTL():
        vikingCis = VIKINGdf[VIKINGdf['cisTrans'] == 'Cis'].copy()

        totalCisCount = len(vikingCis)

        #Create a column to track if replicated
        vikingCis['Replicated'] = 'No'

        for index in vikingCis.index:
            #Check if replicated in Pietzner
            if vikingCis.loc[index, 'Replication p-value'] == vikingCis.loc[index, 'Replication p-value']:
                vikingCis.at[index, 'Replicated'] = 'Yes'
                continue
            
            #Check if replicated in Surapaneni
            rsid = vikingCis.loc[index, 'SNP']
            HUGO = vikingCis.loc[index, 'HUGO']

            #Get rows with the same gene name and rsid
            targetIndex = df[(df['HUGO'] == HUGO) & (df['VIKING rsid'] == rsid)].index.tolist()[0]

            if df.loc[targetIndex, 'Replication AASK SomamerID'] == 'NA':
                continue

            replicationSomamerIDList = df.loc[targetIndex, 'Replication AASK SomamerID'].split('|')
            replicationLDList = df.loc[targetIndex, 'Replication AASK LD'].split('|')

            for i in range(len(replicationSomamerIDList)):
                if replicationSomamerIDList[i] == df.loc[targetIndex, 'VIKING SomamerID']:
                    if float(replicationLDList[i]) >= 0.6:
                        vikingCis.at[index, 'Replicated'] = 'Yes'
                        break
        
        #Count the number of replicated pQTL
        replicatedCisCount = len(vikingCis[vikingCis['Replicated'] == 'Yes'])

        print('Number of replicated cis pQTL in VIKING:', replicatedCisCount, 'out of', totalCisCount)
        
        # --- Define Output File Path ---
        outputFile = os.path.join(
            projectRoot,
            'Scripts/Publication/supportingFiles/test/supplementary1ReplicatedCis.xlsx'
        )

        # --- Export the Replicated pQTL ---
        vikingCis.to_excel(outputFile, index=False)
        return 

    cisReplicationProteins = calculateSurapaneniCis()
    cisMatchedProteins = calculateCisOverlap(cisReplicationProteins)
    calculateLDpass()
    calculateTotalReplicatedcispQTL()

    return

calculateReplicationStats()