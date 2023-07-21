from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import itertools
import time
import re
import pandas as pd
import logging
from enum import Enum, auto

class Bcolors(Enum):
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class GeneStruct():
    def __init__(self, name, location, product):
        self.name = name
        self.location = location
        self.product = product

    @classmethod
    def stripGeneEntries(cls, name, location, product):
        if location.start + 1 >= location.end:
            location = str(FeatureLocation(location.start, location.end))
        else:
            location = str(FeatureLocation(location.start + 1, location.end))
        caractersToReplace = ["[", "]", "'"]
        for caracter in caractersToReplace:
            location = location.replace(caracter, "")
            if product:
                product = product.replace(caracter, "")
        return cls(name, location, product)

class FeatureStruct():
    def __init__(self, previousGene, controlRegion, nextGene):
        self.previousGene = previousGene
        self.controlRegion = controlRegion
        self.nextGene = nextGene

    @staticmethod
    def extractGeneInfoToOneLineStr(gene):
        name = str(gene.type)
        location = gene.location
        if "product" in gene.qualifiers:
            product = str(gene.qualifiers["product"])
        else:
            product = None
        return GeneStruct.stripGeneEntries(name, location, product)

    @classmethod
    def initFromSeqFeatureClass(cls, geneList):
        mapObject = map(cls.extractGeneInfoToOneLineStr, geneList)
        geneList = list(mapObject)
        previousGene = geneList[0]
        controlRegion = geneList[1]
        nextGene = geneList[2]
        return cls(previousGene, controlRegion, nextGene)

class RowStruct():
    def __init__(self, previousFeature, controlRegion, nextFeature, accesionName, taxonomy, organism):
        self.previousFeature = previousFeature
        self.controlRegion = controlRegion
        self.nextFeature = nextFeature
        self.accesionName = accesionName
        self.controlRegion = controlRegion
        self.taxonomy = taxonomy
        self.organism = organism

    def __str__(self) -> str:
        return f"{self.accesionName}\t{self.organism}\t{self.taxonomy}\t\
        {self.previousFeature}\t{self.controlRegion}\t{self.nextFeature}"

class GeneBankReader():
    def __init__(self):
        self.logger = self._initLogger()

    def _initLogger(self):
        logging.basicConfig(format="%(asctime)s-%(levelname)s-%(filename)s:%(lineno)d - %(message)s", level=logging.DEBUG)
        return logging.getLogger(__name__)

    def getAnnotations(self, annotations):
        # taxonomy = str(annotations["taxonomy"]).replace("[", "").replace("]", "")
        return f"{annotations['accessions'][0]}.{annotations['sequence_version']}", annotations["organism"], annotations["taxonomy"]

    def pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def getSurroundingGeneIndexes(self, controlRegionIndex):
        if controlRegionIndex == -1:
            nextGeneIndex = 1
            previousGeneIndex = controlRegionIndex - 1
        elif controlRegionIndex == 1:
            previousGeneIndex = -1
            nextGeneIndex = controlRegionIndex + 1
        else:
            nextGeneIndex = controlRegionIndex + 1
            previousGeneIndex = controlRegionIndex - 1
        return (previousGeneIndex, nextGeneIndex)

    def ifNotCorrectGeneType(self, geneType):
        return geneType == "source" or geneType == "gene" or geneType == "misc_feature" or geneType == "repeat_region" or  geneType == "D-loop"

    def ifLocationNotInControlRegion(self, geneLocation, crLocation):
        return geneLocation[0] not in crLocation and geneLocation[-1] not in crLocation

    def getFirstGene(self, features):
        for gene in features:
            if not self.ifNotCorrectGeneType(gene.type):
                return gene
        if len(features) > 1:
            return features[1]
        else:
            return features[0]

    def getLastGene(self, features):
        lastGeneIndex = len(features) - 1
        newGeneIndex = lastGeneIndex
        for number in range(21):
            newGeneIndex = newGeneIndex - number
            if not self.ifNotCorrectGeneType(features[newGeneIndex].type):
                return features[newGeneIndex]
        return features[lastGeneIndex]

    def getNextGeneNoLocation(self, nextGeneIndex, features):
        originalNextGeneIndex = nextGeneIndex
        for number in range(11):
            nextGeneIndex = originalNextGeneIndex + number
            if nextGeneIndex >= len(features):
                nextGene = self.getFirstGene(features)
                return nextGene
            nextGene = features[nextGeneIndex]
            if not self.ifNotCorrectGeneType(nextGene.type):
                return nextGene
        return features[originalNextGeneIndex]

    def getPreviousNoLocation(self, previousGeneIndex, features):
        originalPreviousGeneIndex = previousGeneIndex
        for number in range(11):
            previousGeneIndex = originalPreviousGeneIndex - number
            if previousGeneIndex < 1:
                previousGene = self.getLastGene(features)
                return previousGene
            previousGene = features[previousGeneIndex]
            if not self.ifNotCorrectGeneType(previousGene.type):
                return previousGene
        return features[originalPreviousGeneIndex]

    def getNextGeneWithLocation(self, controlRegionLocation, nextGeneIndex, features):
        rangeOfLocationsCR = range(controlRegionLocation[0], controlRegionLocation[-1])
        originalNextGeneIndex = nextGeneIndex
        for number in range(11):
            nextGeneIndex = originalNextGeneIndex + number
            if nextGeneIndex >= len(features):
                nextGene = self.getFirstGene(features=features)
                return nextGene
            nextGene = features[nextGeneIndex]
            if self.ifNotCorrectGeneType(nextGene.type):
                continue
            nextGeneLocationRange = list(map(int, re.findall(r'\d+', str(nextGene.location))))
            if len(nextGeneLocationRange) > 2 :
                continue
            if self.ifLocationNotInControlRegion(nextGeneLocationRange, rangeOfLocationsCR):
                return nextGene
            else:
                continue
        return self.getNextGeneNoLocation(originalNextGeneIndex, features)

    def getPreviousGeneWithLocation(self, controlRegionLocation, previousGeneIndex, features):
        rangeOfLocationsCR = range(controlRegionLocation[0], controlRegionLocation[-1])
        originalPreviousGeneIndex = previousGeneIndex
        for number in range(11):
            previousGeneIndex = originalPreviousGeneIndex - number
            if previousGeneIndex < 1:
                previousGene = self.getLastGene(features=features)
                return previousGene
            previousGene = features[previousGeneIndex]
            if self.ifNotCorrectGeneType(previousGene.type):
                continue
            previousGeneLocationRange = list(map(int, re.findall(r'\d+', str(previousGene.location))))
            if len(previousGeneLocationRange) > 2 :
                continue
            if self.ifLocationNotInControlRegion(previousGeneLocationRange, rangeOfLocationsCR):
                return previousGene
            else:
                continue
        return self.getNextGeneNoLocation(originalPreviousGeneIndex, features)

    def ifCrInGeneQualifiers(self, gene):
        for featureValue in gene.qualifiers.values():
            if "control region" in featureValue[0] or "D-loop" in featureValue[0]:
                return True
        return False

    def getFeatures(self, features):
        listOfFeatureStructs = []
        for gene in features:
            geneIndex = features.index(gene)
            previousGeneIndex, nextGeneIndex = self.getSurroundingGeneIndexes(geneIndex)
            controlRegionLocation = list(map(int, re.findall(r'\d+', str(gene.location))))
            controlRegionLocation[0] += 1
            if "control region" in gene.type or "D-loop" in gene.type or "C_region" in gene.type:
                if len(controlRegionLocation) > 2:
                    previousGene = self.getPreviousNoLocation(previousGeneIndex, features)
                    nextGene = self.getNextGeneNoLocation(nextGeneIndex, features)
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
                    self.logger.info(f"Previous Gene - {previousGene.type}; CR - {gene.type}; Next Gene - {nextGene.type} added")
                else:
                    previousGene = self.getPreviousGeneWithLocation(controlRegionLocation, previousGeneIndex, features)
                    nextGene = self.getNextGeneWithLocation(controlRegionLocation, nextGeneIndex, features)
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
                    self.logger.info(f"Previous - {previousGene.type}; CR - {gene.type}; Next - {nextGene.type} added")
            else:
                if self.ifCrInGeneQualifiers(gene):
                    previousGene = self.getPreviousGeneWithLocation(controlRegionLocation, previousGeneIndex, features)
                    nextGene = self.getNextGeneWithLocation(controlRegionLocation, nextGeneIndex, features)
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
                    self.logger.info(f"Previous - {previousGene.type}; CR - {gene.type}; Next - {nextGene.type} added")
        return listOfFeatureStructs

    def extractGeneInfoToOneLineStr(self, gene):
        name = gene.type
        location = gene.location
        if "product" in gene.qualifiers:
            product = gene.qualifiers["product"]
        else:
            product = None
        return f"{name} {location} {product}"

    def readGeneBankFile(self, geneBankPath):
        geneDict = {"ACCESSION": [],
                        "ORGANISM": [],
                        "TAXONOMY": [],
                        "PREVIOUS_GENE_NAME": [],
                        "PREVIOUS_GENE_LOCATION": [],
                        "PREVIOUS_GENE_NAME_PRODUCT": [],
                        "CONTROL_REGION_NAME": [],
                        "CONTROL_REGION_LOCATION": [],
                        "NEXT_GENE_NAME": [],
                        "NEXT_GENE_LOCATION": [],
                        "NEXT_GENE_NAME_PRODUCT": []}
        with open(geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            for record in record.records:
                featureStructList = self.getFeatures(record.features)
                if len(featureStructList) == 0:
                    continue
                accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                for featureStruct in featureStructList:
                    geneDict["ACCESSION"].append(accessions)
                    geneDict["ORGANISM"].append(organism)
                    geneDict["TAXONOMY"].append(taxonomy)
                    geneDict["PREVIOUS_GENE_NAME"].append(featureStruct.previousGene.name)
                    geneDict["PREVIOUS_GENE_LOCATION"].append(featureStruct.previousGene.location)
                    geneDict["PREVIOUS_GENE_NAME_PRODUCT"].append(featureStruct.previousGene.product)
                    geneDict["CONTROL_REGION_NAME"].append(featureStruct.controlRegion.name)
                    geneDict["CONTROL_REGION_LOCATION"].append(featureStruct.controlRegion.location)
                    geneDict["NEXT_GENE_NAME"].append(featureStruct.nextGene.name)
                    geneDict["NEXT_GENE_LOCATION"].append(featureStruct.nextGene.location)
                    geneDict["NEXT_GENE_NAME_PRODUCT"].append(featureStruct.nextGene.product)
        geneDf = pd.DataFrame.from_dict(geneDict)
        geneDf = geneDf.set_index("ACCESSION")
        self.logger.info(Bcolors.OKGREEN.value + f"Data Frame out of {geneBankPath.split(r'/')[-1]} created" + Bcolors.ENDC.value)
        return geneDf

    def mergeTwoDataFrames(self, df1, df2):
        return pd.concat([df1, df2])

    def saveDataFrameToFile(self, savePath, geneDf):
        geneDf.to_excel(savePath)
        self.logger.info(Bcolors.OKGREEN.value + f"File {savePath.split(r'/')[-1]} successfully saved" + Bcolors.ENDC.value)

    def saveDataFrameToCsv(self, savePath, geneDf):
        geneDf.to_csv(savePath)

    def getDataFrameOfGenemesWithoutCR(self, geneBankPath):
        geneDict = {"ACCESSION": [], "ORGANISM": [], "TAXONOMY": []}
        recordCountWithCR = 0
        recordCountTotal = 0
        recordCountNotCR = 0
        with open(geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            for record in record.records:
                recordCountTotal += 1
                features = record.features
                hasCR = False
                for gene in features:
                    if "control region" in gene.type or "D-loop" in gene.type or "C_region" in gene.type:
                        hasCR = True
                        recordCountWithCR += 1
                        break
                    else:
                        for featureValue in gene.qualifiers.values():
                            if "control region" in featureValue[0] or "D-loop" in featureValue[0] or "C_region" in featureValue[0]:
                                hasCR = True
                                recordCountWithCR += 1
                                break
                            else:
                                continue
                if hasCR:
                    continue
                else:
                    recordCountNotCR += 1
                    accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                    geneDict["ACCESSION"].append(accessions)
                    geneDict["ORGANISM"].append(organism)
                    geneDict["TAXONOMY"].append(taxonomy)
        geneDf = pd.DataFrame.from_dict(geneDict)
        geneDf = geneDf.set_index("ACCESSION")
        recordStats = (f"Total Record Cournt: {recordCountTotal}", f"Record Cournt With No D-loop: {recordCountNotCR}", f"Record Count With D-loop: {recordCountWithCR}")
        self.logger.info(Bcolors.OKGREEN.value + f"Data Frame out of {geneBankPath.split(r'/')[-1]} created" + Bcolors.ENDC.value)
        return geneDf, recordStats

    def normalizeNames(self, data_frame):
        new_data_frame = data_frame.replace("small subunit ribosomal RNA", "s-rRNA").replace("12S ribosomal RNA", "s-rRNA")
        new_data_frame = new_data_frame.replace("large subunit ribosomal RNA", "l-rRNA").replace("16S ribosomal RNA", "l-rRNA")
        return new_data_frame

if __name__ == "__main__":
    geneBankReader = GeneBankReader()
    mitochondrion1_df = geneBankReader.readGeneBankFile("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff")
    mitochondrion2_df = geneBankReader.readGeneBankFile("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    mergedDf = geneBankReader.mergeTwoDataFrames(mitochondrion1_df, mitochondrion2_df)
    mergedDf = geneBankReader.normalizeNames(mergedDf)
    geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_mitochondrion.xlsx", mergedDf)

    mitochondrion1_df_withoutCR = geneBankReader.getDataFrameOfGenemesWithoutCR("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff")
    mitochondrion2_df_withoutCR = geneBankReader.getDataFrameOfGenemesWithoutCR("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")

    print(mitochondrion1_df_withoutCR[1], "mit1 stats")
    print(mitochondrion2_df_withoutCR[1], "mit2 stats")
    mergedDf_no_CR = geneBankReader.mergeTwoDataFrames(mitochondrion1_df_withoutCR[0], mitochondrion2_df_withoutCR[0])
    geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/results/no_CR.xlsx", mergedDf_no_CR)

    # TESTING
    # mitochondrion1_df = geneBankReader.readGeneBankFile("/home/rszczygielski/bioinf/magisterka/geneBank/test/new_seq_for_testing.gbff")
    # geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/results/TEST_seq.xlsx", mitochondrion1_df)
