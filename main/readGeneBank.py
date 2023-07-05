from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import itertools
import time
import re
import pandas as pd

class GeneStruct():
    def __init__(self, name, location, product):
        self.name = name
        self.location = location
        self.product = product

    @classmethod
    def stripGeneEntries(cls, name, location, product):
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
        # print(geneList)
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

    def getAnnotations(self, annotations):
        # taxonomy = str(annotations["taxonomy"]).replace("[", "").replace("]", "")
        return f"{annotations['accessions'][0]}.{annotations['sequence_version']}", annotations["organism"], annotations["taxonomy"]

    def pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def getSelectedPreviousGene(self, controlRegionLocation, features, previousGeneIndex, controlRegion):
        range_of_locations = range(controlRegionLocation[0], controlRegionLocation[-1]) # getting range values for the given location of CR
        # print(range_of_locations, "CR")
        originalPreviousGeneIndex = previousGeneIndex
        for number in range(11): # looping through 10 previos genes and comparing locations
            previousGeneIndex = originalPreviousGeneIndex - number #creating new gene index
            if previousGeneIndex < 1:
                previousGeneIndex = len(features) - 1
                previousGene = features[previousGeneIndex]
                return previousGene # returning gene while crossing borrder line
            previousGene = features[previousGeneIndex]
            if previousGene.type == "source" or previousGene.type == "gene" or previousGene.type == "repeat_region" or previousGene.type == "D-loop":
                    continue # skipping unnecessary genes
            previousGeneLocation = list(map(int, re.findall(r'\d+', str(previousGene.location)))) # getting range values for the given location of gene
            # print(previousGeneLocation, "previous")
            if len(previousGeneLocation) > 2 :
                continue # skipping examples with multiple coordinates
            if previousGeneLocation[0] not in range_of_locations and previousGeneLocation[-1] not in range_of_locations:
                return previousGene # returning the gene when coordiantes are not corssing
            else:
                continue
        return features[originalPreviousGeneIndex] #after checking 10 gens returing the gene which is the closest to CR

    def getSelectedNextGene(self, controlRegionLocation, features, nextGeneIndex, controlRegion):
        range_of_locations = range(controlRegionLocation[0], controlRegionLocation[-1])
        originalNextGeneIndex = nextGeneIndex
        for number in range(11):
            nextGeneIndex = originalNextGeneIndex + number
            if nextGeneIndex >= len(features):
                nextGeneIndex = 1
                nextGene = features[nextGeneIndex]
                return nextGene
            nextGene = features[nextGeneIndex]
            if nextGene.type == "source" or nextGene.type == "gene" or nextGene.type == "repeat_region" or  nextGene.type == "D-loop":
                continue
            nextGeneLocation = list(map(int, re.findall(r'\d+', str(nextGene.location))))
            if len(nextGeneLocation) > 2 :
                continue
            if nextGeneLocation[0] not in range_of_locations and nextGeneLocation[-1] not in range_of_locations:
                return nextGene
            else:
                continue
        return features[originalNextGeneIndex]

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

    def getFeatures(self, features):
        listOfFeatureStructs = []
        for gene in features:
            geneIndex = features.index(gene)
            previousGeneIndex, nextGeneIndex = self.getSurroundingGeneIndexes(geneIndex)
            controlRegionLocation = list(map(int, re.findall(r'\d+', str(gene.location))))
            controlRegionLocation[0] += 1
            if "control region" in gene.type or "D-loop" in gene.type or "C_region" in gene.type:
                if len(controlRegionLocation) > 2:
                    if features[nextGeneIndex].type == "repeat_region" or features[nextGeneIndex].type == "source" or features[nextGeneIndex].type == "D-loop":
                        nextGeneIndex  += 1
                    if features[previousGeneIndex].type == "repeat_region" or features[previousGeneIndex].type == "source" or features[previousGeneIndex].type == "D-loop":
                        previousGeneIndex  -= 1
                    nextGene = features[nextGeneIndex]
                    previousGene = features[previousGeneIndex]
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
                else:
                    previousGene = self.getSelectedPreviousGene(controlRegionLocation, features, previousGeneIndex, gene)
                    nextGene = self.getSelectedNextGene(controlRegionLocation, features, nextGeneIndex, gene)
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
            else:
                for featureValue in gene.qualifiers.values():
                    if "source" in featureValue or "control region" in featureValue or "D-loop" in featureValue:
                        previousGene = self.getSelectedPreviousGene(controlRegionLocation, features, previousGeneIndex, gene)
                        nextGene = self.getSelectedNextGene(controlRegionLocation, features, nextGeneIndex, gene)
                        # print(type(gene.location.start))
                        listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
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
        return geneDf

    def mergeTwoDataFrames(self, df1, df2):
        return pd.concat([df1, df2])

    def saveDataFrameToFile(self, savePath, geneDf):
        geneDf.to_excel(savePath)

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
                    # print(gene.type)
                    if "control region" in gene.type or "D-loop" in gene.type or "C_region" in gene.type:
                        hasCR = True
                        recordCountWithCR += 1
                        break
                    else:
                        for featureValue in gene.qualifiers.values():
                            if "control region" in featureValue or "D-loop" in featureValue:
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
        recordStats = (recordCountTotal, recordCountNotCR, recordCountWithCR)
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

    # TESTING
    mitochondrion1_df_withoutCR = geneBankReader.getDataFrameOfGenemesWithoutCR("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff")
    mitochondrion2_df_withoutCR = geneBankReader.getDataFrameOfGenemesWithoutCR("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")

    mergedDf_no_CR = geneBankReader.mergeTwoDataFrames(mitochondrion1_df_withoutCR[0], mitochondrion2_df_withoutCR[0])
    geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/results/no_CR.xlsx", mergedDf_no_CR)