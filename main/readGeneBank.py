from Bio import SeqIO
import itertools
import re
import pandas as pd

class GeneStruct():
    def __init__(self, name, location, product):
        self.name = name
        self.location = location
        self.product = product

class FeatureStruct():
    def __init__(self, previousGene, controlRegion, nextGene):
        self.previousGene = previousGene
        self.controlRegion = controlRegion
        self.nextGene = nextGene

    @staticmethod
    def extractGeneInfoToOneLineStr(gene):
        name = str(gene.type)
        location = str(gene.location)
        if "product" in gene.qualifiers:
            product = str(gene.qualifiers["product"])
        else:
            product = None
        return GeneStruct(name, location, product)

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

    def getAnnotations(self, annotations):
        return f"{annotations['accessions'][0]}.{annotations['sequence_version']}", annotations["organism"], annotations["taxonomy"]

    def pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def getOneOfTenPreviosGenes(self, controlRegionLocation, features, previousGeneIndex, controlRegion):
        if len(controlRegionLocation) > 2: #chceking if there are multiple locations asigned to control region
            if features[previousGeneIndex].type == "repeat_region" or features[previousGeneIndex].type == "source" or features[previousGeneIndex].type == "D-loop":
                previousGeneIndex  -= 1
            previousGene = features[previousGeneIndex]
            return previousGene
        range_of_locations = range(controlRegionLocation[0]+1, controlRegionLocation[-1])
        for number in range(11): #looping through 10 previos genes and comparing locations
            previousGeneIndex = previousGeneIndex - number
            previousGene = features[previousGeneIndex]
            previousGeneLocation = list(map(int, re.findall(r'\d+', str(previousGene.location))))
            if len(previousGeneLocation) > 2 :
                continue
            if previousGeneLocation[0] not in range_of_locations and previousGeneLocation[-1] not in range_of_locations:
                if previousGene.type == "source" or previousGene.type == "gene" or previousGene.type == "repeat_region" or previousGene.type == "D-loop":
                    continue
                return previousGene
            else:
                continue
        previousGene = features[previousGeneIndex + 10]
        return previousGene

    def getPreviousGene(self, features, controlRegion):
        controlRegionIndex = features.index(controlRegion)
        if controlRegionIndex == 1:
            previousGeneIndex = -1
        else:
            previousGeneIndex = controlRegionIndex - 1
        controlRegionLocation = list(map(int, re.findall(r'\d+', str(controlRegion.location))))
        previousGene = self.getOneOfTenPreviosGenes(controlRegionLocation, features, previousGeneIndex, controlRegion)
        # if previousGene.type == "misc_feature":
            # print(previousGene)
            # print(controlRegion)
        return previousGene

    def getOneOfTenNextGenes(self, controlRegionLocation, features, nextGeneIndex, controlRegion):
        if len(controlRegionLocation) > 2:
            if features[nextGeneIndex].type == "repeat_region" or features[nextGeneIndex].type == "source" or features[nextGeneIndex].type == "D-loop":
                nextGeneIndex  += 1
            nextGene = features[nextGeneIndex]
            return nextGene
        range_of_locations = range(controlRegionLocation[0]+1, controlRegionLocation[-1])
        for number in range(11):
            nextGeneIndex = nextGeneIndex + number
            if nextGeneIndex >= len(features):
                nextGeneIndex = 1
            nextGene = features[nextGeneIndex]
            nextGeneLocation = list(map(int, re.findall(r'\d+', str(nextGene.location))))
            if len(nextGeneLocation) > 2 :
                continue
            if nextGeneLocation[0] not in range_of_locations and nextGeneLocation[-1] not in range_of_locations:
                if nextGene.type == "source" or nextGene.type == "gene" or nextGene.type == "repeat_region" or  nextGene.type == "D-loop":
                    continue
                return nextGene
            else:
                continue
        nextGene = features[nextGeneIndex - 10]
        return nextGene

    def getNextGene(self, features, controlRegion):
        controlRegionIndex = features.index(controlRegion)
        if controlRegionIndex == -1:
            nextGeneIndex = 1
        else:
            nextGeneIndex = controlRegionIndex + 1
        controlRegionLocation = list(map(int, re.findall(r'\d+', str(controlRegion.location))))
        nextGene = self.getOneOfTenNextGenes(controlRegionLocation, features, nextGeneIndex, controlRegion)
        return nextGene

    def getFeatures(self, features):
        listOfFeatureStructs = []
        for gene in features:
            if "control region" in gene.type or "D-loop" in gene.type or "C_region" in gene.type:
                previousGene = self.getPreviousGene(features, gene)
                nextGene = self.getNextGene(features, gene)
                listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, gene, nextGene]))
            else:
                for featureValue in gene.qualifiers.values():
                    if "source" in featureValue or "control region" in featureValue or "D-loop" in featureValue:
                        # print(gene)
                        previousGene = self.getPreviousGene(features, gene)
                        nextGene = self.getNextGene(features, gene)
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


if __name__ == "__main__":
    geneBankReader = GeneBankReader()
    mitochondrion1_df = geneBankReader.readGeneBankFile("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff")
    mitochondrion2_df = geneBankReader.readGeneBankFile("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    mergedDf = geneBankReader.mergeTwoDataFrames(mitochondrion1_df, mitochondrion2_df)
    geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion_2.xlsx", mitochondrion2_df)
    geneBankReader.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/main_mitochondrion.xlsx", mergedDf)