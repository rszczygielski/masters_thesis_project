from Bio import SeqIO
import itertools

class FeatureStruct():
    def __init__(self, previousGene, controlRegion, nextGene):
        self.previousGene = previousGene
        self.controlRegion = controlRegion
        self.nextGene = nextGene
    
    @staticmethod
    def extractGeneInfoToOneLineStr(gene):
        name = gene.type
        location = gene.location
        if "product" in gene.qualifiers:
            product = gene.qualifiers["product"]
        else:
            product = None
        return f"{name} {location} {product}"

    @classmethod
    def initFromSeqFeatureClass(cls, geneList):
        geneList = map(cls.extractGeneInfoToOneLineStr, geneList)
        geneList = list(geneList)
        previousGene = geneList[0]
        controlRegion = geneList[1]
        nextGene = geneList[2]
        return cls(previousGene, controlRegion, nextGene)

class RowStruct():
    def __init__(self, previousFeature, controlRegion, nextFeature, accesionName, taxonomy, organism, ):
        self.previousFeature = previousFeature
        self.controlRegion = controlRegion
        self.nextFeature = nextFeature
        self.accesionName = accesionName
        self.controlRegion = controlRegion
        self.taxonomy = taxonomy
        self.organism = organism

    def __str__(self) -> str:
        return f"{self.accesionName}\t{self.organism}\t{self.taxonomy}\t\
        {self.previousFeature}\t{self.controlRegion}\t{self.nextFeature}\t"

class GeneBankReader():
    def __init__(self, geneBankPath):
        self.geneBankPath = geneBankPath

    def getAnnotations(self, annotations):
        return annotations["accessions"], annotations["organism"], annotations["taxonomy"]
    
    def pairwise(self, iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def getFeatures(self, features):
        controlRegionIndex = None
        for previousGene, controlRegion in self.pairwise(features):
            for featureValue in controlRegion.qualifiers.values():
                if "control region" in featureValue or "D-loop" in featureValue:
                    controlRegionIndex = features.index(controlRegion)
                    break
        if controlRegionIndex:
            if len(features) - 1 == controlRegionIndex:
                nextGene = features[0]
            else:
                nextGene = features[controlRegionIndex + 1]
            # return self.extractGeneInfoToOneLineStr(previousGene), self.extractGeneInfoToOneLineStr(controlRegion), self.extractGeneInfoToOneLineStr(nextGene)
            return FeatureStruct.initFromSeqFeatureClass([previousGene, controlRegion, nextGene])
        else:
            # return None, None, None
            return None

    def extractGeneInfoToOneLineStr(self, gene):
        name = gene.type
        location = gene.location
        if "product" in gene.qualifiers:
            product = gene.qualifiers["product"]
        else:
            product = None
        return f"{name} {location} {product}"

    def readGeneBankFile(self):
        listOfRows = []
        with open(self.geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            for record in record.records:
                # previousGene, controlRegion, nextGene = self.getFeatures(record.features)
                featureStruct = self.getFeatures(record.features)
                # print(featureStruct.controlRegion)
                if featureStruct == None:
                    continue
                accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                rowStruct = RowStruct(featureStruct.previousGene, featureStruct.controlRegion, featureStruct.nextGene, accessions, taxonomy, organism)
                # print(rowStruct)
                listOfRows.append(rowStruct.__str__())
        return listOfRows

    def saveRowToFile(self):
        listOfRows = self.readGeneBankFile()
        with open("Organisms_CR.txt", "w") as saveFile:
            saveFile.write("ACCESSION\tORGANISM\tTAXONOMY\tPREVIOUS_GENE\tCONTROL_REGION\tNEXT_GENE\n")
            for row in listOfRows:
                saveFile.write(f"{row}\n")
                print(f"{row} added to file")
                
if __name__ == "__main__":
    geneBankReader = GeneBankReader("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    geneBankReader.saveRowToFile()