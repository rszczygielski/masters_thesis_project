from Bio import GenBank
from Bio import SeqIO
import itertools

class FeatureStruct():
    def __init__(self, name, locasion, qualifiers):
        self.name = name
        self.location = locasion
        self.qualifiers = qualifiers

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
            return self.extractGeneInfoToOneLineStr(previousGene), self.extractGeneInfoToOneLineStr(controlRegion), self.extractGeneInfoToOneLineStr(nextGene)
        else:
            return None, None, None

    def extractGeneInfoToOneLineStr(self, gene):
        name = gene.type
        location = gene.location
        if "product" in gene.qualifiers:
            product = gene.qualifiers["product"]
        else:
            product = None
        return f"{name} {location} {product}"

    def readGeneBankFile(self):
        with open(self.geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            for record in record.records:
                previousGene, controlRegion, nextGene = self.getFeatures(record.features)
                if previousGene == None:
                    continue
                accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                rowStruct = RowStruct(previousGene, controlRegion, nextGene, accessions, taxonomy, organism)
                self.saveRowToFile(rowStruct.__str__())
    
    def saveRowToFile(slef, row):
        with open("Organisms_CR.txt", "a") as saveFile:
            saveFile.write("ACCESSION\tORGANISM\tTAXONOMY\tPREVIOUS_GENE\tCONTROL_REGION\tNEXT_GENE\n")
            saveFile.write(f"{row}\n")
                

if __name__ == "__main__":
    geneBankReader = GeneBankReader("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    geneBankReader.readGeneBankFile()