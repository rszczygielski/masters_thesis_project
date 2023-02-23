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
        self.accesionName = accesionName[0]
        self.controlRegion = controlRegion
        self.taxonomy = taxonomy
        self.organism = organism

    def __str__(self) -> str:
        return f"{self.accesionName}\t{self.organism}\t{self.taxonomy}\t\
        {self.previousFeature}\t{self.controlRegion}\t{self.nextFeature}"

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

    def getNextGenFromFeatures(self, features, controlRegionIndex):
        if len(features) - 1 == controlRegionIndex:
            nextGene = features[0]
            if features[0].type == "source":
                nextGene = features[1]
                if features[1].type == "gene":
                    nextGene = features[2]
        else:
            nextGene = features[controlRegionIndex + 1]
        return nextGene

    def getFeatures(self, features):
        listOfFeatureStructs = []
        for previousGene, controlRegion in self.pairwise(features):
            controlRegionIndex = None
            if "control region" in controlRegion.type or "D-loop" in controlRegion.type:
                controlRegionIndex = features.index(controlRegion)
                nextGene =  self.getNextGenFromFeatures(features, controlRegionIndex)
                listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, controlRegion, nextGene]))
                continue
            for featureValue in controlRegion.qualifiers.values():
                if "control region" in featureValue or "D-loop" in featureValue:
                    controlRegionIndex = features.index(controlRegion)
                if controlRegionIndex:
                    nextGene =  self.getNextGenFromFeatures(features, controlRegionIndex)
                    listOfFeatureStructs.append(FeatureStruct.initFromSeqFeatureClass([previousGene, controlRegion, nextGene]))
        return listOfFeatureStructs

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
                featureStructList = self.getFeatures(record.features)
                if len(featureStructList) == 0:
                    continue
                accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                for featureStruct in featureStructList:
                    rowStruct = RowStruct(featureStruct.previousGene, featureStruct.controlRegion, featureStruct.nextGene, accessions, taxonomy, organism)
                    listOfRows.append(rowStruct.__str__())
        return listOfRows

    def saveRowToFile(self, saveName):
        listOfRows = self.readGeneBankFile()
        with open(f"../results/{saveName}", "w") as saveFile:
            saveFile.write("ACCESSION\tORGANISM\tTAXONOMY\tPREVIOUS_GENE\tCONTROL_REGION\tNEXT_GENE\n")
            for row in listOfRows:
                saveFile.write(f"{row}\n")
                print(f"{row} added to file")

if __name__ == "__main__":
    geneBankReader = GeneBankReader("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff")
    geneBankReader.saveRowToFile("Organisms_mitochondion_1.txt")
    # geneBankReader = GeneBankReader("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    # geneBankReader.saveRowToFile("Organisms_mitochondion_2.txt")

    # pominąć source
    # brac wszystkie regiony kontrolne które są
    # duplikowac wiersz - dodatkowy rekord dla tego samego genomu z kolejnym regionem kontolnym - jak jest więcej niż jeden region

    # posegregować,  wyciągnać nazwy genów
    # sortuję najpierw po genach potem po taksonomii


    # wyjątki w taksonomi - ten sam gen pod inną nazwą
    # Mammalia
    # ile było przypadków oflanowania regionu kontrolnego okreslonymi genami - ile było przypadków kiedy było przypadków region kontrolny a gen jakiś w postaci previous a drugi ile przypadków region kontolny dla next gen
    # tableka zbiorcza - poszczególne jednostek takosnomicznych, a kolumny to kombinacje genów