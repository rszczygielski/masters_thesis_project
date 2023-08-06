from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import itertools
import time
import re
import pandas as pd
import logging
# from IPython.display import display
from enum import Enum, auto
from readGeneBank import GeneStruct

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

class ControlRegionFinder():
    def __init__(self) -> None:
        self.data_frame_dict = {}
        self.logger = self._initLogger()

    def _initLogger(self):
        logging.basicConfig(format="%(asctime)s-%(levelname)s-%(filename)s:%(lineno)d - %(message)s",
                            level=logging.DEBUG)
        return logging.getLogger(__name__)

    def read_gene_file(self, file_path):
        gene_data_frame = pd.read_excel(file_path)
        file_name = file_path.split("/")[-1].split(".")[0]
        self.data_frame_dict[file_name] = gene_data_frame
        self.logger.info(Bcolors.OKGREEN.value + f"{file_name} loaded" + Bcolors.ENDC.value)

    def get_joined_data(self, column_for_joining="TAXONOMY"):
        df_inner = pd.merge(self.data_frame_dict["all_previous_gene_mit"], self.data_frame_dict["all_next_gene_mit"],
                            on=column_for_joining, how='inner')
        df_inner = df_inner.drop(columns=['NEXT_GENE_NAME_PRODUCT', 'PREVIOUS_GENE_NAME_PRODUCT',
                                          'FROM_HOW_MANY_TAXONOMY_x','FROM_HOW_MANY_TAXONOMY_y',
                                          'ACCESSION_y']).rename(columns={"ACCESSION_x": "ACCESSION",
                                                                          "LOCATION_x": "PREVIOUS_GENE_LOCATION",
                                                                          "LOCATION_y": "NEXT_GENE_LOCATION"})
        df_inner = df_inner.set_index("ACCESSION")
        return df_inner

    def getAnnotations(self, annotations):
        return f"{annotations['accessions'][0]}.{annotations['sequence_version']}",\
            annotations["organism"], str(annotations["taxonomy"])

    def get_missing_crs(self, features, previous_gene_location, next_gene_location):
        missing_crs_list = []
        previous_gene_location = int(previous_gene_location.split(":")[1])
        next_gene_location = int(next_gene_location.split(":")[0])
        for gene in features:
            if gene.type == "source":
                continue
            gene_location = list(map(int, re.findall(r'\d+', str(gene.location))))
            gene_location[0] += 1
            if next_gene_location > previous_gene_location:
                location_range = range(previous_gene_location, next_gene_location)
                if gene_location[0] in location_range and gene_location[-1] in location_range:
                    missing_crs_list.append(GeneStruct.singleInit(gene))
                    print("previous", previous_gene_location, "next", next_gene_location)
                    print(gene_location[0], gene_location[1], gene.type, GeneStruct.singleInit(gene).product)
            else:
                if len(gene_location) > 2:
                    continue
                if gene_location[1] < next_gene_location or gene_location[0] > previous_gene_location:
                    print("previous", previous_gene_location, "next", next_gene_location)
                    print(gene_location[0], gene_location[1], gene.type, GeneStruct.singleInit(gene).product)
                    missing_crs_list.append(GeneStruct.singleInit(gene))
        return missing_crs_list

    def read_genebank_file(self, geneBankPath):
        model_df = self.get_joined_data()
        geneDict = {"ACCESSION": [],
                    "ORGANISM": [],
                    "TAXONOMY": [],
                    "CONTROL_REGION_NAME": [],
                    "CONTROL_REGION_PRODUCT": [],
                    "CONTROL_REGION_LOCATION": []}
        with open(geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            # print(self.data_frame_dict["no_CR"]["TAXONOMY"])

            for record in record.records:
                accessions, organism, taxonomy = self.getAnnotations(record.annotations)
                if accessions not in list(self.data_frame_dict["no_CR"]["ACCESSION"]):
                    continue
                for model_taxonomy in list(model_df["TAXONOMY"]):
                    if model_taxonomy not in taxonomy:
                        continue
                    # print("TEST")
                    previous_gene_location = model_df.where(model_df["TAXONOMY"] == model_taxonomy).dropna().PREVIOUS_GENE_LOCATION.values[0]
                    # print(previous_gene_location)
                    if len(model_df.where(model_df["TAXONOMY"] == model_taxonomy).dropna().PREVIOUS_GENE_LOCATION.values) > 1:
                        print("TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
                    next_gene_location = model_df.where(model_df["TAXONOMY"] == model_taxonomy).dropna().NEXT_GENE_LOCATION.values[0]
                    if len(model_df.where(model_df["TAXONOMY"] == model_taxonomy).dropna().NEXT_GENE_LOCATION.values) > 1:
                        print("TEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
                    missing_control_region_list = self.get_missing_crs(record.features, previous_gene_location,
                                                                    next_gene_location)
                    for control_region_struct in missing_control_region_list:
                        geneDict["ACCESSION"].append(accessions)
                        geneDict["ORGANISM"].append(organism)
                        geneDict["TAXONOMY"].append(taxonomy)
                        geneDict["CONTROL_REGION_NAME"].append(control_region_struct.name)
                        geneDict["CONTROL_REGION_PRODUCT"].append(control_region_struct.product)
                        geneDict["CONTROL_REGION_LOCATION"].append(control_region_struct.location)
        geneDf = pd.DataFrame.from_dict(geneDict)
        geneDf = geneDf.set_index("ACCESSION")
        self.logger.info(Bcolors.OKGREEN.value + f"Data Frame out of {geneBankPath.split(r'/')[-1]} created" + Bcolors.ENDC.value)
        return geneDf

    def saveDataFrameToFile(self, savePath, geneDf):
        geneDf.to_excel(savePath)
        self.logger.info(Bcolors.OKGREEN.value + f"File {savePath.split(r'/')[-1]} successfully saved" + Bcolors.ENDC.value)

    def mergeTwoDataFrames(self, df1, df2):
        return pd.concat([df1, df2])

if __name__ == "__main__":
    control_region_finder = ControlRegionFinder()
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_results/all_previous_gene_mit.xlsx")
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_results/all_next_gene_mit.xlsx")
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/no_CR.xlsx")

    # print(control_region_finder.data_frame_list[1])
    joined_data = control_region_finder.get_joined_data()
    # print(joined_data.columns)
    # if "NC_049120.1" in joined_data.index:
    #     print(joined_data.loc["NC_049120.1"])
    mit1_df = control_region_finder.read_genebank_file("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    mit2_df = control_region_finder.read_genebank_file("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.2.genomic.gbff")
    merged_df = control_region_finder.mergeTwoDataFrames(mit1_df, mit2_df)
    control_region_finder.saveDataFrameToFile("/home/rszczygielski/bioinf/magisterka/geneBank/results/finded_CR.xlsx", merged_df)
    # print(joined_data)
    # print(joined_data.where(joined_data["ACCESSION_x"] == joined_data["ACCESSION_y"]))



