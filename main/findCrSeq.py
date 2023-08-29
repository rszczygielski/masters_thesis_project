from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import itertools
import time
import re
import pandas as pd
from pandas import DataFrame
import logging
# from IPython.display import display
from enum import Enum, auto
from readGeneBank import GeneStruct
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

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

    def join_gene_data(self, column_for_joining="TAXONOMY"):
        df_inner = pd.merge(self.data_frame_dict["all_previous_gene_mit"], self.data_frame_dict["all_next_gene_mit"],
                            on=column_for_joining, how='inner')
        df_inner = df_inner.drop(columns=['NEXT_GENE_NAME_PRODUCT', 'PREVIOUS_GENE_NAME_PRODUCT',
                                          'FROM_HOW_MANY_TAXONOMY_x','FROM_HOW_MANY_TAXONOMY_y',
                                          'ACCESSION_y']).rename(columns={"ACCESSION_x": "ACCESSION",
                                                                          "LOCATION_x": "PREVIOUS_GENE_LOCATION",
                                                                          "LOCATION_y": "NEXT_GENE_LOCATION"})
        df_inner = df_inner.set_index("ACCESSION")
        return df_inner

    def get_df(self, df_name):
        return self.data_frame_dict[df_name]

    def getAnnotations(self, annotations):
        return f"{annotations['accessions'][0]}.{annotations['sequence_version']}",\
            str(annotations["taxonomy"])

    def get_cr_sequence(self, sequence, cr_location):
        location = list(map(int, re.findall(r'\d+', str(cr_location))))
        location[0] -= 1
        isComplemetent = True if "(-)" in cr_location else False
        sequence_model = SeqFeature(FeatureLocation(location[0], location[1]), type="control_region")
        control_region_sequence = sequence_model.extract(sequence)
        # print(dir(control_region_sequence.defined))
        if not control_region_sequence.defined:
            return False
        if isComplemetent:
         control_region_sequence =  control_region_sequence.complement()
        return str(control_region_sequence)

    def find_all_cr(self, geneBankPath, model_df):
        seq_dict = {}
        with open(geneBankPath) as geneBankFile:
            record = SeqIO.parse(geneBankFile, 'genbank')
            for record in record.records:
                accessions, taxonomy = self.getAnnotations(record.annotations)
                taxonomy = taxonomy.replace("['Eukaryota', 'Metazoa', ", "[")
                list_of_model_accessions = list(model_df["ACCESSION"])
                list_of_no_cr_accessions = list(self.data_frame_dict["no_CR"]["ACCESSION"])
                if accessions not in list_of_model_accessions or accessions in list_of_no_cr_accessions:
                    continue
                model_info_df = model_df.where(model_df["ACCESSION"] == accessions).dropna()
                if len(model_info_df) != 1:
                    continue
                cr_location = model_info_df.CONTROL_REGION_LOCATION
                surrounding_genes = model_info_df.SURROUNDING_PAIRS.values[0]
                control_region_sequence = self.get_cr_sequence(record.seq, cr_location.values[0])
                if not control_region_sequence:
                    continue
                entry_name = f"{accessions}|{taxonomy}|{surrounding_genes}"
                seq_dict[entry_name] = control_region_sequence
        return seq_dict

    # def get_all_grouped_seq(self, sequence, full_taxonomy, grouped_df:DataFrame):
    #     grouped_dict = {}
    #     for accession, row in grouped_df.iterrows():
    #         control_region_sequence = self.get_cr_sequence(sequence, row.CONTROL_REGION_LOCATION)
    #         if not control_region_sequence:
    #                 continue
    #         full_taxonomy = full_taxonomy.replace("['Eukaryota', 'Metazoa', ", "[")
    #         entry_name = f"{accession}\t{full_taxonomy}\t{row.SURROUNDING_PAIRS}"
    #         grouped_dict[entry_name] = control_region_sequence
    #     return grouped_dict

    # def search_for_taxonomy_in_model(self, model_df, full_taxonomy, sequence):
    #     for short_taxonomy in list(model_df["TAXONOMY"]):
    #         if short_taxonomy not in full_taxonomy:
    #             continue
    #         grouped_df = model_df.where(model_df["TAXONOMY"].str.contains(short_taxonomy, regex=False)).dropna()
    #         print(len(grouped_df), "grouped")
    #         grouped_df = grouped_df.set_index("ACCESSION")
    #         sequence_dict = self.get_all_grouped_seq(sequence, full_taxonomy, grouped_df)
    #         print(len(sequence_dict))
    #         short_taxonomy = short_taxonomy.replace("['Eukaryota', 'Metazoa', ", "[")
    #         return sequence_dict, short_taxonomy

    # def find_grouped_cr(self, geneBankPath, model_df):
    #     grouped_seq_dict = {}
    #     with open(geneBankPath) as geneBankFile:
    #         record = SeqIO.parse(geneBankFile, 'genbank')
    #         for record in record.records:
    #             accessions, full_taxonomy = self.getAnnotations(record.annotations)
    #             list_of_model_accessions = list(model_df["ACCESSION"])
    #             list_of_no_cr_accessions = list(self.data_frame_dict["no_CR"]["ACCESSION"])
    #             if accessions not in list_of_model_accessions or accessions in list_of_no_cr_accessions:
    #                 continue
    #             grouped_taxonomy_dict, short_taxonomy = self.search_for_taxonomy_in_model(model_df, full_taxonomy, record.seq)
    #             grouped_seq_dict[short_taxonomy] = grouped_taxonomy_dict
    #     return grouped_seq_dict

    # def read_fasta_to_dict(self, path):
    #     fasta_dict = {}
    #     key_entry = ""
    #     seq = ""
    #     with open(path, "r") as fasta_file:
    #         for line in fasta_file.readlines():
    #             if ">" in line

    def replace_taxonomy_group_for_short_one(self, dict_with_full_taxonomy, short_model_df):
        main_taxonomy_seq_dict = {}
        for full_key, sequence in dict_with_full_taxonomy.items():
            splited_key = full_key.split("|")
            full_taxonomy = splited_key[1]
            # print(full_taxonomy)
            for short_taxonomy in list(short_model_df["TAXONOMY"]):
                if short_taxonomy not in full_taxonomy:
                    continue
                if short_taxonomy not in main_taxonomy_seq_dict:
                    main_taxonomy_seq_dict[short_taxonomy] = [(full_key, sequence)]
                else:
                    main_taxonomy_seq_dict[short_taxonomy].append((full_key, sequence))
                break
        print(main_taxonomy_seq_dict)
        print(len(main_taxonomy_seq_dict))
        return main_taxonomy_seq_dict

    def createNewDirectory(self, path):
        """Method creates a new directory to save all multiple alignments
        """
        if not os.path.isdir(path):
            os.makedirs(path)

    def save_grouped_seq_to_fasta(self, save_direcotry, seq_dict):
        for file_name, taxonomy_sequece_list in seq_dict.items():
            fasta_string = ""
            for sequence_tuple in taxonomy_sequece_list:
                fasta_string += f">{sequence_tuple[0]}\n"
                fasta_string += f"{sequence_tuple[1]}\n"
            file_path = os.path.join(save_direcotry, file_name)
            with open(file_path, "w") as fastaFile:
                fastaFile.write(fasta_string)

    def saveAllCrFasta(self, savePath, seqDict):
        fasta_string_to_save = ""
        for info, sequence in seqDict.items():
            fasta_string_to_save += f">{info}\n"
            fasta_string_to_save += f"{sequence}\n"
        with open(f"{savePath}/all_seq.fasta", "w") as fastaFile:
            fastaFile.write(fasta_string_to_save)
            self.logger.info(Bcolors.OKGREEN.value + f"File {savePath.split(r'/')[-1]} successfully saved" + Bcolors.ENDC.value)
        return

    def mergeTwoDataFrames(self, df1, df2):
        return pd.concat([df1, df2])

if __name__ == "__main__":
    control_region_finder = ControlRegionFinder()
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/no_CR.xlsx")
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_mitochondrion_no_location.xlsx")
    control_region_finder.read_gene_file("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_results_no_location/surrounding_pairs_mit_no_location.xlsx")
    # model_df_joined = control_region_finder.join_gene_data()
    # model_df = control_region_finder.get_df("main_mitochondrion_no_location")
    model_df = control_region_finder.get_df("main_mitochondrion_no_location")
    short_model_df = control_region_finder.get_df("surrounding_pairs_mit_no_location")

    # seqDict = control_region_finder.find_grouped_cr("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff", model_df)
    seqDict = control_region_finder.find_all_cr("/home/rszczygielski/bioinf/magisterka/geneBank/mitochondrion.1.genomic.gbff", model_df)
    control_region_finder.saveAllCrFasta("/home/rszczygielski/bioinf/magisterka/geneBank/results", seqDict)
    grouped_seq_dict = control_region_finder.replace_taxonomy_group_for_short_one(seqDict, short_model_df)
    control_region_finder.save_grouped_seq_to_fasta("/home/rszczygielski/bioinf/magisterka/fasta_files/grouped_taxonomy", grouped_seq_dict)
# wywalić Eucariota i metazoa

