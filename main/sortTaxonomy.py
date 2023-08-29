import time
from IPython.display import display
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import pandas as pd
import logging
import datetime
from readGeneBank import Bcolors

class SortTaxonomy():
    def __init__(self, path):
        self.taxonomy_excel = pd.read_excel(path, index_col="ACCESSION")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("[", "")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("]", "")
        self.logger = self._initLogger()

    def _initLogger(self):
        logging.basicConfig(format="%(asctime)s-%(levelname)s-%(filename)s:%(lineno)d - %(message)s", level=logging.DEBUG)
        return logging.getLogger(__name__)

    def timeit(func):
        def timeit_wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            total_time = end_time - start_time
            # self.logger.info(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
            print(f'Function {func.__name__}{args} {kwargs} Took {str(datetime.timedelta(seconds=total_time))} seconds')
            return result
        return timeit_wrapper

    def get_further_taxonomy_family(self, taxonomy_string, new_family):
        if new_family in taxonomy_string:
            return taxonomy_string
        return f"{taxonomy_string}, {new_family}"

    def get_splited_taxonomy_row(self, row):
        full_taxonomy_string = row["TAXONOMY"]
        return full_taxonomy_string.split(", ")

    def get_sorted_taxonomy_data_frame(self, column_name):
        temporary_result_dict = {}
        main_data_dict = {
            "ACCESSION": [],
            "TAXONOMY": [],
            "SURROUNDING_PAIRS": [],
            "PREVIOUS_GENE_PRODUCT_NAME": [],
            "PREVIOUS_GENE_LOCATION": [],
            "NEXT_GENE_PRODUCT_NAME": [],
            "NEXT_GENE_LOCATION": [],
            "CONTROL_REGION_LOCATION": [],
            "FROM_HOW_MANY_TAXONOMY": []
        }
        for accession, row in self.taxonomy_excel.iterrows():
            full_taxonomy_list = self.get_splited_taxonomy_row(row)
            main_family = f"{full_taxonomy_list[0]}"
            for family in full_taxonomy_list[1:]:
                main_family = self.get_further_taxonomy_family(main_family, family)
                if main_family in temporary_result_dict:
                    break
                else:
                    sorted_taxonomy_df = self.taxonomy_excel.where(self.taxonomy_excel["TAXONOMY"]\
                                                                .str.contains(main_family, regex=False)).dropna()
                    number_of_rows = sorted_taxonomy_df.shape[0]
                    unique_genes = sorted_taxonomy_df[column_name].unique()
                    if len(unique_genes) == 1:
                        temporary_result_dict[main_family] = unique_genes
                        if number_of_rows == 1:
                            main_data_dict["TAXONOMY"].append(row.TAXONOMY.replace("'Eukaryota', 'Metazoa', ", "["))
                        else:
                            main_data_dict["TAXONOMY"].append(main_family.replace("'Eukaryota', 'Metazoa', ", "["))
                        main_data_dict["PREVIOUS_GENE_PRODUCT_NAME"].append(row.PREVIOUS_GENE_PRODUCT_NAME)
                        main_data_dict["NEXT_GENE_PRODUCT_NAME"].append(row.NEXT_GENE_PRODUCT_NAME)
                        main_data_dict["NEXT_GENE_LOCATION"].append(row.NEXT_GENE_LOCATION)
                        main_data_dict["PREVIOUS_GENE_LOCATION"].append(row.PREVIOUS_GENE_LOCATION)
                        main_data_dict["ACCESSION"].append(accession)
                        main_data_dict["SURROUNDING_PAIRS"].append(row.SURROUNDING_PAIRS)
                        main_data_dict["CONTROL_REGION_LOCATION"].append(row.CONTROL_REGION_LOCATION)
                        main_data_dict["FROM_HOW_MANY_TAXONOMY"].append(number_of_rows)
                        self.logger.info(f"Row added to dictionary{accession}")
                        break
        main_df = pd.DataFrame.from_dict(main_data_dict)
        main_df = main_df.set_index("ACCESSION")
        self.logger.info(Bcolors.OKGREEN.value + f"Data Frame out of {column_name} created" + Bcolors.ENDC.value)
        return main_df

    @timeit
    def save_data_frame_to_exl_both_genes(self, path):
        previos_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("PREVIOUS_GENE_PRODUCT_NAME")
        previos_gene_df.to_excel(f"{path}/all_previous_gene_mit.xlsx")
        self.logger.info(Bcolors.OKGREEN.value + "Previos Gene saved" + Bcolors.ENDC.value)
        next_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("NEXT_GENE_PRODUCT_NAME")
        next_gene_df.to_excel(f"{path}/all_next_gene_mit.xlsx")
        self.logger.info(Bcolors.OKGREEN.value + "Next Gene saved" + Bcolors.ENDC.value)

    @timeit
    def save_data_frame_to_exl(self, path, file_name):
        surrounding_pairs_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("SURROUNDING_PAIRS")
        surrounding_pairs_gene_df.to_excel(f"{path}/{file_name}.xlsx")
        self.logger.info(Bcolors.OKGREEN.value + "SURROUNDING_PAIRS saved" + Bcolors.ENDC.value)

if __name__ == "__main__":
    sortTaxonomy = SortTaxonomy(f"/home/rszczygielski/bioinf/magisterka/geneBank/results/main_mitochondrion_no_location.xlsx")
    # sortTaxonomy.save_data_frame_to_exl_both_genes("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_results/")
    sortTaxonomy.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_results_no_location/", "surrounding_pairs_mit_no_location")
    # sortTaxonomy.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/test")

