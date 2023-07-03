import time
from IPython.display import display
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import pandas as pd

class SortTaxonomy():
    def __init__(self, path):
        self.taxonomy_excel = pd.read_excel(path, index_col="ACCESSION")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("[", "")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("]", "")

    def timeit(func):
        def timeit_wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = func(*args, **kwargs)
            end_time = time.perf_counter()
            total_time = end_time - start_time
            print(f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
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
        main_result_dict = {}
        main_df = pd.DataFrame(columns=['TAXONOMY', column_name, "FROM_HOW_MANY_TAXONOMY"], index=["ACCESSION"])
        for accession, row in self.taxonomy_excel.iterrows():
            full_taxonomy_list = self.get_splited_taxonomy_row(row)
            main_family = f"{full_taxonomy_list[0]}"
            for family in full_taxonomy_list[1:]:
                main_family = self.get_further_taxonomy_family(main_family, family)
                if main_family in main_result_dict:
                    break
                else:
                    sorted_taxonomy_df = self.taxonomy_excel.where(self.taxonomy_excel["TAXONOMY"]\
                                                                .str.contains(main_family, regex=False)).dropna()
                    number_of_rows = sorted_taxonomy_df.shape[0]
                    unique_number_of_genes = sorted_taxonomy_df[column_name].unique()
                    if len(unique_number_of_genes) == 1:
                        main_result_dict[main_family] = unique_number_of_genes
                        main_df.loc[accession] = [main_family, unique_number_of_genes, number_of_rows]
                        break
        main_df = main_family.set_index("ACCESSION")
        return main_df

    # @timeit
    def save_data_frame_to_exl(self, path):
        previos_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("PREVIOUS_GENE_NAME_PRODUCT")
        previos_gene_df.to_excel(f"{path}/selected_previous_gene_mit.xlsx")
        print(f"\033[92mPrevios Gene saved\033[0m")
        next_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("NEXT_GENE_NAME_PRODUCT")
        next_gene_df.to_excel(f"{path}/selected_next_gene_mit.xlsx")
        print(f"\033[92mNext Gene saved\033[0m")

if __name__ == "__main__":
    sortTaxonomy = SortTaxonomy(f"/home/rszczygielski/bioinf/magisterka/geneBank/results/main_mitochondrion.xlsx")
    sortTaxonomy.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_selected_taxonomy/not_selected_mit")
