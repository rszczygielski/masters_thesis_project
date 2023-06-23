from IPython.display import display
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"
import pandas as pd

class SortTaxonomy():
    def __init__(self, path):
        self.taxonomy_excel = pd.read_excel(path, index_col="ACCESSION")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("[", "")
        self.taxonomy_excel["TAXONOMY"] = self.taxonomy_excel["TAXONOMY"].str.replace("]", "")

    def get_further_taxonomy_family(self, taxonomy_string, new_family):
        if new_family in taxonomy_string:
            return taxonomy_string
        return f"{taxonomy_string}, {new_family}"

    def get_splited_taxonomy_row(self, row):
        full_taxonomy_string = row["TAXONOMY"]
        return full_taxonomy_string.split(", ")

    def get_sorted_taxonomy_data_frame(self, column_name):
        main_result_dict = {}
        main_df = pd.DataFrame(columns=['TAXONOMY', column_name], index=["ACCESSION"])
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
                    unique_number_previus_genes = sorted_taxonomy_df[column_name].unique()
                    if len(unique_number_previus_genes) == 1:
                        print(main_family, "taxonomy added")
                        main_result_dict[main_family] = unique_number_previus_genes
                        main_df.loc[accession] = [main_family, unique_number_previus_genes]
                        break
        return main_df

    def save_data_frame_to_exl(self, path, type_of_file):
        previos_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("PREVIOUS_GENE")
        previos_gene_df.to_excel(f"{path}/selected_previous_gene_mit{type_of_file}.xlsx")
        print(f"\033[92mPrevios Gene saved\033[0m")
        next_gene_df = sortTaxonomy.get_sorted_taxonomy_data_frame("NEXT_GENE")
        next_gene_df.to_excel(f"{path}/selected_next_gene_mit{type_of_file}.xlsx")
        print(f"\033[92mNext Gene saved\033[0m")

if __name__ == "__main__":
    type_of_file = "2"
    sortTaxonomy = SortTaxonomy(f"/home/rszczygielski/bioinf/magisterka/geneBank/results/General_info_mitochondrion_{type_of_file}.xlsx")
    sortTaxonomy.save_data_frame_to_exl(f"/home/rszczygielski/bioinf/magisterka/geneBank/results/main_selected_taxonomy/mitochondrion_{type_of_file}/", type_of_file)
