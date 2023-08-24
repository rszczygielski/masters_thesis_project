import pandas as pd

class Analysis():
    def __init__(self, excel_main_file_path) -> None:
        self.mian_data_frame = pd.read_excel(excel_main_file_path)
        self.description_table = pd.DataFrame(columns=['DESCRIPTION', 'TAXONOMY_GROUP'])

    def get_grouped_df(self, taxonomy_group_list):
        columns = ["SURROUNDING_PAIRS", "REPETITIONS_NUMBER", "TAXONOMY_GROUP"]
        main_df = pd.DataFrame(columns=columns)
        for taxonomy_group in taxonomy_group_list:
            selecte_group_df = self.mian_data_frame.where(self.mian_data_frame["TAXONOMY"].str.contains(taxonomy_group)).dropna()
            grouped_df = selecte_group_df.groupby('SURROUNDING_PAIRS').size().reset_index(name='REPETITIONS_NUMBER')
            grouped_df["TAXONOMY_GROUP"] = taxonomy_group
            # print(grouped_df)
            main_df = main_df.merge(grouped_df, "outer")
        main_df = main_df.sort_values(by=['TAXONOMY_GROUP', 'REPETITIONS_NUMBER'], ascending=[True, False])
        return main_df

    def add_description(self, df, description_dict):
        for group, description in description_dict.items():
            df.loc[df['TAXONOMY_GROUP'] == group, 'DESCRIPTION'] = description
        return df

    def save_data_frame_to_exl(self, path, file_name, data_frame):
        data_frame.to_excel(f"{path}/{file_name}.xlsx")

    def ad_description_to_description_table(self, df):
        taxnomy_description = df[['DESCRIPTION', 'TAXONOMY_GROUP']].drop_duplicates(subset=['TAXONOMY_GROUP'])
        self.description_table = self.description_table.merge(taxnomy_description, how="outer")


if __name__ == "__main__":
    analisis_instance = Analysis("/home/rszczygielski/bioinf/magisterka/geneBank/results/main_mitochondrion_no_location.xlsx")

    # MAMMALS
    mammals_df = analisis_instance.get_grouped_df(['Rodentia','Chiroptera','Soricomorpha','Primates','Carnivora','Artiodactyla','Diprotodontia',
                                                'Lagomorpha','Didelphimorphia','Cetacea','Dasyuromorphia','Afrosoricida','Erinaceomorpha',
                                                'Cingulata','Peramelemorphia','Scandentia','Perissodactyla','Macroscelidea','Pilosa',
                                                'Monotremata','Sirenia','Proboscidea'])
    mammals_df = analisis_instance.add_description(mammals_df, {
    'Rodentia': 'Gryzonie to grupa ssaków charakteryzujących się obecnością siekaczy i zębów trzonowych o stałym wzroście.',
    'Chiroptera': 'Nietoperze to ssaki latające, posiadające skrzydła zbudowane z błony skórnej rozpiętej między kośćmi.',
    'Soricomorpha': 'Ssaki owadożerne to małe zwierzęta o specjalnie przystosowanych zębach do zjedzenia owadów i bezkręgowców.',
    'Primates': 'Naczelne to grupa ssaków o rozwiniętych mózgach, zróżnicowanym społeczeństwie i zdolnościach manualnych.',
    'Carnivora': 'Drapieżniki to ssaki mięsożerne, które posiadają przystosowania do polowań i konsumpcji mięsa.',
    'Artiodactyla': 'Parzystokopytne to ssaki o parzystej liczbie palców, często hodowane ze względu na mięso lub produkty uboczne.',
    'Diprotodontia': 'Diprotodonta to rzęd ssaków australijskich, charakteryzujący się wydłużonymi szczękami i trzonowcami.',
    'Lagomorpha': 'Zające i świstaki to ssaki charakteryzujące się dużymi tylnymi nogami i długimi uszami.',
    'Didelphimorphia': 'Oposy to amerykańskie ssaki torbacze, posiadające charakterystyczną torbę na brzuchu.',
    'Cetacea': 'Walenie to grupa ssaków wodnych, które zaliczają się zarówno do delfinów, jak i waleni.',
    'Dasyuromorphia': 'Tasmania to ssaki drapieżne związane głównie z Australią i Tasmanią.',
    'Afrosoricida': 'Afrosory to grupa ssaków żyjących w Afryce, często przypominających wyglądem owadożerne.',
    'Erinaceomorpha': 'Jeżowce to niewielkie ssaki, często występujące w Europie, charakteryzujące się kolczastymi grzbietami.',
    'Cingulata': 'Panzery to ssaki z rzędu pancerników, mające charakterystyczne pancerze na grzbiecie.',
    'Peramelemorphia': 'Bandicoots to australijskie ssaki, przypominające wyglądem króliki lub szczury.',
    'Scandentia': 'Wiewiórkowce to małe, nadrzewne ssaki o długim ogonie i małych oczach.',
    'Perissodactyla': 'Nieparzystokopytne to grupa ssaków o nieparzystej liczbie palców na nogach, w tym konie i nosorożce.',
    'Macroscelidea': 'Erytrogliksy to ssaki nadrzewne, przypominające wyglądem małe gryzonie.',
    'Pilosa': 'Leniwece to ssaki z Ameryki Południowej, charakteryzujące się wolnym tempem życia i niską przemianą materii.',
    'Monotremata': 'Jajorodne to ssaki, które składają jaja zamiast rodzić żywe młode, m.in. dziobak i kolczatka.',
    'Sirenia': 'Sirenie to morskie ssaki o dużych ciałach, w tym manaty i dugongi.',
    'Proboscidea': 'Słoniowate to ogromne lądowe ssaki, wyróżniające się długim trąbą i dużymi kłami.',
    })
    analisis_instance.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/data_analysis/main_tables", "mammals_table", mammals_df)
    analisis_instance.ad_description_to_description_table(mammals_df)

    # REPTILIA
    reptilia_df = analisis_instance.get_grouped_df(['Varanus', 'Python', 'Crocodilus', 'Naja', 'Gecko', 'Agkistrodon', 'Alligator', 'Vipera', 'Iguana', 'Chamaeleo'])
    reptilia_df = analisis_instance.add_description(reptilia_df, {
    'Varanus': 'Rodzaj gadów z rodziny waranów.',
    'Python': 'Rodzaj węża z rodziny pytonów.',
    'Crocodilus': 'Rodzaj krokodyla z rodziny krokodyli.',
    'Naja': 'Rodzaj jadowitego węża z rodziny zdradnicowatych.',
    'Gecko': 'Rodzaj jaszczurki z rodziny gekonowatych.',
    'Agkistrodon': 'Rodzaj węża z rodziny mocno jadowitych żmijowatych.',
    'Alligator': 'Rodzaj drapieżnego gadziny z rodziny aligatorowatych.',
    'Vipera': 'Rodzaj jadowitego węża z rodziny żmijowatych.',
    'Iguana': 'Rodzaj jaszczurki z rodziny legwanów.',
    'Chamaeleo': 'Rodzaj jaszczurki z rodziny kameleonów.'})
    analisis_instance.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/data_analysis/main_tables", "reptilia_table", reptilia_df)
    analisis_instance.ad_description_to_description_table(reptilia_df)

    # AVES
    aves_df = analisis_instance.get_grouped_df(['Passer', 'Corvus', 'Turdus', 'Struthio', 'Aquila', 'Falcon', 'Anas', 'Gallus', 'Columba', 'Accipiter'])
    aves_df = analisis_instance.add_description(aves_df, {
    'Passer': 'Rodzaj ptaka z rodziny wróbli.',
    'Corvus': 'Rodzaj ptaka z rodziny krukowatych.',
    'Turdus': 'Rodzaj ptaka z rodziny drozdowatych.',
    'Struthio': 'Rodzaj ptaka z rodziny strusi.',
    'Aquila': 'Rodzaj ptaka drapieżnego z rodziny jastrzębiowatych.',
    'Falcon': 'Rodzaj ptaka drapieżnego z rodziny sokołowatych.',
    'Anas': 'Rodzaj ptaka z rodziny kaczkowatych.',
    'Gallus': 'Rodzaj ptaka z rodziny kurowatych, obejmujący kury domowe.',
    'Columba': 'Rodzaj ptaka z rodziny gołębiowatych.',
    'Accipiter': 'Rodzaj ptaka drapieżnego z rodziny jastrzębiowatych.'})
    analisis_instance.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/data_analysis/main_tables", "aves_table", aves_df)
    analisis_instance.ad_description_to_description_table(aves_df)

    # amphibiorum AMPHIBIORUM
    amphibiorum_df = analisis_instance.get_grouped_df( ['Anura', 'Caudata', 'Gymnophiona'])
    amphibiorum_df = analisis_instance.add_description(amphibiorum_df, {
        'Anura': 'Rząd płazów bezogonowych, obejmujący m.in. żaby i ropuchy.',
        'Caudata': 'Rząd płazów ogoniastych, obejmujący m.in. salamandry i młoteczki.',
        'Gymnophiona': 'Rząd płazów beznogich, obejmujący m.in. węże krajowe i glebne.' })
    analisis_instance.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/data_analysis/main_tables", "amphibiorum_table", amphibiorum_df)
    analisis_instance.ad_description_to_description_table(amphibiorum_df)

    analisis_instance.save_data_frame_to_exl("/home/rszczygielski/bioinf/magisterka/geneBank/data_analysis/main_tables","description_table",
                                             analisis_instance.description_table)


