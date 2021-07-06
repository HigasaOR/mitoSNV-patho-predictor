import re
from collections import Counter
import json
from bs4 import BeautifulSoup


def prep_phylo_occ_json(file_path):
    ''' This function generates a json file.
        The json file contains information about the levels of
        occurence of variants in phylotree. The phylotree is
        stored in a html file.
    '''

    # prepare regex
    pattern = re.compile(r"[AaTtCcGg]+[0-9]+[AaTtCcGg]+")

    webpage_html = ""
    with open(file_path, 'r', encoding='utf-8') as f:
        webpage_html = f.read()

    soup = BeautifulSoup(webpage_html, 'html.parser')

    table = soup.find("table")
    tree_table_rows = table.find_all("tr", {"class": "xl10617826"})

    occ_dict = {}

    for row in tree_table_rows:

        cells = row.find_all("td")

        # # of precedent spaces (indicating which level of tree)
        num_pre_space = 0
        sp_flag = True

        for cell in cells:

            text_list_in_cell = cell.find_all(text=True)
            text_in_cell = ""
            for text in text_list_in_cell:
                if text.strip() == "":
                    continue
                text_in_cell += " " + text.strip()
            new_text_list_in_cell = text_in_cell.strip().split()
            text_in_cell = " ".join(new_text_list_in_cell)

            if text_in_cell != "":
                # set flag to False if no more precedent spaces
                sp_flag = False

                for text in new_text_list_in_cell:
                    res = pattern.match(text.strip())
                    if res is None:
                        nothing_to_do_for_now = True
                    else:
                        variant = res.group(0).upper()
                        # print("var : ", variant)
                        if variant in occ_dict:
                            pos_list = occ_dict[variant]
                            pos_list.append(num_pre_space + 1)
                            occ_dict[variant] = pos_list
                        else:
                            occ_dict[variant] = [num_pre_space + 1]

            # count number of precedent spaces
            if sp_flag == True:
                num_pre_space += 1

    occ_str_dict = {}  # dict of customized occurence string

    for variant, pos_list in occ_dict.items():

        # count position occurence
        pos_counter = Counter(pos_list)
        occ_string = ""  # occurence string

        for pos, pcnt in pos_counter.items():
            occ_string += str(pos) + "(" + str(pcnt) + "), "

        occ_string = occ_string[:-2]
        occ_str_dict[variant] = occ_string

    # print(occ_str_dict)

    with open('./data/phylotree_dict.json', 'w') as f:
        json.dump(occ_str_dict, f, sort_keys=True, indent=4)


if __name__ == "__main__":

    # file_path = "./mtDNA_tree_Build_17.html"
    file_path = "./downloads/mtDNA.html"  # into utf-8

    prep_phylo_occ_json(file_path)
