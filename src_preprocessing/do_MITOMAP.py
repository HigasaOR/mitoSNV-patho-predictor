import re
import json
import pandas as pd


def reorder_disease():
    read_path = "./download/disease.cgi.tsv"
    dis_df = pd.read_table(read_path, sep="\t", dtype=str)

    # print(dis_df)

    dis_df.sort_values(by=['pos'])
    dis_df.drop(columns=['id'], inplace=True)

    save_path = "./download/disease.tsv"
    dis_df.to_csv(save_path, sep='\t', index=False)


def reorder_poly():
    read_path = "./download/polymorphisms.cgi.tsv"
    poly_df = pd.read_table(read_path, sep="\t", dtype=str)

    # print(dis_df)

    poly_df.sort_values(by=['pos'])
    poly_df.drop(columns=['id'], inplace=True)

    save_path = "./download/polymorphisms.tsv"
    poly_df.to_csv(save_path, sep='\t', index=False)


def count_gb_freq():
    read_path = "./download/polymorphisms.tsv"
    poly_df = pd.read_table(read_path, sep="\t", dtype=str)

    freq_dict = {}

    reg = re.compile(r'^[A-Z]$')

    for _, row in poly_df.iterrows():

        if (not reg.match(row["ref"])) or (not reg.match(row["alt"])):
            continue

        str_var = row["pos"].zfill(5) + row["ref"] \
            + ">" + row["alt"]

        if str_var not in freq_dict:
            freq_dict[str_var] = float(row["gbfreq"])

    with open('gbfreq_count.json', 'w') as fp:
        json.dump(freq_dict, fp, sort_keys=True, indent=4)


def extract_mitotip_score():
    read_path = "./download/mitotip_scores.tsv"
    mitip_df = pd.read_table(read_path, sep="\t", dtype=str)

    score_dict = {}

    reg = re.compile(r'^[A-Z]$')

    for _, row in mitip_df.iterrows():

        if not reg.match(row["Alt"]):
            str_var = row["Position"].zfill(5) + row["rCRS"] \
                + ">" + row["rCRS"]
        else:
            str_var = row["Position"].zfill(5) + row["rCRS"] \
                + ">" + row["Alt"]

        if str_var not in score_dict:
            score_dict[str_var] = float(row["MitoTIP_Score"])

    with open('./data/mitoTIP_score.json', 'w') as fp:
        json.dump(score_dict, fp, sort_keys=True, indent=4)

    # N = 16100
    # for i in range(N):
    #     print(i)


if __name__ == "__main__":

    # reordering
    # reorder_disease()
    reorder_poly()

    # extract GenBank frequencies
    count_gb_freq()

    # extract MitoTIP scores
    extract_mitotip_score()
