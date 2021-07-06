import pandas as pd
import json

from sklearn.svm import SVC
import lightgbm as lgb
import joblib
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix


def extract_labeled_from_MitImpact():

    # Extract MitImpact data with ClinVar significance.

    read_path = "../download/MitImpact_db_3.0.6.tsv"
    df = pd.read_table(read_path, sep="\t", dtype=str)

    df.insert(loc=0, column="variant", value="")

    print("add variant string...")
    for index, row in df.iterrows():
        str_pos = row["Start"].zfill(5)
        str_ref, str_alt = row["Ref"], row["Alt"]

        str_var = str_pos + str_ref + ">" + str_alt
        df.loc[index, "variant"] = str_var

        if index % 100 == 0:
            print(f"row {index} done")

    print("saving file...")
    save_path = "../data/MitImpact_all.tsv"
    df.to_csv(save_path, sep='\t', index=False)
    print("done")

# def delete_columns_MitImpact():


def if_value_missing(s):
    if s == ".":
        return "1"
    return "0"


def change_to_zero_if_missing(s):
    if s == ".":
        return "0"
    else:
        return s


def process_dataframe():

    read_path = "../data/MitImpact.tsv"
    df = pd.read_table(read_path, sep="\t", dtype=str)

    col_dict = {}
    with open('../data/col_MitImpact.json') as fp:
        col_dict = json.load(fp)

    cols_del = col_dict["columns_del"] + \
        col_dict["columns_del_2nd"] + col_dict["columns_del_3rd"]

    print("dropping columns...")
    df.drop(columns=cols_del, inplace=True)

    cols_loss = col_dict["columns_loss"]
    cols_one_hot = col_dict["columns_one_hot"]

    df = pd.get_dummies(df, columns=cols_one_hot)

    print("padding...")
    for col in cols_loss:
        print(col)
        df[col] = df[col].apply(change_to_zero_if_missing)

    # only one column cannot be processed;
    # its missing too many values so I am going to just delete the row.
    # df = df[df["PhyloP_100V"] != '.']

    print("saving file...")
    save_path = "../data/MitoImpact_all_preprocessed.tsv"
    df.to_csv(save_path, sep='\t', index=False)
    print("done")


if __name__ == "__main__":

    # extract_labeled_from_MitImpact()
    process_dataframe()
