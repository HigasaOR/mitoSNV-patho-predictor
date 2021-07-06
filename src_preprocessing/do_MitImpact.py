import pandas as pd
import json


def extract_labeled_from_MitImpact():
    ''' Extract MitImpact data with ClinVar significance.
    '''

    # load ClinVar label
    label_dict = {}
    with open("../data/label_clinvar.json") as fp:
        label_dict = json.load(fp)

    # read MitImpact database file
    read_path = "../download/MitImpact_db_3.0.6.tsv"
    df = pd.read_table(read_path, sep="\t", dtype=str)

    df["label"] = 0  # Add label meanwhile.
    df.insert(loc=1, column="variant", value="")

    loc_list = []

    for index, row in df.iterrows():
        str_pos = row["Start"].zfill(5)
        str_ref, str_alt = row["Ref"], row["Alt"]

        str_var = str_pos + str_ref + ">" + str_alt

        if str_var in label_dict:
            loc_list.append(index)
            df.loc[index, "label"] = label_dict[str_var]
            df.loc[index, "variant"] = str_var

    df_labeled = df.iloc[loc_list]  # .sort_values("Start")

    save_path = "../data/MitImpact_ClinVar.tsv"
    df_labeled.to_csv(save_path, sep='\t', index=False)


def if_value_missing(s):
    ''' retunr 1 (means "is missing") if s(the value of cell)
        is missing(a '.' in the database file); 0 otherwise
    '''

    if s == ".":
        return "1"
    return "0"


def change_to_zero_if_missing(s):
    ''' return 0 if the value is missing for padding
    '''

    if s == ".":
        return "0"
    else:
        return s


def process_dataframe():
    ''' process the MitImpact dataframe (with ClinVar significance)
    '''

    # read dataframe
    read_path = "../data/MitoImpact_ClinVar.tsv"
    df = pd.read_table(read_path, sep="\t", dtype=str)

    # load the json file specifying how to deal with the columns
    col_dict = {}
    with open('./data/col_MitImpact.json') as fp:
        col_dict = json.load(fp)

    # which columns to delete
    cols_del = col_dict["columns_del"] + \
        col_dict["columns_del_2nd"] + col_dict["columns_del_3rd"]

    # delete the columns
    df.drop(columns=cols_del, inplace=True)

    # which columns have missing values and which columns need to be
    # transformed to one-hot values (e.g., strings)
    cols_loss = col_dict["columns_loss"]
    cols_one_hot = col_dict["columns_one_hot"]

    # one-hot encoding
    df = pd.get_dummies(df, columns=cols_one_hot)

    # pad zero if value missing
    for col in cols_loss:
        # str_col_miss = col + "_missing"
        # insert_col = df[col].apply(if_value_missing)
        # insert_pos = df.columns.get_loc(col) + 1
        # df.insert(insert_pos, str_col_miss, insert_col)
        df[col] = df[col].apply(change_to_zero_if_missing)

    df = df[[col for col in df.columns if col != "label"]+["label"]]

    save_path = "../data/MitoImpact_preprocessed.tsv"
    df.to_csv(save_path, sep='\t', index=False)


if __name__ == "__main__":

    # preprocessing
    extract_labeled_from_MitImpact()
    process_dataframe()
