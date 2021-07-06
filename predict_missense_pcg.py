import pandas as pd
import json

from sklearn.svm import SVC
import lightgbm as lgb
import joblib
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix


if __name__ == "__main__":

    read_path = "./data/MitImpact_all_preprocessed.tsv"
    df = pd.read_table(read_path, sep="\t")

    gbfreq_dict = {}
    with open('./data/gbfreq_count.json') as fp:
        gbfreq_dict = json.load(fp)

    phylotree_score_dict = {}
    with open('./data/phylotree_score.json') as fp:
        phylotree_score_dict = json.load(fp)

    gnomad_freq_dict = {}
    with open('./data/gnomad_freq.json') as fp:
        gnomad_freq_dict = json.load(fp)

    df.insert(loc=1, column="gb_freq", value=0.0)

    df.insert(loc=2, column="AF_hom", value=0.0)
    df.insert(loc=3, column="AF_het", value=0.0)

    df.insert(loc=4, column="phylotree_score", value=0.0)

    print("adding scores...")
    for index, row in df.iterrows():
        str_var = row["variant"]

        if str_var in gbfreq_dict:
            df.loc[index, "gb_freq"] = gbfreq_dict[str_var]

        if str_var in gnomad_freq_dict:
            df.loc[index, "AF_hom"] = gnomad_freq_dict[str_var][0]
            df.loc[index, "AF_het"] = gnomad_freq_dict[str_var][1]

        if str_var in phylotree_score_dict:
            df.loc[index, "phylotree_score"] = phylotree_score_dict[str_var]
        else:
            df.loc[index, "phylotree_score"] = 25.0

        if index % 100 == 0:
            print(f"row {index} done")

    # save file since the above really const much time
    # maybe wont need to to again, ust load the dataframe from here
    save_path = "./data/MitoImpact_all_p_added.tsv"
    df.to_csv(save_path, sep='\t', index=False)

    # read_path = "./data_1/MitoImpact_all_added.tsv"
    # df = pd.read_table(read_path, sep="\t")

    df_feature = df.drop(columns=["variant"])

    print(df_feature.shape)

    # load model
    model = joblib.load('missense_pcg_model.pkl')

    # predict
    proba_list = model.predict_proba(df_feature)

    # add result to dataframe
    df_pred = pd.DataFrame(proba_list, columns=[
                           'prob_benign', 'prob_pathogenic'])

    df_j = pd.concat([df, df_pred], axis=1)

    # save dataframe
    save_path = "missense_pcg_pred.tsv"
    df_j.to_csv(save_path, sep='\t', index=False)
