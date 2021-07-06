import pandas as pd
import json

from sklearn.svm import SVC
import lightgbm as lgb
import joblib
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix


if __name__ == "__main__":

    # Training

    # load preprocessed MitImpact dataframe
    read_path = "./data/MitImpact_preprocessed.tsv"
    df = pd.read_table(read_path, sep="\t")

    # GenBank frequencies
    gbfreq_dict = {}
    with open('./data/gbfreq_count.json') as fp:
        gbfreq_dict = json.load(fp)

    # phylotree info. (see Files section from README for more information)
    phylotree_score_dict = {}
    with open('./data/phylotree_score.json') as fp:
        phylotree_score_dict = json.load(fp)

    # gnomAD hom./het. frequencies
    gnomad_freq_dict = {}
    with open('./data/gnomad_freq.json') as fp:
        gnomad_freq_dict = json.load(fp)

    # create default feature columns
    df.insert(loc=1, column="gb_freq", value=0.0)

    # allele frequency of homoplasmy/heteroplasmy
    df.insert(loc=2, column="AF_hom", value=0.0)
    df.insert(loc=3, column="AF_het", value=0.0)

    df.insert(loc=4, column="phylotree_score", value=0.0)

    # fill in feature values
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

    # # save dataframe
    # save_path = "./data/tmp.tsv"
    # df.to_csv(save_path, sep='\t', index=False)

    # delete position (SNV itself should not be a feature)
    df_train = df.drop(columns=["variant"])

    # Split dataset
    X_train, X_valid, Y_train, Y_valid = \
        train_test_split(
            df_train.drop(columns=["label"]), df_train["label"],
            test_size=0.1, random_state=1010101
        )

    # print(df.shape)
    # print(df.groupby("label").count())

    print(X_train.shape)
    print(Y_train.shape)

    # model = SVC(kernel='linear').fit(X_train, Y_train)
    model = lgb.LGBMClassifier().fit(X_train, Y_train)

    print("Training acc.:", model.score(X_train, Y_train))
    print("Validating acc.:", model.score(X_valid, Y_valid))

    # save model
    joblib.dump(model, 'missense_pcg_model.pkl')

    # plot confusion matrix
    plt.rcParams.update({'font.size': 16})
    plot_confusion_matrix(model, df_train.drop(columns=["label"]),
                          df_train["label"], values_format='')
    plt.savefig('confusion_matrix_missense_pcg_lgbm.png')
    plt.clf()
