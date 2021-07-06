import json
import pandas as pd
from collections import Counter
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
import lightgbm as lgb
import joblib

import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix


def find_gene_function_type(pos):
    ''' This function finds the gene coding part of the given position.
        It is not well-designed with the perspective of performance 
        as it opens a file to search every time.
    '''

    df = pd.read_csv("./data/Mit_gene_function_list.csv")

    gt = None  # gene type

    # iterate through the data to find out which region the given position
    # is included
    for index, row in df.iterrows():
        start_pos, end_pos = row['Starting'], row['Ending']
        if start_pos <= pos and pos <= end_pos:
            gt = row['Gene type']
            break

    return gt


def extract_tRNA_label():
    ''' Extracts pathogenicity labels of tRNA coding region from ClinVar
        label data.
    '''

    label_dict = {}
    with open('./data/label_clinvar.json') as fp:
        label_dict = json.load(fp)

    label_tRNA_dict = {}

    # if the label is of tRNA coding region, add to the extracted dict.
    for var in label_dict:
        pos = int(var[:5])
        str_type = find_gene_function_type(pos)
        if str_type == "tRNA":
            label_tRNA_dict[var] = label_dict[var]

    # dump the extracted dict. to a file
    with open('./data/label_clinvar_tRNA.json', 'w') as fp:
        json.dump(label_tRNA_dict, fp, sort_keys=True, indent=4)


if __name__ == "__main__":

    print("Start processing...")

    # Extract tRNA coding gene label.
    # The default extracted label file is stored at
    # data/label_clinvar_tRNA.json
    extract_tRNA_label()

    # load label file to a dict.
    label_tRNA_dict = {}
    with open('./data/label_clinvar_tRNA.json') as fp:
        label_tRNA_dict = json.load(fp)

    # see #labels
    # print(Counter(label_tRNA_dict.values()))

    # make the dict. a pd dataframe
    df = pd.DataFrame(list(label_tRNA_dict.items()),
                      columns=["variant", "label"])

    # load features

    # MitoTIP score
    mitip_score_dict = {}
    with open('./data/mitoTIP_score.json') as fp:
        mitip_score_dict = json.load(fp)

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
    df.insert(loc=1, column="mitoTIP_score", value=12.0)
    df.insert(loc=2, column="gb_freq", value=0.0)

    df.insert(loc=3, column="AF_hom", value=0.0)
    df.insert(loc=4, column="AF_het", value=0.0)

    df.insert(loc=5, column="phylotree_score", value=0.0)

    # fill in feature values
    for index, row in df.iterrows():
        str_var = row["variant"]
        if str_var in mitip_score_dict:
            df.loc[index, "mitoTIP_score"] = mitip_score_dict[str_var]
        if str_var in gbfreq_dict:
            df.loc[index, "gb_freq"] = gbfreq_dict[str_var]

        if str_var in gnomad_freq_dict:
            df.loc[index, "AF_hom"] = gnomad_freq_dict[str_var][0]
            df.loc[index, "AF_het"] = gnomad_freq_dict[str_var][1]

        if str_var in phylotree_score_dict:
            df.loc[index, "phylotree_score"] = phylotree_score_dict[str_var]
        else:
            df.loc[index, "phylotree_score"] = 25.0

    # Now the dataframe has features and label included.
    # (rows: SNVs; columns: features + label)

    # Training
    print("Start training...")

    # delete position (SNV itself(ex. 3412G>A) should not be a feature)
    df_train = df.drop(columns=["variant"])

    # Split dataset
    X_train, X_valid, Y_train, Y_valid = \
        train_test_split(
            df_train.drop(columns=["label"]), df_train["label"],
            test_size=0.1, random_state=1010101
        )

    print("#Training data:", X_train.shape[0])
    print("#Validating data:", X_valid.shape[0])

    # using LGBMClassifier as the machine learning model
    model = lgb.LGBMClassifier().fit(X_train, Y_train)
    # model = SVC(kernel='linear').fit(X_train, Y_train)

    # save model
    model_path = 'tRNA_model.pkl'
    joblib.dump(model, model_path)
    print(f"Model saved at {model_path}")

    print("Training acc.:", model.score(X_train, Y_train))
    print("Validating acc.:", model.score(X_valid, Y_valid))

    # Plot confusion matrix
    plt.rcParams.update({'font.size': 16})
    plot_confusion_matrix(model, X_valid, Y_valid, values_format='')
    plt.savefig('confusion_matrix_tRNA_lgbm_valid.png')
    plt.clf()
