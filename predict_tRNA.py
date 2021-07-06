import pandas as pd
import json
import lightgbm as lgb
import joblib


if __name__ == "__main__":

    df = pd.read_csv("./data/tRNA_var_list.tsv", sep='\t')

    mitip_score_dict = {}
    with open('./data_0/mitoTIP_score.json') as fp:
        mitip_score_dict = json.load(fp)

    gbfreq_dict = {}
    with open('./data_0/gbfreq_count.json') as fp:
        gbfreq_dict = json.load(fp)

    phylotree_score_dict = {}
    with open('./data_0/phylotree_score.json') as fp:
        phylotree_score_dict = json.load(fp)

    gnomad_freq_dict = {}
    with open('./data_0/gnomad_freq.json') as fp:
        gnomad_freq_dict = json.load(fp)

    df.insert(loc=3, column="mitoTIP_score", value=12.0)
    df.insert(loc=4, column="gb_freq", value=0.0)

    df.insert(loc=5, column="AF_hom", value=0.0)
    df.insert(loc=6, column="AF_het", value=0.0)

    df.insert(loc=7, column="phylotree_score", value=0.0)

    for index, row in df.iterrows():
        str_var = str(row["pos"]).zfill(5) + row["ref"] + ">" + row["alt"]
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

    df_feature = df.drop(columns=["pos", "ref", "alt"])

    model = joblib.load('tRNA_model_2.pkl')
    proba_list = model.predict_proba(df_feature)
    df_pred = pd.DataFrame(proba_list, columns=[
                           'prob_benign', 'prob_pathogenic'])

    df_j = pd.concat([df, df_pred], axis=1)
    save_path = "tRNA_pred.tsv"
    df_j.to_csv(save_path, sep='\t', index=False)
