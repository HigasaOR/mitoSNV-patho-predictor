import re
import json
import pandas as pd


def extract_mito_snp_from_clinvar():

    read_path = "./data/variant_summary.txt"
    output_path = "./data/ClinVar_mito_snp_sig.tsv"

    df = pd.read_table(read_path, sep="\t", dtype=str)

    mt_snp_df = df.loc[
        (df['Chromosome'] == "MT")
        & (df["Type"] == "single nucleotide variant")
        & (df["ClinicalSignificance"] != "not provided")
        & (df["ClinicalSignificance"] != "Uncertain significance")
        & (df["ClinicalSignificance"] != "Conflicting interpretations of pathogenicity")
        & (df["ClinicalSignificance"] != "drug response")
    ]

    mt_snp_df.to_csv(output_path, sep='\t', index=False)


def find_str_var_in_name(str_name):
    found = re.search(r"[0-9]+[ATCG]>[ATCG]", str_name)
    if found is None:
        raise Exception("No variant string information in Name.")
    return found.group()


def make_label():

    read_path = "./data/ClinVar_mito_snp_sig.tsv"

    df = pd.read_table(read_path, sep="\t", dtype=str)

    str_patho = "Pathogenic"
    # str_likely_patho = "Likely pathogenic"
    str_patho_li_patho = "Pathogenic/Likely pathogenic"
    str_benign = "Benign"
    str_likely_benign = "Likely benign"

    cal_dict = {}

    for i, row in df.iterrows():

        found = re.search(r"[0-9]+[ATCG]>[ATCG]", row["Name"])
        if found is None:
            continue

        str_var = row["PositionVCF"].zfill(5) + row["ReferenceAlleleVCF"] \
            + ">" + row["AlternateAlleleVCF"]

        str_clin_sig = row["ClinicalSignificance"]

        if str_var in cal_dict:
            continue

        if str_clin_sig == str_patho or str_clin_sig == str_patho_li_patho:
            cal_dict[str_var] = 1
        elif str_clin_sig == str_benign or str_clin_sig == str_likely_benign:
            cal_dict[str_var] = 0

    with open('./data/label_clinvar.json', 'w') as fp:
        json.dump(cal_dict, fp, sort_keys=True, indent=4)


def add_emory_label():
    read_path = "./data/EmoryMito.tsv"
    df = pd.read_table(read_path, sep="\t", dtype=str)

    label_dict = {}
    with open('./data/label_clinvar.json') as fp:
        label_dict = json.load(fp)

    str_patho = "Pathogenic"
    str_benign = "Benign"
    str_likely_benign = "Likely benign"
    str_vous = "VOUS"

    new_add_count = 0

    for _, row in df.iterrows():

        str_emv = row["EmVClass"]
        if str_emv == str_vous:
            continue

        s = row["Transcript, Nucleotide"]
        found = re.search(r"[0-9]+[ATCG]>[ATCG]", s)
        if found is None:
            continue
        str_tmp = found.group()
        str_var = str_tmp[:-3].zfill(5) + str_tmp[-3:]

        label_emory = None
        if str_emv == str_benign or str_emv == str_likely_benign:
            label_emory = 0
        elif str_emv == str_patho:
            label_emory = 1

        label_clinvar = None

        if str_var in label_dict:
            label_clinvar = label_dict[str_var]
        else:
            new_add_count += 1
            label_dict[str_var] = label_emory

        if (label_clinvar is not None) and \
                (label_clinvar != label_emory):
            print("Inconsistancy:", str_var)
            print("CLINVAR annotation is chosen.")

    with open('./data/label_clinvar_emory.json', 'w') as fp:
        json.dump(label_dict, fp, sort_keys=True, indent=4)


def generate_str_var(pos, ref, alt):
    return pos.zfill(5) + ref + ">" + alt


if __name__ == "__main__":

    # Extract only mito variants data from clinvar
    # that has clinical significance.
    extract_mito_snp_from_clinvar()

    # Make label json file.
    make_label()

    # add label from EmoryMito
    # add_emory_label()
