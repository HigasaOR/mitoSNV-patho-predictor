import re
import json
import pandas as pd


def extract_gnomad_freq():
    read_path = "../download/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv"
    gnomad_df = pd.read_table(read_path, sep="\t", dtype=str)

    freq_dict = {}
    reg = re.compile(r'^[A-Z]$')

    for _, row in gnomad_df.iterrows():
        if reg.match(row["ref"]) and reg.match(row["alt"]):
            str_var = row["position"].zfill(5) + row["ref"] \
                + ">" + row["alt"]
            if str_var not in freq_dict:
                freq_dict[str_var] = [
                    float(row["AF_hom"]), float(row["AF_het"])]

    with open('../data/gnomad_freq.json', 'w') as fp:
        json.dump(freq_dict, fp, sort_keys=True, indent=4)


if __name__ == "__main__":

    extract_gnomad_freq()
