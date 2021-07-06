import itertools
import pandas as pd


def get_pos_list():

    pos_list = []
    df = pd.read_csv("../data/Mit_gene_function_list.csv")

    for index, row in df.iterrows():
        if not isinstance(row['Gene type'], str) or ("tRNA" not in row['Gene type']):
            continue

        start_pos, end_pos = int(row['Starting']), int(row['Ending'])
        pos_list += [i for i in range(start_pos, end_pos+1)]

    return pos_list


if __name__ == "__main__":

    char_list = ['A', 'T', 'C', 'G']

    pos_list = get_pos_list()
    to_predict_set = [
        pos_list,
        char_list,
        char_list
    ]

    var_list = [i for i in itertools.product(*to_predict_set) if i[1] != i[2]]

    df = pd.DataFrame(var_list, columns=['pos', 'ref', 'alt'])
    save_path = "../data/tRNA_var_list.tsv"
    df.to_csv(save_path, sep='\t', index=False)
