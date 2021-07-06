
import json


def make_new_str(str_var):
    new_str = str_var[1:-1].zfill(5) + str_var[0] + '>' + str_var[-1]
    return new_str


def phylo_tree_score(str_lvls):

    # print(str_lvls)

    tot_count = 0
    tot_level = 0

    lvls_list = str_lvls.split(',')
    for str_lvl in lvls_list:
        lvl_pair = str_lvl.strip().split('(')
        level = int(lvl_pair[0])
        count = int(lvl_pair[1][:-1])
        tot_level += level * count
        tot_count += count

    return float(tot_level)/tot_count


if __name__ == "__main__":

    phylo_dict = {}
    with open('./phylotree_dict.json') as fp:
        phylo_dict = json.load(fp)

    make_new_str("A10005G")

    new_dict = {}
    score_dict = {}
    for str_var in phylo_dict:
        new_str = make_new_str(str_var)
        new_dict[new_str] = phylo_dict[str_var]
        score_dict[new_str] = phylo_tree_score(phylo_dict[str_var])

    with open('./data/new_phylotree_dict.json', 'w') as fp:
        json.dump(new_dict, fp, sort_keys=True, indent=4)

    with open('./data/phylotree_score.json', 'w') as fp:
        json.dump(score_dict, fp, sort_keys=True, indent=4)
