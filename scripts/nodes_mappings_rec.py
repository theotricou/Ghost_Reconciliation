#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr
import sys, os
import glob
import re
import pandas as pd


ale = tr("aletree", format = 1)
ext = tr("SampledSpeciesTree.nwk", format = 1)

map = dict()

map[ale.name] = ext.name

for i in ext:
    ext_node = i
    ale_node = ale.search_nodes(name=i.name)[0]
    while ale_node.name not in map.keys() :
        map[ale_node.name] = ext_node.name
        ext_node = ext_node.up
        ale_node = ale_node.up

all_uTs = []
uTs = glob.glob("*sampledtree.nwk.ale.uTs")
for i in uTs:
    if os.stat(i).st_size != 0:
        all_uTs.append(pd.read_csv(i, sep="\t", header = None))

full_uTs = pd.concat(all_uTs, ignore_index=True)
full_uTs[[0, 1]] = full_uTs[[0, 1]].astype(str)


def is_to_strict_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0]
    if n_to in n_from.get_descendants():
        return(True)
    else:
        return(False)

def is_to_time_descendant(name_from ,name_to, tree):
    n_from = tree.search_nodes(name = name_from)[0]
    n_to = tree.search_nodes(name = name_to)[0].up
    if n_from.get_distance(tree) >= n_to.get_distance(tree):
        return(False)
    elif n_from.get_distance(tree) <= n_to.get_distance(tree):
        return(True)
    else:
        return("ERROR")


res_df = pd.DataFrame()
res_df["ALE"] = [i.name for i in ale.traverse()]
res_df["Node"] = [map[i.name] for i in ale.traverse()]
res_df["br_length"] = res_df.apply(lambda x: ale.search_nodes(name = x["ALE"])[0].dist, axis = 1)
res_df["dist_to_root"] = res_df.apply(lambda x: round(ale.search_nodes(name = x["ALE"])[0].get_distance(ale),6), axis = 1)
res_df["N_transfers_donor"] = res_df.apply(lambda x: round(full_uTs.loc[full_uTs[0] == x["ALE"]][2].sum(), 6), axis = 1)
res_df["N_transfers_recip"] = res_df.apply(lambda x: round(full_uTs.loc[full_uTs[1] == x["ALE"]][2].sum(), 6), axis = 1)
full_uTs["to_direct_descendant"] = full_uTs.apply(lambda x: is_to_strict_descendant(x[0], x[1], ale), axis = 1)
full_uTs["to_descendant"] = full_uTs.apply(lambda x: is_to_time_descendant(x[0], x[1], ale), axis = 1)
res_df["to_direct_descendant"] = res_df.apply(lambda x: round(full_uTs.loc[(full_uTs[0] == x["ALE"]) & (full_uTs.to_direct_descendant)][2].sum(), 6), axis = 1)
res_df["to_descendant"] = res_df.apply(lambda x: round(full_uTs.loc[(full_uTs[0] == x["ALE"]) & (full_uTs.to_descendant)][2].sum(), 6), axis = 1)
res_df.to_csv(sys.argv[1], sep=' ', index = False)

# GNU Ghost
