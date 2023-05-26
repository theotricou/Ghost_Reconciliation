#!/usr/bin/env python3
# by Theo Tricou

from ete3 import Tree as tr

# read a Complete tree with extnct speceis
# the corresponding extant tree
# and the ALE tree with ale style node ID

com = tr("T/CompleteTree.nwk", format = 1)
ext = tr("SAMPLE_1/SampledSpeciesTree.nwk", format = 1)
ale = tr("aletree", format = 1)

# Creates a dictionary for ale and extant tree nodes mapping
dict_node = {}

for i in ext:
    cor_node = ale.search_nodes(name=i.name)[0]
    dict_node[i.name] = cor_node.name
    while i.name != ext.name:
        i = i.up
        if not i.name in dict_node:
            cor_node = cor_node.up
            dict_node[i.name] = cor_node.name

# Dictionary with mapp information between Complete tree's nodes and
# Extant tree's nodes.
dict_map = {}

# Here, all extinct nodes that are deeper that the mca of all extant nodes,
# if any, are added to the dictionary.
cor_node = com.search_nodes(name=ext.name)[0]
if ext.name != com.name:
    while cor_node:
        dict_map[cor_node.name] = [cor_node.name, ext.name, ext.name]
        cor_node = cor_node.up
    for i in list(set(com.get_descendants()) - set(com.search_nodes(name=ext.name)[0].get_descendants())):
        if not i.name in dict_map:
            dict_map[i.name] = [i.name, ext.name, "None"]
else:
    dict_map[cor_node.name] = [cor_node.name, ext.name, ext.name]

# Then node ancestors that are more shallow than the mca of all extant nodes.
for i in ext.iter_descendants():
    cor_node = com.search_nodes(name=i.name)[0]
    dict_map[cor_node.name] = [cor_node.name, i.name, i.name]
    cor_node = cor_node.up
    while not cor_node.name in dict_node:
        if not cor_node.name in dict_map:
            dict_map[cor_node.name] = [cor_node.name, i.name, i.name]
        cor_node = cor_node.up


ext_bl  = {}

anc = com.search_nodes(name=ext.name)[0]
extroot_to_comroot = anc.get_distance(com)
ext_bl[com.name] = 0

for i in com.iter_descendants():
    # compute the length of branch that are contempory ghost to extant/sampled braches
    anc = com.search_nodes(name=i.name)[0]
    dist_to_root = anc.get_distance(com)
    distup_to_root = anc.up.get_distance(com)
    if dist_to_root > extroot_to_comroot and distup_to_root > extroot_to_comroot:
        ext_bl[i.name] = round(anc.dist, 6)
    elif dist_to_root > extroot_to_comroot and distup_to_root < extroot_to_comroot:
        ext_bl[i.name] = round((dist_to_root - distup_to_root) - (extroot_to_comroot - distup_to_root), 6)
    else:
        ext_bl[i.name] = 0
    # Then all remaining extinct nodes
    if not i.name in dict_map:
        node = i
        while not node.name in dict_map:
            node = node.up
        dict_map[i.name] = [i.name, dict_map[node.name][1], "None"]

fout=open("Nodes_mapping.txt","w")
fout.write(" ".join(["Node", "Cor_Donor", "Cor_Recip", "ALE_ID"])+"\n")
for i in sorted(dict_map.keys()):
    fout.write(" ".join([" ".join(dict_map[i]), dict_node[dict_map[i][1]]])+"\n")

fout.close()


for i in dict_node:
    dict_node[i] = [dict_node[i], 0, 0]

for i in dict_map:
    dict_node[dict_map[i][1]][1] += ext_bl[dict_map[i][0]]
    dict_node[dict_map[i][1]][2] += 1




# GNU Ghost
